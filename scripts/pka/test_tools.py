import json
import tempfile
import unittest
from pathlib import Path

import calibrate_empirical
import compare_engines
import compare_propka
import evaluate
import external_set
import merge_sets
import pkadr
import score_literature


def record(index, pdb, res="ASP", pka="4.0", name="P", similar="", warning=""):
    return {
        "Index": str(index),
        "PDB": pdb,
        "Chain": "A",
        "ResID_in_PDB": "10",
        "ResName": res,
        "Expt_pKa": pka,
        "pKa_Classification": "Main",
        "Mut_Pos": "",
        "Protein_Name": name,
        "Sequence_Identity_>_30%": similar,
        "Alternative_PDBs": "",
        "Warning": warning,
    }


class ToolTests(unittest.TestCase):
    def test_cluster_split_never_separates_cluster(self):
        self.assertEqual(
            evaluate.cluster_split("cluster-1", "seed", 0.8),
            evaluate.cluster_split("cluster-1", "seed", 0.8),
        )

    def test_metrics(self):
        rows = [
            {
                "experimental_pka": "4",
                "predicted_pka": "5",
                "model_pka": "4",
                "comp_id": "ASP",
                "burial_fraction": "0",
                "cluster_id": "x",
                "split": "train",
            }
        ]
        self.assertEqual(evaluate.metrics(rows)["mae"], 1.0)
        self.assertEqual(evaluate.metrics(rows)["rmse"], 1.0)

    def test_state_flips_and_shift_strata_use_reference_pka(self):
        def row(site, experimental, predicted):
            return {
                "comp_id": site,
                "experimental_pka": str(experimental),
                "predicted_pka": str(predicted),
                "model_pka": "0",  # engine model pKa must not matter
                "burial_fraction": "0",
            }

        rows = [
            row("HIS", 6.6, 6.0),  # no flip; state right (both deprotonated)
            row("HIS", 7.5, 7.2),  # flip vs reference 6.5; predicted right
            row("CYS", 6.0, 8.0),  # flip vs reference 9.0; predicted wrong
        ]
        m = evaluate.metrics(rows)
        self.assertEqual(m["state_flips"], 2)
        self.assertEqual(m["state_flips_correct"], 1)
        self.assertAlmostEqual(float(m["state_accuracy"] or 0), 2 / 3)
        self.assertAlmostEqual(
            float(m["reference_null_mae"] or 0), (0.1 + 1.0 + 3.0) / 3
        )
        self.assertEqual([evaluate.shift_bin(r) for r in rows], ["<0.5", "1-2", ">=2"])

    def test_propka_summary_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.pka"
            path.write_text(
                "SUMMARY OF THIS PREDICTION\n ASP  42 A     4.25       3.80\n"
            )
            self.assertEqual(compare_propka.parse_propka(path)[("ASP", 42, "A")], 4.25)

    def test_propka_summary_parser_handles_large_numbers_and_arg(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.pka"
            path.write_text(
                "SUMMARY OF THIS PREDICTION\n"
                "   ASP1014 D     3.10       3.80\n"
                "   ARG  12 A    12.40      12.50\n"
            )
            values = compare_propka.parse_propka(path)
            self.assertEqual(values[("ASP", 1014, "D")], 3.10)
            self.assertEqual(values[("ARG", 12, "A")], 12.40)


class PkadrTests(unittest.TestCase):
    def test_similar_entries_share_a_cluster(self):
        rows, _ = pkadr.normalize(
            [
                record(1, "1AAA", name="Alpha"),
                record(2, "2BBB", name="Beta", similar="1AAA.A"),
                record(3, "3CCC", name="Gamma"),
            ]
        )
        clusters = {row["pdb_id"]: row["cluster_id"] for row in rows}
        self.assertEqual(clusters["1AAA"], clusters["2BBB"])
        self.assertNotEqual(clusters["1AAA"], clusters["3CCC"])

    def test_censored_values_and_unsupported_residues_are_excluded(self):
        rows, summary = pkadr.normalize(
            [
                record(1, "1AAA", pka="<2.0"),
                record(2, "1AAA", res="ARG"),
                record(3, "1AAA", res="N-term"),
            ]
        )
        self.assertEqual([row["site"] for row in rows], ["N_TERM"])
        self.assertEqual(summary["excluded"]["censored_or_nonnumeric_pka"], 1)
        self.assertEqual(summary["excluded"]["unsupported_residue"], 1)

    def test_snase_is_a_stress_holdout(self):
        rows, _ = pkadr.normalize([record(1, "3BDC", name="SNase")])
        row = {**rows[0], "model_pka": "4", "predicted_pka": "4"}
        self.assertEqual(evaluate.assign_split(row, "seed", 0.8), "stress:snase")

    def test_external_set_drops_pkadr_structures_and_homologs(self):
        pkadr_records = [{"PDB": "1AAA", "Alternative_PDBs": "1BBB"}]
        clusters = {"1AAA": {0}, "1BBB": {1}, "2CCC": {1, 5}, "3DDD": {7}}
        pkad3 = [
            {
                "PDB_ID": "1BBB",
                "Chain": "A",
                "Res_Name": "D",
                "resid": "5",
                "Expt_pKa": "3.0",
                "note": "/",
            },
            {
                "PDB_ID": "2CCC",
                "Chain": "A",
                "Res_Name": "E",
                "resid": "6",
                "Expt_pKa": "4.0",
                "note": "/",
            },
            {
                "PDB_ID": "9ZZZ",
                "Chain": "A",
                "Res_Name": "K",
                "resid": "7",
                "Expt_pKa": "10.0",
                "note": "/",
            },
            {
                "PDB_ID": "3DDD",
                "Chain": "A",
                "Res_Name": "H",
                "resid": "8",
                "Expt_pKa": "6.0",
                "note": "/",
            },
            {
                "PDB_ID": "3DDD",
                "Chain": "A",
                "Res_Name": "H",
                "resid": "8",
                "Expt_pKa": "7.0",
                "note": "/",
            },
            {
                "PDB_ID": "3DDD",
                "Chain": "A",
                "Res_Name": "D",
                "resid": "9",
                "Expt_pKa": "2.0",
                "note": "T26H",
            },
        ]
        blind = [
            {
                "PDB_ID": "3DDD",
                "Chain": "A",
                "Res_Name": "LYS",
                "Res_ID": "10",
                "Expt_pKa": "10.5",
            }
        ]
        rows, summary = external_set.build(pkad3, blind, pkadr_records, clusters)
        self.assertEqual(summary["pkadr_structure"], 1)
        self.assertEqual(summary["homolog_of_pkadr"], 1)
        self.assertEqual(summary["not_in_cluster_file"], 1)
        self.assertEqual(
            [(r["auth_seq_id"], r["site"], r["experimental_pka"]) for r in rows],
            [("10", "LYS", "10.500"), ("8", "HIS", "6.500")],
        )
        self.assertTrue(all(r["holdout_family"] == "external" for r in rows))
        self.assertEqual(evaluate.assign_split(rows[0], "seed", 0.8), "stress:external")

    def test_compare_engines_reports_term_deltas_and_missing_sites(self):
        def site(seq, predicted, burial):
            return {
                "id": {
                    "chain_auth": "A",
                    "auth_seq_id": seq,
                    "ins_code": None,
                    "site": "ASP",
                },
                "pka": {
                    "predicted": predicted,
                    "terms": {
                        "model": 3.8,
                        "burial": burial,
                        "hydrogen_bond": 0.0,
                        "electrostatic": 0.0,
                    },
                },
            }

        oracle = {
            "models": [{"model": 1, "sites": [site(1, 4.0, 0.2), site(2, 3.0, 0.0)]}]
        }
        candidate = {"models": [{"model": 1, "sites": [site(1, 4.05, 0.25)]}]}
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for name, report in (("o.json", oracle), ("c.json", candidate)):
                path = Path(tmp) / name
                path.write_text(json.dumps(report))
                paths.append(path)
            rows, only_oracle, only_candidate = compare_engines.compare(
                compare_engines.load_sites(paths[0]),
                compare_engines.load_sites(paths[1]),
            )
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]["delta"], 0.05)
        self.assertAlmostEqual(rows[0]["delta_burial"], 0.05)
        self.assertEqual(only_oracle, [(1, "A", 2, "", "ASP")])
        self.assertEqual(only_candidate, [])
        ok, _ = compare_engines.summarize(rows, 0.02, 0.1, 0.99)
        self.assertFalse(ok)

    def test_predict_appends_engine_and_extra_args(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            structures = root / "structures"
            structures.mkdir()
            (structures / "1ABC.cif").write_text("data_1ABC\n")
            argv_log = root / "argv.json"
            # Stand-in for zhyd: record argv, write a one-site report.
            fake = root / "fake_zhyd.py"
            fake.write_text(
                "import json, sys\n"
                f"open({str(argv_log)!r}, 'w').write(json.dumps(sys.argv[1:]))\n"
                "out = sys.argv[sys.argv.index('-o') + 1]\n"
                "site = {'id': {'chain_auth': 'A', 'auth_seq_id': 10, 'site': 'ASP'},"
                " 'pka': {'predicted': 4.1, 'burial_fraction': 0.0, 'terms': {"
                "'model': 3.8, 'burial': 0.3, 'hydrogen_bond': 0.0,"
                " 'electrostatic': 0.0, 'coupling': 0.0}}, 'confidence': 'high'}\n"
                "json.dump({'models': [{'sites': [site]}]}, open(out, 'w'))\n"
            )
            fake.chmod(0o755)
            launcher = root / "zhyd"
            launcher.write_text(f'#!/bin/sh\nexec python3 {fake} "$@"\n')
            launcher.chmod(0o755)
            rows, _ = pkadr.normalize([record(1, "1ABC")])
            joined, status = pkadr.predict(
                rows,
                structures,
                launcher,
                root / "reports",
                engine="empirical",
                extra_args=["--empirical-mode", "zhyd"],
            )
            argv = json.loads(argv_log.read_text())
        self.assertEqual(status["matched"], 1)
        self.assertEqual(joined[0]["predicted_pka"], "4.1")
        self.assertEqual(
            argv[-4:], ["--engine", "empirical", "--empirical-mode", "zhyd"]
        )

    def test_nested_cv_reports_held_out_fold_errors(self):
        rows = [
            {"site": "CYS", "cluster_id": c, "experimental_pka": "8.0"}
            for c in "abcdefghijklmnopqrst"
        ]
        default = (9.0, 1.0)
        table = {
            (9.0, 1.0): [9.0] * 20,
            (8.0, 1.0): [8.0] * 20,
        }
        errors = {p: calibrate_empirical.abs_errors(v, rows) for p, v in table.items()}
        nested = calibrate_empirical.nested_cv(
            errors, rows, "CYS", [8.0, 9.0], [1.0], default, 2, "seed"
        )
        self.assertEqual(nested["default_mae"], 1.0)
        self.assertEqual(nested["nested_cv_mae"], 0.0)
        self.assertTrue(all(c["point"] == [8.0, 1.0] for c in nested["chosen"]))

    def test_merge_joins_clusters_through_shared_sequence_clusters(self):
        base = {
            "chain": "A",
            "auth_seq_id": "1",
            "site": "ASP",
            "experimental_pka": "3",
        }
        a = [
            {
                **base,
                "entry_id": "1",
                "pdb_id": "1AAA",
                "cluster_id": "x",
                "holdout_family": "",
            },
            {
                **base,
                "entry_id": "2",
                "pdb_id": "1SNA",
                "cluster_id": "s",
                "holdout_family": "snase",
            },
        ]
        b = [
            {
                **base,
                "entry_id": "1",
                "pdb_id": "2BBB",
                "cluster_id": "y",
                "holdout_family": "external",
            }
        ]
        clusters = {"1AAA": {3}, "2BBB": {3, 4}, "1SNA": {9}}
        rows, summary = merge_sets.merge([("r", a), ("e", b)], clusters, {"snase"})
        self.assertEqual([r["entry_id"] for r in rows], ["r:1", "r:2", "e:1"])
        self.assertEqual(rows[0]["cluster_id"], rows[2]["cluster_id"])
        self.assertNotEqual(rows[0]["cluster_id"], rows[1]["cluster_id"])
        self.assertEqual([r["holdout_family"] for r in rows], ["", "snase", ""])
        self.assertEqual(summary["fitting_rows"], 2)

    def test_empirical_calibration_selects_by_cluster_cv(self):
        # Two CYS rows in clusters a and b and one ASP row in cluster c.
        rows = [
            {"site": "CYS", "cluster_id": "a", "experimental_pka": "8.0"},
            {"site": "CYS", "cluster_id": "b", "experimental_pka": "8.2"},
            {"site": "ASP", "cluster_id": "c", "experimental_pka": "4.0"},
        ]
        default = (9.0, 1.0)
        models = [8.0, 9.0]
        scales = [0.5, 1.0]
        # Predictions per (model pKa, scale): lowering the model pKa helps both
        # CYS rows; the scale only helps row a and hurts row b.
        table = {
            (9.0, 1.0): [9.0, 9.2, 4.0],
            (8.0, 1.0): [8.0, 8.2, 4.0],
            (9.0, 0.5): [8.5, 9.6, 4.0],
            (8.0, 0.5): [7.9, 8.6, 4.0],
        }
        errors = {p: calibrate_empirical.abs_errors(v, rows) for p, v in table.items()}
        self.assertEqual(
            calibrate_empirical.fit(errors, list(table), [0, 1, 2], default),
            (8.0, 1.0),
        )
        # Leaving cluster a out, the fit on b and c alone is still (8.0, 1.0).
        cv, chosen = calibrate_empirical.loco_cv(
            errors, list(table), ["a", "b", "c"], default
        )
        self.assertEqual(chosen["a"], (8.0, 1.0))
        self.assertEqual(cv, [0.0, 0.0, 0.0])
        result = calibrate_empirical.select(
            errors, rows, "CYS", models, scales, default
        )
        self.assertEqual(result["selected"], "model_pka")
        self.assertEqual(result["variants"]["default"]["cv_mae"], round(2.0 / 3.0, 4))
        self.assertEqual(result["variants"]["model_pka"]["cv_mae"], 0.0)
        # Ties prefer the point closest to the default.
        flat = {p: [0.0, 0.0, 0.0] for p in table}
        self.assertEqual(
            calibrate_empirical.fit(flat, list(flat), [0, 1, 2], default), default
        )
        overrides = calibrate_empirical.overrides(
            "v", "CYS", (8.0, 0.5), {"model_pka": {"TYR": 9.5}}
        )
        # Every kind starts from the parity table, whatever the binary's
        # built-in recalibrated table is; fixed kinds and the point win.
        self.assertEqual(overrides["model_pka"]["CYS"], 8.0)
        self.assertEqual(overrides["model_pka"]["TYR"], 9.5)
        self.assertEqual(overrides["model_pka"]["HIS"], 6.5)
        self.assertEqual(overrides["desolvation_scale"]["CYS"], 0.5)
        self.assertEqual(overrides["desolvation_scale"]["HIS"], 1.0)
        self.assertEqual(
            set(overrides["model_pka"]), set(calibrate_empirical.DEFAULT_MODEL_PKA)
        )

    def test_literature_scoring_handles_upper_bounds(self):
        def row(pdb, seq, value, qualifier="point"):
            return {
                "entry_id": f"{pdb}-{seq}",
                "pdb_id": pdb,
                "chain": "A",
                "auth_seq_id": str(seq),
                "site": "HIS",
                "experimental_pka": value,
                "value_qualifier": qualifier,
            }

        self.assertIs(
            score_literature.experimental_state(row("X", 1, "6.5"), 7.0), False
        )
        self.assertIs(
            score_literature.experimental_state(row("X", 1, "3.5", "upper_bound"), 7.0),
            False,
        )
        self.assertIsNone(
            score_literature.experimental_state(row("X", 1, "8.0", "upper_bound"), 7.0)
        )

        def site(seq, pka, confidence="high"):
            return {
                "id": {"chain_auth": "A", "auth_seq_id": seq, "site": "HIS"},
                "pka": {"predicted": pka},
                "confidence": confidence,
            }

        with tempfile.TemporaryDirectory() as tmp:
            report = {"models": [{"sites": [site(1, 6.0), site(2, 7.5, "low")]}]}
            (Path(tmp) / "X.json").write_text(json.dumps(report))
            rows = [
                row("X", 1, "6.5"),
                row("X", 2, "3.5", "upper_bound"),
                row("X", 3, "6.0"),
            ]
            result = score_literature.score(rows, Path(tmp), 7.0)
        self.assertEqual(result["points"], 1)
        self.assertAlmostEqual(result["mae"], 0.5)
        self.assertEqual((result["bounds_met"], result["bounds"]), (0, 1))
        self.assertEqual((result["state_agree"], result["states"]), (1, 2))
        self.assertEqual(result["low_confidence"], ["X-2"])
        self.assertEqual(result["missing"], ["X-3"])


if __name__ == "__main__":
    unittest.main()
