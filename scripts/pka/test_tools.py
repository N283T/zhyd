import json
import tempfile
import unittest
from pathlib import Path

import compare_engines
import compare_propka
import evaluate
import pkadr


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


if __name__ == "__main__":
    unittest.main()
