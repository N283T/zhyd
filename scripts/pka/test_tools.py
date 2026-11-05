import tempfile
import unittest
from pathlib import Path

import compare_propka
import evaluate


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


if __name__ == "__main__":
    unittest.main()
