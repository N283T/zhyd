#!/usr/bin/env python3
"""Evaluate canonical pKa CSV files with leakage-resistant cluster splits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

REQUIRED = {
    "comp_id",
    "experimental_pka",
    "predicted_pka",
    "model_pka",
    "burial_fraction",
    "cluster_id",
}


def cluster_split(cluster_id: str, seed: str, train_fraction: float) -> str:
    digest = hashlib.sha256(f"{seed}\0{cluster_id}".encode()).digest()
    fraction = int.from_bytes(digest[:8], "big") / 2**64
    return "train" if fraction < train_fraction else "test"


def metrics(rows: list[dict[str, str]]) -> dict[str, float | int | None]:
    if not rows:
        return {"n": 0, "mae": None, "rmse": None, "shift_direction_accuracy": None}
    errors = [
        float(row["predicted_pka"]) - float(row["experimental_pka"]) for row in rows
    ]
    directions: list[bool] = []
    for row in rows:
        model = float(row["model_pka"])
        experimental_shift = float(row["experimental_pka"]) - model
        predicted_shift = float(row["predicted_pka"]) - model
        if abs(experimental_shift) >= 0.1:
            directions.append((experimental_shift > 0) == (predicted_shift > 0))
    return {
        "n": len(rows),
        "mae": sum(abs(value) for value in errors) / len(errors),
        "rmse": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "shift_direction_accuracy": (sum(directions) / len(directions))
        if directions
        else None,
    }


def evaluate(rows: list[dict[str, str]]) -> dict[str, object]:
    report: dict[str, object] = {"overall": metrics(rows)}
    for dimension, key_fn in (
        ("split", lambda row: row["split"]),
        ("residue", lambda row: row["comp_id"].upper()),
        (
            "exposure",
            lambda row: "buried" if float(row["burial_fraction"]) >= 0.5 else "surface",
        ),
    ):
        groups: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            groups[key_fn(row)].append(row)
        report[dimension] = {
            key: metrics(value) for key, value in sorted(groups.items())
        }
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", default="zhyd-pka-split-v1")
    parser.add_argument("--train-fraction", type=float, default=0.8)
    args = parser.parse_args()
    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(f"missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    for row in rows:
        assigned = cluster_split(row["cluster_id"], args.seed, args.train_fraction)
        if row.get("split") and row["split"] != assigned:
            raise SystemExit(
                f"row has split {row['split']!r}, expected cluster-derived {assigned!r}"
            )
        row["split"] = assigned
    payload = {
        "schema_version": "zhyd.pka-evaluation.v1",
        "split": {
            "unit": "protein_or_sequence_cluster",
            "seed": args.seed,
            "train_fraction": args.train_fraction,
        },
        "metrics": evaluate(rows),
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
