"""Score final predictions: complete molecules from recognition plus every assembled table."""
import argparse
from pathlib import Path

import pandas as pd
import yaml

from src.evaluate import evaluate_all
from src.load_data import load_dataset


def final_predictions(results_dir):
    """Complete molecules from validation, plus all assembled_*.csv outputs in results_dir."""
    results_dir = Path(results_dir)
    v = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    parts = [v[v.status == "valid"][["paper_id", "canonical_smiles"]]
             .rename(columns={"canonical_smiles": "predicted_smiles"})]
    for path in sorted(results_dir.glob("assembled_*.csv")):
        a = pd.read_csv(path, dtype={"paper_id": str}).dropna(subset=["assembled_smiles"])
        parts.append(a[["paper_id", "assembled_smiles"]].rename(columns={"assembled_smiles": "predicted_smiles"}))
        print(f"Including {len(a)} assembled molecules from {path.name}")
    return pd.concat(parts, ignore_index=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Score the pipeline's final predictions.")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config)); d = cfg["data"]
    gt, _ = load_dataset(cfg)
    preds = final_predictions(args.results)
    _, summary, per_paper = evaluate_all(gt, preds, d["id_column"], d["name_column"], d["smiles_column"])
    per_paper.to_csv(Path(args.results) / "final_per_paper.csv")
    print(per_paper)
    print(summary)