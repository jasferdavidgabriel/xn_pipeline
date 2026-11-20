"""Build a visual review sheet of crops for one paper and status, plus a CSV to record verdicts."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image


def build_review(results_dir, crop_root, paper_id, status, n, cols=4):
    results_dir, crop_root = Path(results_dir), Path(crop_root)
    crops = pd.read_csv(results_dir / "crops.csv", dtype={"paper_id": str})
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})

    # Join each prediction to its crop image, then keep the requested paper and status
    review = validated.merge(crops[["crop_id", "crop_path", "confidence"]], on="crop_id")
    review = review[(review["paper_id"] == str(paper_id)) & (review["status"] == status)].head(n)
    if review.empty:
        print("No crops match that paper and status.")
        return

    rows = -(-len(review) // cols)  # ceiling division
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3.6))
    for ax in axes.flat:
        ax.axis("off")
    for ax, (_, r) in zip(axes.flat, review.iterrows()):
        ax.imshow(Image.open(crop_root / r["crop_path"]))
        smiles = str(r["canonical_smiles"])
        smiles = smiles if len(smiles) <= 38 else smiles[:35] + "..."
        ax.set_title(f"{r['crop_id']}  (conf {r['confidence']:.2f})\n{smiles}", fontsize=8)

    name = f"review_paper{paper_id}_{status}"
    plt.tight_layout()
    plt.savefig(results_dir / f"{name}.png", dpi=150)

    # A sheet for recording what you see in each crop
    sheet = review[["crop_id", "canonical_smiles"]].copy()
    sheet["verdict"] = ""  # scaffold_with_R | abbreviation | not_a_molecule | partial_or_merged
    sheet["notes"] = ""
    sheet.to_csv(results_dir / f"{name}.csv", index=False)
    print(f"Saved {name}.png and {name}.csv with {len(review)} crops")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visual review of crops by paper and status.")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--crop-root", default="../1.2_Individual_Benchmarking_Code",
                        help="Folder the crop_path values are relative to")
    parser.add_argument("--paper", default="2")
    parser.add_argument("--status", default="wildcard")
    parser.add_argument("--n", type=int, default=12)
    args = parser.parse_args()
    build_review(args.results, args.crop_root, args.paper, args.status, args.n)