"""Import Sanika's Phase 1 run (Config A baseline) into the Phase 2 pipeline's CSV format.

Produces results/A/crops.csv and results/A/raw_predictions.csv, matching the
interfaces agreed in the 1.4 report, so the shared stages can score them.
"""
import argparse
from pathlib import Path

import pandas as pd


def make_crop_id(paper_id, page, crop_number):
    """Shared crop ID format, e.g. paper 1, page 3, crop 2 -> 1_p3_002."""
    return f"{paper_id}_p{page}_{int(crop_number):03d}"


def import_baseline(source_dir, output_dir):
    source_dir, output_dir = Path(source_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    preds = pd.read_csv(source_dir / "course_predictions.csv")
    preds["paper_id"] = preds["paperID"].astype(str)
    preds["crop_id"] = [
        make_crop_id(p, page, n)
        for p, page, n in zip(preds["paper_id"], preds["page"], preds["crop_number"])
    ]

    # crops.csv: in Config A the whole page is the "figure", so figure_id = paper + page
    crops = pd.DataFrame({
        "crop_id": preds["crop_id"],
        "figure_id": [f"{p}_p{page}" for p, page in zip(preds["paper_id"], preds["page"])],
        "paper_id": preds["paper_id"],
        "page": preds["page"],
        "x1": preds["x1"], "y1": preds["y1"], "x2": preds["x2"], "y2": preds["y2"],
        "confidence": preds["detector_confidence"],
        "crop_path": preds["image_path"],
    })
    crops.to_csv(output_dir / "crops.csv", index=False)

    # raw_predictions.csv: cleaned SMILES from MolGlyph, one row per crop
    raw = pd.DataFrame({
        "crop_id": preds["crop_id"],
        "paper_id": preds["paper_id"],
        "raw_smiles": preds["predicted_smiles"],
        "inference_time": preds["recognition_seconds"],
    })
    raw.to_csv(output_dir / "raw_predictions.csv", index=False)

    print(f"Imported {len(crops)} crops from {preds['paper_id'].nunique()} papers into {output_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Import Sanika's Phase 1 baseline results.")
    parser.add_argument("--source", default="../1.2_Individual_Benchmarking_Code/course_results",
                        help="Folder containing Sanika's course_predictions.csv")
    parser.add_argument("--output", default="results/A", help="Where to write Config A CSVs")
    args = parser.parse_args()
    import_baseline(args.source, args.output)