"""Rebuild crop images from the PDFs, using the box coordinates saved in crops.csv.

Boxes are stored in pixels at box_dpi (200 for the Phase 1 baseline). Rendering at a
different dpi scales every box and margin to match, so the same regions are cropped
at a different resolution.
"""
import argparse
from pathlib import Path

import pandas as pd
import pymupdf
import yaml
from PIL import Image


def rebuild_crops(results_dir, pdf_dir, dpi=200, margin=10, box_dpi=200):
    results_dir = Path(results_dir)
    crop_dir = results_dir / "crops"
    crop_dir.mkdir(parents=True, exist_ok=True)
    crops = pd.read_csv(results_dir / "crops.csv", dtype={"paper_id": str})
    scale = dpi / box_dpi
    pad = int(round(margin * scale))

    for paper_id, paper_rows in crops.groupby("paper_id"):
        doc = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")
        for page_no, page_rows in paper_rows.groupby("page"):
            pix = doc[int(page_no) - 1].get_pixmap(dpi=dpi)   # pages are 1 based in the CSV
            page_img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
            for idx, r in page_rows.iterrows():
                box = (
                    max(0, int(r["x1"] * scale) - pad), max(0, int(r["y1"] * scale) - pad),
                    min(pix.width, int(r["x2"] * scale) + pad), min(pix.height, int(r["y2"] * scale) + pad),
                )
                crop_path = crop_dir / f"{r['crop_id']}.png"
                page_img.crop(box).save(crop_path)
                crops.loc[idx, "crop_path"] = str(crop_path)
        doc.close()

    crops.to_csv(results_dir / "crops.csv", index=False)
    print(f"Rebuilt {len(crops)} crops at {dpi} DPI into {crop_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rebuild crop images from PDFs and saved boxes.")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--dpi", type=int, default=200, help="Render resolution for the crops")
    args = parser.parse_args()
    config = yaml.safe_load(open(args.config))
    rebuild_crops(args.results, config["data"]["pdf_dir"], dpi=args.dpi)