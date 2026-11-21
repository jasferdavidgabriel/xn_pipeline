"""Read short text cells in an image table by targeted OCR, one cell at a time.

Rows come from crops with a given status (for example, the R-group fragments, one per
row). Columns are x ranges in 200 DPI pixels, each with the characters it may contain.
Run with --preview first (no OCR, runs anywhere) to check the cells are well placed.
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
import yaml

from src.pieces import load_crops
from PIL import Image, ImageOps

UPSCALE = 2           # cells are already rendered at high DPI, so less enlarging is needed
BOX_DPI = 200         # crop boxes and column positions are in 200 DPI pixels


def ink_height_ratio(image, threshold=128, side_trim=0.15):
    """Height of the dark pixels in the cell's center, as a fraction of the cell height."""
    w, h = image.size
    center = image.crop((int(w * side_trim), 0, int(w * (1 - side_trim)), h))
    ys = np.where((np.array(center.convert("L")) < threshold).any(axis=1))[0]
    return 0.0 if len(ys) == 0 else (ys.max() - ys.min() + 1) / h


def read_columns(results_dir, pdf_dir, paper_id, page, columns, row_status="fragment",
                 half_height=25, preview=False, use_gpu=True, ocr_dpi=400):
    """columns: {name: [x_from, x_to, allowlist]} in 200 DPI pixels. Returns one row per table row.

    Cells are rendered at ocr_dpi and contrast-stretched, because thin serif text at
    200 DPI renders as light gray, which OCR misreads.
    """
    crops = load_crops(results_dir, paper_id, page)
    rows = crops[crops.status == row_status].sort_values("yc")
    pix = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")[page - 1].get_pixmap(dpi=ocr_dpi)
    page_img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    s = ocr_dpi / BOX_DPI
    
    cells = []   # (row crop_id, column name, image)
    for _, r in rows.iterrows():
        for name, (x0, x1, _) in columns.items():
            box = (int(x0 * s), int((r.yc - half_height) * s), int(x1 * s), int((r.yc + half_height) * s))
            cells.append((r.crop_id, name, ImageOps.autocontrast(page_img.crop(box).convert("L"))))
            
    out_stem = Path(results_dir) / f"table_cells_paper{paper_id}_page{page}"
    if preview:
        n_cols = len(columns)
        fig, axes = plt.subplots(len(rows), n_cols, figsize=(n_cols * 2.2, len(rows) * 1.0))
        for ax, (crop_id, name, img) in zip(np.array(axes).flat, cells):
            ax.imshow(img); ax.set_title(f"{crop_id} {name}", fontsize=7); ax.axis("off")
        plt.tight_layout()
        plt.savefig(f"{out_stem}_preview.png", dpi=120)
        print(f"Preview saved to {out_stem}_preview.png")
        return None

    import easyocr   # only needed for real OCR runs (Colab)
    reader = easyocr.Reader(["en"], gpu=use_gpu, verbose=False)
    records = {}
    for crop_id, name, img in cells:
        big = img.resize((img.width * UPSCALE, img.height * UPSCALE), Image.LANCZOS)
        result = reader.recognize(np.array(big.convert("L")), allowlist=columns[name][2])
        text, conf = (result[0][1].strip(), round(float(result[0][2]), 3)) if result else ("", 0.0)
        rec = records.setdefault(crop_id, {"crop_id": crop_id})
        rec[name], rec[f"{name}_conf"], rec[f"{name}_ink"] = text, conf, round(ink_height_ratio(img), 3)
    df = pd.DataFrame(records.values())
    df.to_csv(f"{out_stem}.csv", index=False)
    print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Targeted OCR of text columns in an image table.")
    parser.add_argument("--results", default="results/A300")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", required=True)
    parser.add_argument("--page", type=int, required=True)
    parser.add_argument("--columns", required=True, help='JSON: {"name": [x_from, x_to, "allowlist"], ...}')
    parser.add_argument("--preview", action="store_true", help="Save cell images only, no OCR")
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    read_columns(args.results, cfg["data"]["pdf_dir"], args.paper, args.page, json.loads(args.columns),
                 preview=args.preview, use_gpu=not args.cpu)