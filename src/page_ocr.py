"""Full-page OCR, for text that lives inside images (like the labels in a synthesis scheme).

Saves every text piece found, with its position in 200 DPI pixels (the same units as the crop
boxes), so later stages can match labels to the structures near them.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
import yaml
from PIL import Image, ImageOps


def ocr_page(results_dir, pdf_dir, paper_id, page, dpi=400, use_gpu=True):
    import easyocr
    pix = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")[page - 1].get_pixmap(dpi=dpi)
    img = ImageOps.autocontrast(Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("L"))
    reader = easyocr.Reader(["en"], gpu=use_gpu, verbose=False)
    s = 200 / dpi
    rows = []
    for box, text, conf in reader.readtext(np.array(img)):
        xs, ys = [p[0] for p in box], [p[1] for p in box]
        rows.append({"text": text, "conf": round(float(conf), 3),
                     "x1": round(min(xs) * s), "y1": round(min(ys) * s),
                     "x2": round(max(xs) * s), "y2": round(max(ys) * s)})
    df = pd.DataFrame(rows)
    df.to_csv(Path(results_dir) / f"page_ocr_paper{paper_id}_page{page}.csv", index=False)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full-page OCR with positions.")
    parser.add_argument("--results", default="results/A300")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", required=True)
    parser.add_argument("--page", type=int, required=True)
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    df = ocr_page(args.results, cfg["data"]["pdf_dir"], args.paper, args.page, use_gpu=not args.cpu)
    looks_like_label = df.text.str.contains("=") | df.text.str.fullmatch(r"1[129]|2[01]|39|4[01]")
    print(df[looks_like_label].sort_values(["y1", "x1"]).to_string(index=False))