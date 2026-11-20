"""Read short text in a table image (X atoms, dashes) with targeted OCR on known cell positions.

Full page OCR skips isolated single characters (C, N, digits). Cropping each small region,
enlarging it, and restricting the allowed characters makes them readable.
"""
import argparse
from pathlib import Path

import easyocr
import numpy as np
import pandas as pd
import pymupdf
import yaml
from PIL import Image

from src.assemble import COL_GAP_PX, ROW_GAP_PX, cluster

UPSCALE = 3            # enlarge small crops so OCR can see single characters
X_BAND_WIDTH_PX = 90   # how far right of the R4 column to look for the X letter
ROW_HALF_HEIGHT_PX = 22


def ocr_crop(reader, image, allowlist):
    """Recognize the whole crop as one piece of text (skips text detection, which drops tiny glyphs)."""
    big = image.resize((image.width * UPSCALE, image.height * UPSCALE), Image.LANCZOS)
    results = reader.recognize(np.array(big.convert("L")), allowlist=allowlist)
    if not results:
        return "", 0.0
    _, text, conf = results[0]
    return text.strip(), round(float(conf), 3)


def read_table_text(results_dir, pdf_dir, paper_id, page, use_gpu=True, dpi=200):
    results_dir = Path(results_dir)
    crops = pd.read_csv(results_dir / "crops.csv", dtype={"paper_id": str})
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    m = crops.merge(validated[["crop_id", "status"]], on="crop_id")
    frags = m[(m.paper_id == str(paper_id)) & (m.page == page) & (m.status == "fragment")].copy()
    frags["xc"], frags["yc"] = (frags.x1 + frags.x2) / 2, (frags.y1 + frags.y2) / 2
    frags["row"] = cluster(list(frags.yc), ROW_GAP_PX)
    frags["col"] = cluster(list(frags.xc), COL_GAP_PX)

    last_col = frags.col.max()
    r4_right = frags[frags.col == last_col].x2.max()   # right edge of the R4 column
    row_y = frags.groupby("row").yc.mean()

    pix = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")[page - 1].get_pixmap(dpi=dpi)
    page_img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    reader = easyocr.Reader(["en"], gpu=use_gpu, verbose=False)

    # X atom per row: the letter just right of the R4 column
    rows = []
    for row, y in row_y.items():
        band = page_img.crop((int(r4_right + 5), int(y - ROW_HALF_HEIGHT_PX),
                              int(r4_right + 5 + X_BAND_WIDTH_PX), int(y + ROW_HALF_HEIGHT_PX)))
        text, conf = ocr_crop(reader, band, allowlist="CN")
        rows.append({"row": row, "compound_inferred": row + 1,
                     "x_atom": text[:1] if text[:1] in ("C", "N") else "",
                     "x_ocr_text": text, "x_conf": conf})
    rows = pd.DataFrame(rows)
    rows.to_csv(results_dir / f"table_text_paper{paper_id}_page{page}.csv", index=False)

    # Empty cells: check which ones are a dash (no substituent)
    cells_csv = results_dir / f"missing_cells_paper{paper_id}_page{page}.csv"
    if cells_csv.exists():
        cells = pd.read_csv(cells_csv)
        reads = [ocr_crop(reader, Image.open(p), allowlist="-HFClBrO") for p in cells.crop_path]
        cells["cell_ocr_text"] = [t for t, _ in reads]
        cells["cell_ocr_conf"] = [c for _, c in reads]
        cells["is_dash"] = cells.cell_ocr_text.isin(["-", "_"])
        cells.to_csv(cells_csv, index=False)
        print(cells[["cell_id", "cell_ocr_text", "cell_ocr_conf", "is_dash"]].to_string(index=False))

    print(rows.to_string(index=False))
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Targeted OCR of table text (X atoms, dashes).")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", default="2")
    parser.add_argument("--page", type=int, default=5)
    parser.add_argument("--cpu", action="store_true", help="Run OCR on CPU instead of GPU")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    read_table_text(args.results, cfg["data"]["pdf_dir"], args.paper, args.page, use_gpu=not args.cpu)