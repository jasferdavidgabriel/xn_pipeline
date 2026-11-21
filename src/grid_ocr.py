"""OCR the text cells of a table grid found by table_grid.py, merged cells included.

Each merged cell is read once over its whole area, and its text goes to every row it covers.
Cells are rendered at 400 DPI and contrast-stretched, since thin table text renders as light
gray at 200 DPI.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
import yaml
from PIL import Image, ImageOps

BOX_DPI = 200


def read_grid(results_dir, pdf_dir, paper_id, page, columns, ocr_dpi=400, use_gpu=True):
    """columns: {column name: allowed characters}. Returns one row per table band."""
    import easyocr
    results_dir = Path(results_dir)
    grid = json.loads((results_dir / f"grid_paper{paper_id}_page{page}.json").read_text())
    bands = {b["band"]: b for b in grid["bands"]}
    names = list(grid["columns"])
    x_ranges = dict(zip(names, zip(grid["v_lines"], grid["v_lines"][1:])))

    pix = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")[page - 1].get_pixmap(dpi=ocr_dpi)
    page_img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("L")
    s = ocr_dpi / BOX_DPI
    reader = easyocr.Reader(["en"], gpu=use_gpu, verbose=False)

    rows = {}
    for name, allowlist in columns.items():
        x_from, x_to = x_ranges[name]
        for group in grid["columns"][name]:
            top, bottom = bands[group[0]]["y_top"], bands[group[-1]]["y_bottom"]
            box = (int((x_from + 4) * s), int((top + 3) * s), int((x_to - 4) * s), int((bottom - 3) * s))
            img = np.array(ImageOps.autocontrast(page_img.crop(box)))
            found = reader.readtext(img, allowlist=allowlist, detail=1)
            found.sort(key=lambda r: (r[0][0][1], r[0][0][0]))       # top to bottom, then left to right
            if not found:                                              # tiny single tokens can escape detection
                found = reader.recognize(img, allowlist=allowlist)
            text = "".join(t for _, t, _ in found).strip()
            conf = round(min((float(c) for *_, c in found), default=0.0), 3)
            for band in group:
                rec = rows.setdefault(band, {"band": band})
                rec[name], rec[f"{name}_conf"], rec[f"{name}_merged"] = text, conf, len(group) > 1

    df = pd.DataFrame(sorted(rows.values(), key=lambda r: r["band"]))
    df.to_csv(results_dir / f"grid_cells_paper{paper_id}_page{page}.csv", index=False)
    print(df.to_string(index=False))
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OCR text cells of a detected table grid.")
    parser.add_argument("--results", default="results/A300")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", required=True)
    parser.add_argument("--page", type=int, required=True)
    parser.add_argument("--columns", required=True, help='JSON: {"column name": "allowed characters"}')
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    read_grid(args.results, cfg["data"]["pdf_dir"], args.paper, args.page, json.loads(args.columns),
              use_gpu=not args.cpu)