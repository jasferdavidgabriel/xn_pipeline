"""Crop table cells where detection found nothing, so recognition can read them later."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pymupdf
import yaml
from PIL import Image

from src.assemble import COL_GAP_PX, ROW_GAP_PX, cluster

CELL_MARGIN_PX = 10  # small padding so the attachment line and label are fully inside the crop


def crop_missing_cells(results_dir, pdf_dir, paper_id, page, dpi=200):
    results_dir = Path(results_dir)
    crops = pd.read_csv(results_dir / "crops.csv", dtype={"paper_id": str})
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    m = crops.merge(validated[["crop_id", "status"]], on="crop_id")
    frags = m[(m.paper_id == str(paper_id)) & (m.page == page) & (m.status == "fragment")].copy()

    frags["xc"], frags["yc"] = (frags.x1 + frags.x2) / 2, (frags.y1 + frags.y2) / 2
    frags["w"], frags["h"] = frags.x2 - frags.x1, frags.y2 - frags.y1
    frags["row"] = cluster(list(frags.yc), ROW_GAP_PX)
    frags["col"] = cluster(list(frags.xc), COL_GAP_PX)

    col_x, row_y = frags.groupby("col").xc.mean(), frags.groupby("row").yc.mean()
    # Cell size: the larger of the column's typical crop and the table's typical crop
    cell_w = frags.groupby("col").w.median().clip(lower=frags.w.median())
    cell_h = frags.groupby("col").h.median().clip(lower=frags.h.median())
    filled = set(zip(frags.row, frags.col))

    pix = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")[page - 1].get_pixmap(dpi=dpi)
    page_img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    out_dir = results_dir / "cells"
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for row, y in row_y.items():
        for col, x in col_x.items():
            if (row, col) in filled:
                continue
            half_w, half_h = cell_w[col] / 2 + CELL_MARGIN_PX, cell_h[col] / 2 + CELL_MARGIN_PX
            box = (int(x - half_w), int(y - half_h), int(x + half_w), int(y + half_h))
            cell_id = f"{paper_id}_p{page}_r{row:02d}_c{col}"
            path = out_dir / f"{cell_id}.png"
            page_img.crop(box).save(path)
            records.append({"cell_id": cell_id, "paper_id": str(paper_id), "page": page,
                            "row": row, "col": col, "x1": box[0], "y1": box[1],
                            "x2": box[2], "y2": box[3], "crop_path": str(path)})

    cells = pd.DataFrame(records)
    cells.to_csv(results_dir / f"missing_cells_paper{paper_id}_page{page}.csv", index=False)

    # Review sheet so you can check what each empty cell actually contains
    cols = 6
    rows = -(-len(cells) // cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.2, rows * 2.2))
    for ax in axes.flat:
        ax.axis("off")
    for ax, (_, c) in zip(axes.flat, cells.iterrows()):
        ax.imshow(Image.open(c.crop_path))
        ax.set_title(f"row {c.row}, col {c.col}", fontsize=8)
    plt.tight_layout()
    plt.savefig(results_dir / f"missing_cells_paper{paper_id}_page{page}.png", dpi=120)
    print(f"Cropped {len(cells)} empty cells into {out_dir}/")
    return cells


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crop table cells that detection missed.")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", default="2")
    parser.add_argument("--page", type=int, default=5)
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    crop_missing_cells(args.results, cfg["data"]["pdf_dir"], args.paper, args.page)