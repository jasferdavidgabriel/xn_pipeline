"""Find a table's grid lines in a page image: rows, columns, and (later) merged cells.

A horizontal line is a pixel row where a large share of pixels is dark; a vertical line is a
pixel column that is dark for a large share of the table's height. Positions are in 200 DPI
pixels, the same units as the crop boxes.
"""
import argparse
from pathlib import Path

import numpy as np
import pymupdf
import yaml
from PIL import Image, ImageDraw

DPI = 200
DARK = 200          # table rules are often drawn in light gray, so count anything darker than near-white


def runs(indices):
    """Group consecutive pixel indices (a line a few pixels thick) into one center each."""
    groups, start = [], None
    for i, v in enumerate(indices):
        if start is None:
            start = v
        if i == len(indices) - 1 or indices[i + 1] != v + 1:
            groups.append((start + v) // 2)
            start = None
    return groups


def render_gray(pdf_path, page):
    pix = pymupdf.open(pdf_path)[page - 1].get_pixmap(dpi=DPI)
    return Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("L")


def detect_grid(gray, box, min_frac_h=0.25, min_frac_v=0.6, bottom=None):
    """Return (horizontal line ys, vertical line xs) inside box = (x_from, y_from, x_to, y_to).

    Horizontal lines use a low threshold because merged cells interrupt them. Vertical lines use
    a high one: a real column border runs the table's full height, while a stack of identical
    characters (like "2," on every row) is broken at every row gap.
    """
    x0, y0, x1, y1 = box
    dark = np.array(gray)[y0:y1, x0:x1] < DARK
    h_lines = [y + y0 for y in runs(np.where(dark.mean(axis=1) > min_frac_h)[0])]
    if len(h_lines) < 2:
        return h_lines, []
    table = dark[h_lines[0] - y0:h_lines[-1] - y0]
    v_lines = [x + x0 for x in runs(np.where(table.mean(axis=0) > min_frac_v)[0])]
    # Some tables have no bottom rule: the last row ends where the column borders stop
    border_cols = [x - x0 for x in v_lines]
    still_bordered = np.where(dark[:, border_cols].mean(axis=1) > 0.6)[0]
    if len(still_bordered) and still_bordered.max() + y0 > h_lines[-1] + 20:
        h_lines.append(int(still_bordered.max() + y0))
    return h_lines, v_lines

def bands_from_lines(h_lines, min_height=30):
    """Bands between consecutive lines. Thin bands are separators (like the gray bars)."""
    return [{"band": i, "y_top": int(a), "y_bottom": int(b), "separator": bool((b - a) < min_height)}
            for i, (a, b) in enumerate(zip(h_lines, h_lines[1:]))]

def has_border(dark, y, x_from, x_to, frac=0.5):
    """True if a ruled line crosses this column at height y (checked 1 pixel above and below)."""
    strip = dark[y - 1:y + 2, x_from + 3:x_to - 3]   # stay clear of the column's own vertical borders
    return strip.any(axis=0).mean() > frac


def merged_groups(gray, h_lines, v_lines, names, header_bands=1):
    """For each named column, group data rows that share one merged cell."""
    dark = np.array(gray) < DARK
    bands = bands_from_lines(h_lines)
    rows = [b for b in bands[header_bands:] if not b["separator"]]
    columns = {}
    for name, x_from, x_to in zip(names, v_lines, v_lines[1:]):
        groups, current = [], [rows[0]["band"]]
        for prev, cur in zip(rows, rows[1:]):
            touching = prev["y_bottom"] == cur["y_top"]          # false across a gray bar
            if touching and not has_border(dark, cur["y_top"], x_from, x_to):
                current.append(cur["band"])
            else:
                groups.append(current)
                current = [cur["band"]]
        groups.append(current)
        columns[name] = groups
    return bands, columns

def preview(gray, h_lines, v_lines, path):
    img = gray.convert("RGB")
    draw = ImageDraw.Draw(img)
    for y in h_lines:
        draw.line([(0, y), (img.width, y)], fill=(220, 0, 0), width=2)
    for x in v_lines:
        draw.line([(x, h_lines[0]), (x, h_lines[-1])], fill=(0, 90, 220), width=2)
    for i, (a, b) in enumerate(zip(h_lines, h_lines[1:])):
        draw.text((v_lines[0] - 40 if v_lines else 5, (a + b) // 2 - 6), str(i), fill=(220, 0, 0))
    top, bottom = max(h_lines[0] - 40, 0), min(h_lines[-1] + 40, img.height)
    img.crop((0, top, img.width, bottom)).save(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect and preview a table's grid lines.")
    parser.add_argument("--results", default="results/A300")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", required=True)
    parser.add_argument("--page", type=int, required=True)
    parser.add_argument("--box", required=True, help="Table area x_from,y_from,x_to,y_to (200 DPI pixels)")
    parser.add_argument("--names", help="Comma separated column names, left to right")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    gray = render_gray(Path(cfg["data"]["pdf_dir"]) / f"{args.paper}.pdf", args.page)
    h, v = detect_grid(gray, [int(n) for n in args.box.split(",")])
    print(f"{len(h)} horizontal lines: {h}")
    print(f"{len(v)} vertical lines:   {v}")
    out = Path(args.results) / f"grid_paper{args.paper}_page{args.page}.png"
    preview(gray, h, v, out)
    print(f"Preview saved to {out}")
    if args.names:
        import json
        names = args.names.split(",")
        bands, columns = merged_groups(gray, h, v, names)
        for name, groups in columns.items():
            merged = [g for g in groups if len(g) > 1]
            print(f"{name:>8}: {len(groups)} cells; merged: {merged if merged else 'none'}")
        grid = {"h_lines": [int(y) for y in h], "v_lines": [int(x) for x in v],
                "bands": bands, "columns": columns}
        out_json = Path(args.results) / f"grid_paper{args.paper}_page{args.page}.json"
        out_json.write_text(json.dumps(grid, indent=1))
        print(f"Grid saved to {out_json}")