"""Compare two front ends crop by crop and molecule by molecule.

Crops are matched by overlap on the same page (both use 200 DPI page coordinates). Also checks
whether a front end found the table cells that full-page detection missed in paper 2.
"""
import argparse
from pathlib import Path

import pandas as pd


def iou(a, b):
    ix = max(0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0, min(a[3], b[3]) - max(a[1], b[1]))
    inter = ix * iy
    union = (a[2] - a[0]) * (a[3] - a[1]) + (b[2] - b[0]) * (b[3] - b[1]) - inter
    return inter / union if union else 0.0


def boxes(df):
    return list(zip(df.x1, df.y1, df.x2, df.y2))


def match_crops(a, b, threshold=0.3):
    """Per paper: crops in A, in B, matched in both, only in A, only in B."""
    rows = []
    for paper in sorted(set(a.paper_id) | set(b.paper_id), key=int):
        found_a = found_b = 0
        pa, pb = a[a.paper_id == paper], b[b.paper_id == paper]
        for page in set(pa.page) | set(pb.page):
            ba, bb = boxes(pa[pa.page == page]), boxes(pb[pb.page == page])
            found_a += sum(any(iou(x, y) >= threshold for y in bb) for x in ba)
            found_b += sum(any(iou(y, x) >= threshold for x in ba) for y in bb)
        rows.append({"paper": paper, "A crops": len(pa), "B crops": len(pb),
                     "A found by B": found_a, "only in A": len(pa) - found_a, "only in B": len(pb) - found_b})
    return pd.DataFrame(rows)


def cells_found(cells, b):
    """How many missed table cells have their center inside some B crop on the same page."""
    hits = 0
    for c in cells.itertuples():
        cx, cy = (c.x1 + c.x2) / 2, (c.y1 + c.y2) / 2
        same = b[(b.paper_id == str(c.paper_id)) & (b.page == c.page)]
        hits += any(x1 <= cx <= x2 and y1 <= cy <= y2 for x1, y1, x2, y2 in boxes(same))
    return hits


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare two front ends.")
    parser.add_argument("--a", default="results/G")
    parser.add_argument("--b", default="results/B")
    parser.add_argument("--cells", default="results/A300/missing_cells_paper2_page5.csv")
    args = parser.parse_args()
    A, B = Path(args.a), Path(args.b)
    ca = pd.read_csv(A / "crops.csv", dtype={"paper_id": str})
    cb = pd.read_csv(B / "crops.csv", dtype={"paper_id": str})

    print("Crops matched by position (overlap >= 0.3):")
    print(match_crops(ca, cb).to_string(index=False))

    cells = pd.read_csv(args.cells)
    if {"x1", "y1", "x2", "y2"} <= set(cells.columns):
        cells["paper_id"], cells["page"] = cells.get("paper_id", "2"), cells.get("page", 5)
        print(f"\nPaper 2 table cells missed by full-page detection: {len(cells)}; "
              f"found by B: {cells_found(cells, cb)}")
    else:
        print("\nMissed-cell file has no box columns:", list(cells.columns))

    ma = pd.read_csv(A / "final_per_molecule.csv", dtype={"paper_id": str})
    mb = pd.read_csv(B / "final_per_molecule.csv", dtype={"paper_id": str})
    m = ma.merge(mb, on=["paper_id", "molecule_name"], suffixes=("_A", "_B"))
    cols = ["paper_id", "molecule_name", "best_tanimoto_A", "best_tanimoto_B"]
    print("\nFound only by A:"); print(m[m.exact_match_A & ~m.exact_match_B][cols].to_string(index=False))
    print("\nFound only by B:"); print(m[~m.exact_match_A & m.exact_match_B][cols].to_string(index=False))