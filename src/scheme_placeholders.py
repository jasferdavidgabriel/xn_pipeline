"""Resolve placeholders in a synthesis scheme with the scheme's own definition box (paper 4, Scheme 1).

The box draws linkers A to D, each with one attachment point. In the scheme, MolGlyph labels:
- a wildcard with several box letters (e.g. "ABD" for "A/B/C/D"): every box linker, one molecule each;
- a wildcard with one box letter (e.g. "C"): that linker;
- an unlabeled wildcard (<dum>) printed in blue: a placeholder letter too. Every placeholder used
  inside this scheme's structures is C, so it is resolved as C and flagged.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf
import yaml
from PIL import Image

from src.expand_abbreviations import ANNOTATION, expand


def blue_pixel_count(page_rgb, box):
    """Pixels that are clearly blue (blue much stronger than red and green), inside box."""
    a = np.array(page_rgb.crop(box)).astype(int)
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    return int(((b - r > 60) & (b - g > 30) & (b > 120)).sum())


def resolve_scheme(results_dir, pdf_dir, paper_id, page, box_crops, letters="ABCD",
                   blue_letter="C", min_blue=15):
    results_dir = Path(results_dir)
    raw = pd.read_csv(results_dir / "raw_predictions.csv").set_index("crop_id").raw_output
    crops = pd.read_csv(results_dir / "crops.csv", dtype={"paper_id": str})
    smiles_of = pd.read_csv(results_dir / "validated_predictions.csv").set_index("crop_id").canonical_smiles
    linkers = {letter: smiles_of[crop] for letter, crop in zip(letters, box_crops)}

    pix = pymupdf.open(Path(pdf_dir) / f"{paper_id}.pdf")[page - 1].get_pixmap(dpi=200)
    page_rgb = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    here = crops[(crops.paper_id == str(paper_id)) & (crops.page == page) & ~crops.crop_id.isin(box_crops)]

    records = []
    def add(crop_id, smiles, status, notes, flag=""):
        if status == "expanded_complete":
            records.append({"paper_id": str(paper_id), "page": page, "crop_id": crop_id,
                            "assembled_smiles": smiles, "labels": notes, "flags": flag})

    for c in here.itertuples():
        out = raw.get(c.crop_id)
        if not isinstance(out, str) or "<sep>" not in out:
            continue
        labels = [lab for _, lab in ANNOTATION.findall(out.split("<sep>", 1)[1])]
        multi = [lab for lab in labels if len(lab.replace("/", "")) > 1 and set(lab.replace("/", "")) <= set(letters)]
        single = [lab for lab in labels if lab in letters]
        if multi:                                            # e.g. "ABD": one molecule per box linker
            for letter in letters:
                add(c.crop_id, *expand(out, extra={multi[0]: linkers[letter]}), flag=f"linker {letter}")
        elif single:                                         # e.g. "C": that box linker
            add(c.crop_id, *expand(out, extra={single[0]: linkers[single[0]]}))
        elif labels.count("<dum>") == 1:
            box = (int(c.x1), int(c.y1), int(c.x2), int(c.y2))
            if blue_pixel_count(page_rgb, box) >= min_blue:  # an unlabeled blue placeholder
                add(c.crop_id, *expand(out, dum_fragment=linkers[blue_letter]),
                    flag=f"blue placeholder assumed {blue_letter}")

    out = pd.DataFrame(records)
    out.to_csv(results_dir / f"assembled_scheme_paper{paper_id}_page{page}.csv", index=False)
    print(out[["crop_id", "labels", "flags"]].to_string(index=False) if len(out) else "nothing resolved")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Resolve scheme placeholders from the scheme's own box.")
    parser.add_argument("--results", default="results/A300")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    # Paper 4 Scheme 1: the box draws linkers A, B, C, D from top to bottom
    resolve_scheme(args.results, cfg["data"]["pdf_dir"], "4", 3,
                   box_crops=["4_p3_040", "4_p3_041", "4_p3_039", "4_p3_031"])