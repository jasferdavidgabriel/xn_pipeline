"""Assemble compounds defined by text labels in a scheme, like '19 R = Br; R1 = H' (paper 4, Scheme 1).

Labels come from full-page OCR (page_ocr.py). Each label belongs to the scaffold whose vertical span
contains it and that sits immediately to its left, as schemes are laid out. Values come from a small
vocabulary; an empty value (e.g. '12 R==') means the substituent is drawn, i.e. the chain crop.
"""
import argparse
import re
from pathlib import Path

import pandas as pd
from rdkit import Chem

from src.expand_abbreviations import ABBREVIATIONS, ANNOTATION, expand
from src.pieces import load_crops

VALUES = {"BOC": ABBREVIATIONS["Boc"], "BR": "*Br", "H": "*[H]", "BPIN": "*B1OC(C)(C)C(C)(C)O1",
          "SEM": ABBREVIATIONS["SEM"], "DMB": ABBREVIATIONS["DMB"]}
KEY_FIX = {"R": "R", "R1": "R1", "RI": "R1", "RL": "R1", "R'": "R'"}   # OCR often reads subscript 1 as I or l
MOLGLYPH_KEY = {"R": "R", "R[1]": "R1", "R'": "R'"}                 # how MolGlyph labels the same slots


def parse_label(text):
    """'19 R =Br; RI=H' -> ('19', {'R': 'BR', 'R1': 'H'}). An empty value means a drawn substituent."""
    m = re.match(r"^\s*(\d+)\s+(.*)$", str(text))
    if not m:
        return None
    number, body = m.groups()
    assignments = {}
    for part in re.split(r"[;,]", body):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key = KEY_FIX.get(key.strip().upper().replace("’", "'"))
        if key:
            assignments[key] = value.replace("=", "").strip().upper() or "DRAWN"
    return (number, assignments) if assignments else None


def assemble_scheme_labels(results_dir, paper_id, page, chain_crop, x_tolerance=30):
    results_dir = Path(results_dir)
    ocr = pd.read_csv(results_dir / f"page_ocr_paper{paper_id}_page{page}.csv")
    raw = pd.read_csv(results_dir / "raw_predictions.csv").set_index("crop_id").raw_output
    crops = load_crops(results_dir, paper_id, page)
    chain = crops.set_index("crop_id").canonical_smiles[chain_crop]

    records = []
    for o in ocr.itertuples():
        parsed = parse_label(o.text)
        if not parsed:
            continue
        number, assign = parsed
        yc = (o.y1 + o.y2) / 2
        left = crops[(crops.y1 <= yc) & (crops.y2 >= yc) & (crops.x2 <= o.x1 + x_tolerance)]
        if left.empty:
            continue
        scaffold = left.loc[(o.x1 - left.x2).idxmin()]          # the nearest structure to the label's left
        labels = [lab for _, lab in ANNOTATION.findall(str(raw[scaffold.crop_id]).split("<sep>", 1)[-1])]
        extra, problems = {}, []
        for lab in labels:
            key = MOLGLYPH_KEY.get(lab)
            if key is None:
                problems.append(f"unresolved label {lab}")
            elif key not in assign:
                problems.append(f"no value for {key}")
            else:
                value = assign[key]
                extra[lab] = chain if value == "DRAWN" else VALUES.get(value)
                if extra[lab] is None:
                    problems.append(f"unknown value {value}")
        smiles = None
        if not problems:
            out, status, _ = expand(raw[scaffold.crop_id], extra=extra)
            if status == "expanded_complete":
                smiles = Chem.MolToSmiles(Chem.MolFromSmiles(out))   # re-parse to drop explicit H atoms
            else:
                problems.append(status)
        records.append({"paper_id": str(paper_id), "page": page, "compound": number,
                        "label_ocr": o.text, "scaffold_crop": scaffold.crop_id,
                        "assembled_smiles": smiles, "problems": "; ".join(problems)})

    out = pd.DataFrame(records)
    out.to_csv(results_dir / f"assembled_labels_paper{paper_id}_page{page}.csv", index=False)
    print(out[["compound", "label_ocr", "scaffold_crop", "problems"]].to_string(index=False))
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble compounds from text labels in a scheme.")
    parser.add_argument("--results", default="results/A300")
    args = parser.parse_args()
    assemble_scheme_labels(args.results, "4", 3, chain_crop="4_p3_038")