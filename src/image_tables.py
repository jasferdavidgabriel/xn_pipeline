"""Assemble rows of an image table whose text columns were read by targeted OCR (read_columns.py).

OCR on tiny cells is unreliable, so each column combines OCR with the table's own logic:
entries keep their last two digits (a stray leading digit comes from the table border) and a
single unreadable entry is filled from a consecutive run; a dash in the ring column (found by
ink height) means the scaffold without a ring slot. Every inferred value is flagged.
"""
import argparse
from pathlib import Path

import pandas as pd
import yaml

from src.assemble import assemble_row
from src.text_tables import attach_single, fill_ring_atoms

DASH_INK = 0.1   # a dash cell's ink height is ~0.03; any letter is ~0.24


def clean_entries(raw):
    """Return (entries, inferred_flags)."""
    vals = [s[-2:] if isinstance(s, str) and len(s) >= 2 and s[-2:].isdigit() else None for s in raw]
    known = [(i, int(v)) for i, v in enumerate(vals) if v]
    offsets = {v - i for i, v in known}
    if known and len(offsets) == 1:                 # the readable entries form one consecutive run
        start = offsets.pop()
        return [str(start + i) for i in range(len(vals))], [v is None for v in vals]
    return vals, [v is None for v in vals]


def ring_value(text, ink):
    """Aromatic token for the ring column: 'n', 'c', None for a dash, or '?' if unreadable."""
    text = text if isinstance(text, str) else ""
    if ink < DASH_INK:
        return None
    if "N" in text and "C" not in text:
        return "n"
    if "CH" in text:
        return "c"
    return "?"


def assemble_image_table(results_dir, paper_id, page, scaffolds, ring_slot, r_slot):
    """scaffolds: {"A": crop with a ring slot, "B": crop without one}."""
    results_dir = Path(results_dir)
    cells = pd.read_csv(results_dir / f"table_cells_paper{paper_id}_page{page}.csv",
                        dtype={"entry": str, "scaffold": str, "A": str})
    smiles_of = pd.read_csv(results_dir / "validated_predictions.csv").set_index("crop_id").canonical_smiles
    names, inferred = clean_entries(cells.entry.tolist())

    records = []
    for (_, c), name, was_inferred in zip(cells.iterrows(), names, inferred):
        ring = ring_value(c.A, c.A_ink)
        letter = "B" if ring is None else "A"
        ocr_letter = c.scaffold if isinstance(c.scaffold, str) else ""
        flags = []
        if was_inferred:
            flags.append("entry inferred")
        if ocr_letter and ocr_letter != letter:
            flags.append(f"OCR scaffold {ocr_letter} disagrees")
        smiles, status = None, "assembled"
        try:
            if ring == "?":
                raise ValueError("ring column unreadable")
            fragment = smiles_of[c.crop_id]
            if letter == "A":
                filled = fill_ring_atoms(smiles_of[scaffolds["A"]], {ring_slot: ring})
                smiles = assemble_row(filled, {r_slot: fragment}, "C")
            else:
                smiles = attach_single(smiles_of[scaffolds["B"]], fragment)
        except Exception as exc:
            status = f"error: {exc}"
        records.append({"paper_id": str(paper_id), "page": page, "compound": name, "scaffold": letter,
                        "ring": ring or "-", "fragment_crop": c.crop_id, "assembled_smiles": smiles,
                        "flags": "; ".join(flags), "status": status})

    out = pd.DataFrame(records)
    out.to_csv(results_dir / f"assembled_image_paper{paper_id}_page{page}.csv", index=False)
    print(out[["compound", "scaffold", "ring", "fragment_crop", "flags", "status"]].to_string(index=False))
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble an OCR-read image table with two scaffolds.")
    parser.add_argument("--results", default="results/A300")
    args = parser.parse_args()
    # Paper 19 Table 7: scaffold A has ring atom [1*] and R slot [2*]; scaffold B has one R slot
    assemble_image_table(args.results, "19", 7, {"A": "19_p7_002", "B": "19_p7_001"}, ring_slot=1, r_slot=2)