"""Assemble a table read by grid detection (table_grid.py) and grid OCR (grid_ocr.py).

Paper 19 Table 11: R1 text is snapped to a small vocabulary, the fluorine pattern (2,4 or 2,6)
and ring atom X come from OCR, and each R2 drawing goes to the merged R2 cell it sits in.
The scaffold is Table 9's verified scaffold with its fixed propyl and chlorophenyl cut off,
opening slots for R1 and R2.
"""
import argparse
import difflib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem

from src.assemble import assemble_row
from src.pieces import load_crops
from src.text_tables import fill_ring_atoms

R1_VOCAB = {"nC3H7": "*CCC", "benzyl": "*Cc1ccccc1", "-CH3": "*C", "-C2H5": "*CC", "-(CH2)2-CF3": "*CCC(F)(F)F"}
# Table 9 slots: Ra = 3, Rb = 2, Rc = 4, Rd = 5 (verified there by self test); None means H
# Numbered from the sulfonamide side (C1 bears the NH), as in vemurafenib's IUPAC name, where
# "2,4-difluoro" puts both fluorines next to the carbonyl. Table 9 slots: Ra = 3, Rb = 2, Rc = 4, Rd = 5.
FLUORINE = {"2,4": {3: "*F", 2: None, 4: None, 5: "*F"}, "2,6": {3: "*F", 2: "*F", 4: None, 5: None}}
RING = {"CH": "c", "N": "n"}
R1_SLOT, R2_SLOT = 6, 7


def snap(text, vocab, min_margin=0.08):
    """Closest vocabulary entry; flagged if the runner-up is nearly as close."""
    text = text if isinstance(text, str) else ""
    scores = sorted(((difflib.SequenceMatcher(None, text, v).ratio(), v) for v in vocab), reverse=True)
    (best_score, best), (second_score, _) = scores[0], scores[1]
    return best, round(best_score, 3), (best_score - second_score) < min_margin


def open_slots(scaffold_smiles):
    """Cut the fixed propyl (R1) and 4-chlorophenyl (R2) off Table 9's scaffold, leaving slots 6 and 7."""
    mol = Chem.MolFromSmiles(scaffold_smiles)
    aryl = {m[:2] for m in mol.GetSubstructMatches(Chem.MolFromSmarts("c-!@c1ccc(Cl)cc1"))}
    propyl = {m[:2] for m in mol.GetSubstructMatches(Chem.MolFromSmarts("[SX4](=O)(=O)-[CH2][CH2][CH3]"))}
    if len(aryl) != 1 or len(propyl) != 1:
        raise ValueError("expected exactly one chlorophenyl and one propyl to cut")
    (a0, a1), = aryl
    s_idx, c_idx = next(iter(propyl))[0], mol.GetSubstructMatch(Chem.MolFromSmarts("[SX4](=O)(=O)-[CH2][CH2][CH3]"))[3]
    bonds = [mol.GetBondBetweenAtoms(a0, a1).GetIdx(), mol.GetBondBetweenAtoms(s_idx, c_idx).GetIdx()]
    cut = Chem.FragmentOnBonds(mol, bonds, addDummies=True, dummyLabels=[(70, 70), (60, 60)])
    core = max(Chem.GetMolFrags(cut, asMols=True), key=lambda m: m.GetNumAtoms())
    for atom in core.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() in (60, 70):
            atom.SetIsotope(R1_SLOT if atom.GetIsotope() == 60 else R2_SLOT)
    return Chem.MolToSmiles(core)


def assign_r2(grid, fragments):
    """Map each band to the R2 drawing whose center lies inside that band's merged R2 cell."""
    bands = {b["band"]: b for b in grid["bands"]}
    sections, current = [], []                       # runs of data rows between separators
    for b in grid["bands"][1:]:
        if b["separator"]:
            sections.append(current); current = []
        else:
            current.append(b["band"])
    sections.append(current)

    band_r2, group_r2 = {}, {}
    for group in grid["columns"]["R2"]:
        top, bottom = bands[group[0]]["y_top"], bands[group[-1]]["y_bottom"]
        inside = fragments[(fragments.yc > top) & (fragments.yc < bottom)]
        group_r2[tuple(group)] = inside.iloc[0] if len(inside) == 1 else None
    flags = {}
    for group, frag in group_r2.items():
        if frag is None:                             # empty cell: nearest filled group in the same section
            section = next(s for s in sections if group[0] in s)
            candidates = [(min(abs(b - g) for b in other for g in group), f) for other, f in group_r2.items()
                          if f is not None and other[0] in section]
            frag = min(candidates, key=lambda c: c[0])[1] if candidates else None
            for band in group:
                flags[band] = "R2 inferred from neighbor"
        for band in group:
            band_r2[band] = frag
    return band_r2, flags


def assemble_grid_table(results_dir, paper_id, page, scaffold_crop):
    results_dir = Path(results_dir)
    grid = json.loads((results_dir / f"grid_paper{paper_id}_page{page}.json").read_text())
    cells = pd.read_csv(results_dir / f"grid_cells_paper{paper_id}_page{page}.csv")
    smiles_of = pd.read_csv(results_dir / "validated_predictions.csv").set_index("crop_id").canonical_smiles
    core = open_slots(smiles_of[scaffold_crop])
    crops = load_crops(results_dir, paper_id, page)
    band_r2, flags = assign_r2(grid, crops[crops.status == "fragment"])

    records = []
    for c in cells.itertuples():
        notes = [flags[c.band]] if c.band in flags else []
        r1, r1_score, r1_close = snap(c.R1, R1_VOCAB)
        if r1_close:
            notes.append(f"R1 '{c.R1}' ambiguous")
        pattern = "2,4" if "2,4" in str(c.A) else "2,6" if "2,6" in str(c.A) else None
        ring = RING.get(str(c.X))
        frag = band_r2.get(c.band)
        smiles, status = None, "assembled"
        try:
            if pattern is None or ring is None or frag is None:
                raise ValueError(f"unreadable: A={c.A!r} X={c.X!r} R2={'ok' if frag is not None else 'none'}")
            slots = dict(FLUORINE[pattern])
            slots[R1_SLOT], slots[R2_SLOT] = R1_VOCAB[r1], frag.canonical_smiles
            smiles = assemble_row(fill_ring_atoms(core, {1: ring}), slots, "C")
        except Exception as exc:
            status = f"error: {exc}"
        records.append({"paper_id": str(paper_id), "page": page, "band": c.band, "R1": r1,
                        "R1_ocr": c.R1, "A": pattern, "X": c.X,
                        "R2_crop": frag.crop_id if frag is not None else None,
                        "assembled_smiles": smiles, "flags": "; ".join(notes), "status": status})

    out = pd.DataFrame(records)
    out.to_csv(results_dir / f"assembled_grid_paper{paper_id}_page{page}.csv", index=False)
    print(out[["band", "R1_ocr", "R1", "A", "X", "R2_crop", "flags", "status"]].to_string(index=False))
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble a grid-detected table (paper 19 Table 11).")
    parser.add_argument("--results", default="results/A300")
    args = parser.parse_args()
    assemble_grid_table(args.results, "19", 9, scaffold_crop="19_p8_001")