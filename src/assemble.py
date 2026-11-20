"""Assemble complete molecules from a scaffold plus a substituent table drawn as R-group cells.

For one paper and page: group fragment crops into table rows and columns by position,
read each row's compound number and variable atom (X) from the page text, and join the
scaffold with that row's fragments. Every empty cell assumed to be H is flagged.
"""
import argparse
import re
from pathlib import Path

import pandas as pd
import pymupdf
import yaml
from rdkit import Chem, RDLogger
from src.validate import classify_smiles

RDLogger.DisableLog("rdApp.*")

ROW_GAP_PX = 30      # crops within 30 px vertically belong to the same table row
COL_GAP_PX = 40      # crops within 40 px horizontally belong to the same column
H_FRAGMENTS = {"*[H]", "[H]*"}


def cluster(values, gap):
    """Group 1D positions: start a new group when the next sorted value jumps by more than gap."""
    order = sorted(set(values))
    labels, current = {}, 0
    for i, val in enumerate(order):
        if i and val - order[i - 1] > gap:
            current += 1
        labels[val] = current
    return [labels[v] for v in values]


def fill_ring_variable(scaffold_smiles, element):
    """Replace a variable ring atom (like X, written [5*]) with aromatic c or n.

    Returns (smiles, removed_slot). For N, the substituent slot on X is dropped,
    because a pyridine nitrogen carries no substituent.
    """
    mol = Chem.MolFromSmiles(scaffold_smiles)
    ring_dummies = [a.GetIsotope() for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.IsInRing()]
    if not ring_dummies:
        return scaffold_smiles, None
    token = f"[{ring_dummies[0]}*]"
    smiles, removed = scaffold_smiles, None
    if element == "N":
        match = re.search(re.escape(token) + r"\(\[(\d+)\*\]\)", smiles)
        if match:
            removed = int(match.group(1))
            smiles = smiles.replace(match.group(0), token, 1)
    return smiles.replace(token, "n" if element == "N" else "c", 1), removed


def with_map_numbers(smiles, slot=None):
    """Turn wildcard labels into atom map numbers, which molzip uses to pair attachment points."""
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomMapNum(slot if slot is not None else atom.GetIsotope())
            atom.SetIsotope(0)
    return mol


def assemble_row(scaffold_smiles, fragments_by_slot, x_element):
    """Join scaffold + fragments; empty or H slots become hydrogen. Returns canonical SMILES."""
    filled, removed_slot = fill_ring_variable(scaffold_smiles, x_element)
    mol = with_map_numbers(filled)
    for slot, frag in fragments_by_slot.items():
        if slot == removed_slot or frag is None or frag in H_FRAGMENTS:
            continue
        mol = Chem.molzip(mol, with_map_numbers(frag, slot))
    # Any wildcard left over (H cells, empty cells) is removed so it becomes an implicit H
    rw = Chem.RWMol(mol)
    for idx in sorted((a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0), reverse=True):
        rw.RemoveAtom(idx)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out)


def read_row_text(pdf_path, page_no, row_centers, first_col_x, last_col_x, dpi=200):
    """For each row, find the compound number left of the table and the X value (C/N) to its right."""
    scale = dpi / 72  # PDF text positions are in points; crop boxes are pixels at 200 DPI
    words = pymupdf.open(pdf_path)[page_no - 1].get_text("words")
    info = {}
    for row, y in row_centers.items():
        name, x_elem = None, "C"
        for x0, y0, x1, y1, word, *_ in words:
            xc, yc = (x0 + x1) / 2 * scale, (y0 + y1) / 2 * scale
            if abs(yc - y) > ROW_GAP_PX:
                continue
            if word.isdigit() and xc < first_col_x - 20:
                name = word
            elif word in ("C", "N") and xc > last_col_x + 40:
                x_elem = word
        info[row] = (name, x_elem)
    return info


def assemble_table(results_dir, pdf_dir, paper_id, page):
    results_dir = Path(results_dir)
    crops = pd.read_csv(results_dir / "crops.csv", dtype={"paper_id": str})
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    m = crops.merge(validated[["crop_id", "canonical_smiles", "status"]], on="crop_id")
    page_rows = m[(m.paper_id == str(paper_id)) & (m.page == page)]

    scaffold = page_rows[page_rows.status == "scaffold"].sort_values("y1").iloc[0]
    frags = page_rows[page_rows.status == "fragment"].copy()
    frags["xc"], frags["yc"] = (frags.x1 + frags.x2) / 2, (frags.y1 + frags.y2) / 2
    frags["row"] = cluster(list(frags.yc), ROW_GAP_PX)
    frags["col"] = cluster(list(frags.xc), COL_GAP_PX)

    col_x = frags.groupby("col").xc.mean().sort_index()
    slots = sorted(a.GetIsotope() for a in Chem.MolFromSmiles(scaffold.canonical_smiles).GetAtoms()
                   if a.GetAtomicNum() == 0 and not a.IsInRing())
    if len(col_x) != len(slots):
        raise ValueError(f"{len(col_x)} columns found but scaffold has {len(slots)} R slots")
    col_to_slot = dict(zip(col_x.index, slots))  # leftmost column = lowest R number

    row_y = frags.groupby("row").yc.mean()
    text = read_row_text(Path(pdf_dir) / f"{paper_id}.pdf", page, row_y.to_dict(),
                         col_x.iloc[0], col_x.iloc[-1])

    # Prefer targeted OCR results when available: compound number from row order, X from OCR
    table_text_csv = results_dir / f"table_text_paper{paper_id}_page{page}.csv"
    if table_text_csv.exists():
        for r in pd.read_csv(table_text_csv).itertuples():
            if r.row in text:
                x_atom = r.x_atom if r.x_atom in ("C", "N") else "C"
                text[r.row] = (str(r.compound_inferred), x_atom)

    # Optional: predictions for table cells that detection missed (from fill_cells + recognize)
    cell_map = {}
    cells_csv = results_dir / f"missing_cells_paper{paper_id}_page{page}.csv"
    preds_csv = results_dir / f"cell_predictions_paper{paper_id}_page{page}.csv"
    if cells_csv.exists() and preds_csv.exists():
        filled = pd.read_csv(cells_csv).merge(
            pd.read_csv(preds_csv)[["crop_id", "raw_smiles"]],
            left_on="cell_id", right_on="crop_id")
        for c in filled.itertuples():
            canonical, status = classify_smiles(c.raw_smiles)
            if status == "fragment":
                cell_map[(c.row, c.col)] = canonical
        print(f"Using {len(cell_map)} recognized cells to fill gaps")

    records = []
    for row, row_crops in frags.groupby("row"):
        by_slot = {slot: None for slot in slots}
        for _, c in row_crops.iterrows():
            by_slot[col_to_slot[c.col]] = c.canonical_smiles
        filled_from_cells = []
        for col, slot in col_to_slot.items():
            if by_slot[slot] is None and (row, col) in cell_map:
                by_slot[slot] = cell_map[(row, col)]
                filled_from_cells.append(f"R{slot}")
        assumed = [f"R{s}" for s, f in by_slot.items() if f is None]
        name, x_elem = text[row]
        try:
            smiles, error = assemble_row(scaffold.canonical_smiles, by_slot, x_elem), ""
        except Exception as exc:  # keep going; record why this row failed
            smiles, error = None, str(exc)
        records.append({
            "paper_id": str(paper_id), "page": page, "row": row, "compound": name, "x_atom": x_elem,
            "assembled_smiles": smiles, "filled_from_cells": ",".join(filled_from_cells),
            "assumed_h": ",".join(assumed), "scaffold_crop": scaffold.crop_id,
            "fragment_crops": ",".join(row_crops.crop_id), "error": error,
        })

    out = pd.DataFrame(records)
    out_path = results_dir / f"assembled_paper{paper_id}_page{page}.csv"
    out.to_csv(out_path, index=False)
    print(f"Assembled {out.assembled_smiles.notna().sum()} of {len(out)} rows -> {out_path}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble molecules from a scaffold and R-group table.")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", default="2")
    parser.add_argument("--page", type=int, default=5)
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    assemble_table(args.results, cfg["data"]["pdf_dir"], args.paper, args.page)