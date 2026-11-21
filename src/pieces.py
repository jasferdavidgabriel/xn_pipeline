"""Assemble molecules from separately drawn pieces, as in paper 4's PROTAC tables.

Table 3 style: a warhead with one slot (R) combined with a second piece (OMe, or the
linker + ligand piece "A"). Table 1 style: warhead + two-ended linker + ligand.
Compound names are read from the PDF text layer and matched to pieces by position.
"""
import argparse
import re
from pathlib import Path

import pandas as pd
import pymupdf
import yaml
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
SCALE = 200 / 72   # crop boxes are pixels at 200 DPI; PDF words are in points


def load_crops(results_dir, paper_id, page):
    crops = pd.read_csv(Path(results_dir) / "crops.csv", dtype={"paper_id": str})
    v = pd.read_csv(Path(results_dir) / "validated_predictions.csv", dtype={"paper_id": str})
    m = crops.merge(v[["crop_id", "status", "canonical_smiles"]], on="crop_id")
    m = m[(m.paper_id == str(paper_id)) & (m.page == page)].copy()
    m["xc"], m["yc"] = (m.x1 + m.x2) / 2, (m.y1 + m.y2) / 2
    return m


def page_words(pdf_path, page, pattern):
    """Words matching pattern, with centers converted to crop pixels."""
    words = pymupdf.open(pdf_path)[page - 1].get_text("words")
    rx = re.compile(pattern)
    return [((x0 + x1) / 2 * SCALE, (y0 + y1) / 2 * SCALE, w)
            for x0, y0, x1, y1, w, *_ in words if rx.match(w)]


def with_map(smiles, map_no):
    """Label every wildcard of a one-slot piece with map_no, for molzip."""
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetIsotope(0)
            a.SetAtomMapNum(map_no)
    return mol


def oriented_linker(smiles, warhead_map, ligand_map):
    """The linker end bonded to N joins the ligand's carbonyl; the other end joins the warhead."""
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    ends = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(ends) != 2:
        raise ValueError("linker needs exactly two ends")
    n_ends = [a.GetIdx() for a in ends if a.GetNeighbors()[0].GetSymbol() == "N"]
    if len(n_ends) != 1:
        raise ValueError("cannot orient linker: need exactly one end on N")
    for a in ends:
        a.SetIsotope(0)
        a.SetAtomMapNum(ligand_map if a.GetIdx() == n_ends[0] else warhead_map)
    return mol


def acyl_from_ligand_piece(smiles):
    """Cut a 'linker + ligand' piece at its amide N-C(=O) bond and keep the acyl side, *C(=O)Ar."""
    mol = Chem.MolFromSmiles(smiles)
    matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3][CX3](=O)c"))
    if len(matches) != 1:
        raise ValueError("expected exactly one aryl amide to cut")
    bond = mol.GetBondBetweenAtoms(matches[0][0], matches[0][1]).GetIdx()
    pieces = Chem.GetMolFrags(Chem.FragmentOnBonds(mol, [bond]), asMols=True)
    acyl = [p for p in pieces if p.HasSubstructMatch(Chem.MolFromSmarts("[#0]C(=O)c"))]
    return Chem.MolToSmiles(acyl[0])


def zip_all(*mols):
    out = mols[0]
    for m in mols[1:]:
        out = Chem.molzip(out, m)
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out)


def warhead_table(results_dir, pdf_dir, paper_id, page, piece_crop, name_pattern,
                  ligand_pattern=r"^5[2-6]$", small_group="*OC"):
    """Warheads, each combined with small_group (e.g. OMe) and with the big piece.

    The labels under the warheads are part of an image, so names come from the data table
    below the figure: each row pairs a ligand (warhead + small_group) with its PROTAC
    (warhead + piece). Rows run in the same left to right order as the warheads; the
    answer key and compound 44's second route (Table 1) verify that assumption.
    """
    crops = load_crops(results_dir, paper_id, page)
    piece = crops[crops.crop_id == piece_crop].iloc[0]
    warheads = crops[(crops.status == "fragment") & (crops.y2 < piece.y1)].sort_values("xc")

    rows = {}
    for x, y, w in page_words(Path(pdf_dir) / f"{paper_id}.pdf", page, name_pattern):
        if y > piece.y2:                                   # the data table below the figure
            rows.setdefault(round(y / 8), []).append((x, w))
    pairs = []
    for key in sorted(rows):
        names = [w for _, w in sorted(rows[key])]          # left to right within the row
        if len(names) >= 2 and re.match(ligand_pattern, names[0]):
            pairs.append((names[0], names[1]))

    if len(pairs) != len(warheads):
        raise ValueError(f"{len(pairs)} table rows but {len(warheads)} warheads")
    records = []
    for (_, wh), (ligand, protac) in zip(warheads.iterrows(), pairs):
        for name, other, label in ((ligand, small_group, "OMe"), (protac, piece.canonical_smiles, piece_crop)):
            records.append({"paper_id": str(paper_id), "page": page, "compound": name,
                            "assembled_smiles": zip_all(with_map(wh.canonical_smiles, 1), with_map(other, 1)),
                            "pieces": f"{wh.crop_id}+{label}"})
    return pd.DataFrame(records)

def linker_table(results_dir, pdf_dir, paper_id, page, warhead_smiles, ligand_piece_smiles, name_pattern):
    """Warhead + each two-ended linker on the page + the ligand's acyl piece."""
    crops = load_crops(results_dir, paper_id, page)
    linkers = crops[crops.status == "linker"]
    labels = page_words(Path(pdf_dir) / f"{paper_id}.pdf", page, name_pattern)
    acyl = acyl_from_ligand_piece(ligand_piece_smiles)

    records = []
    for _, lk in linkers.iterrows():
        row_labels = [w for w in labels if abs(w[1] - lk.yc) < 40 and w[0] < lk.x1]
        name = min(row_labels, key=lambda w: abs(w[1] - lk.yc))[2] if row_labels else None
        try:
            smiles = zip_all(with_map(warhead_smiles, 1), oriented_linker(lk.canonical_smiles, 1, 2),
                             with_map(acyl, 2))
            status = "assembled"
        except Exception as exc:
            smiles, status = None, f"error: {exc}"
        records.append({"paper_id": str(paper_id), "page": page, "compound": name,
                        "assembled_smiles": smiles, "pieces": lk.crop_id, "status": status})
    return pd.DataFrame(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble paper 4 PROTAC tables from drawn pieces.")
    parser.add_argument("--results", default="results/A300")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    pdf_dir, results = cfg["data"]["pdf_dir"], Path(args.results)

    # Table 3 (page 6): five warheads, each with OMe and with piece A
    t3 = warhead_table(results, pdf_dir, "4", 6, "4_p6_001", r"^(5[2-6]|4[4-9]|CCT400028)$")
    t3.to_csv(results / "assembled_pieces_paper4_page6.csv", index=False)
    print(f"Table 3: {len(t3)} compounds"); print(t3[["compound", "pieces"]].to_string(index=False))

    # Table 1 (page 4): its scaffold failed recognition, so the thalidomide warhead and the
    # ligand piece are taken from Table 3's drawings on page 6 (the same structures).
    v = pd.read_csv(results / "validated_predictions.csv")
    get = lambda crop: v.loc[v.crop_id == crop, "canonical_smiles"].iloc[0]
    t1 = linker_table(results, pdf_dir, "4", 4, get("4_p6_004"), get("4_p6_001"), r"^(JB\d+|4[2-6]|51)$")
    t1.to_csv(results / "assembled_pieces_paper4_page4.csv", index=False)
    print(f"\nTable 1: {len(t1)} compounds"); print(t1[["compound", "pieces", "status"]].to_string(index=False))