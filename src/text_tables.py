"""Assemble molecules from tables whose R-groups are written as text (e.g. paper 19, Table 6).

Reads the table rows straight from the PDF text layer, maps each R text to a SMILES
fragment with a curated dictionary, and attaches it to a single-slot scaffold crop.
"""
import argparse
import re
from pathlib import Path

import pandas as pd
import pymupdf
import yaml
from rdkit import Chem, RDLogger

from src.assemble import assemble_row, with_map_numbers

RDLogger.DisableLog("rdApp.*")

# Text as written in the tables -> SMILES fragment with one attachment point (*).
# Every entry is plain chemistry that can be checked by eye.
SUBSTITUENTS = {
    "-F": "*F",
    "-H": "*[H]",
    "-OCH3": "*OC",
    "-CH3": "*C",
    "-CH2-CH3": "*CC",
    "nC3H7": "*CCC",
    "cyclo-prop": "*C1CC1",
    "nC4H9": "*CCCC",
    "-CH2-iprop": "*CC(C)C",
    "nC5H11": "*CCCCC",
    "-(CH2)2-iprop": "*CCC(C)C",
    "nC6H13": "*CCCCCC",
    "c-hexyl": "*C1CCCCC1",
    "-phenyl": "*c1ccccc1",
    "-benzyl": "*Cc1ccccc1",
    "-ethylphenyl": "*CCc1ccccc1",        # read as 2-phenylethyl in the benzyl series; flagged below
    "-(CH2)2-CF3": "*CCC(F)(F)F",
    "-CH2-C(O)-CH3": "*CC(C)=O",
    "-CH2-C(OH)-CH3": "*CC(C)O",
    "-CH2-(o)F-phenyl": "*Cc1ccccc1F",
    "-CH2-(m)F-phenyl": "*Cc1cccc(F)c1",
    "-CH2-(p)F-phenyl": "*Cc1ccc(F)cc1",
}
AMBIGUOUS = {"-ethylphenyl"}  # text with more than one reasonable reading; reviewed by hand

DASHES = dict.fromkeys(map(ord, "\u2010\u2011\u2012\u2013\u2014\u2212"), "-")
NUMBER = re.compile(r"^([<>]?[\d.,]+|N/D\w*)$")   # numbers, or N/D with an optional footnote letter
RING_ATOMS = {
    "CH": "c", "N": "n",
    "-CH=": "c", "-N=": "n",
    "-C(Cl)=": "c(Cl)", "-C(CN)=": "c(C#N)",
}

# Linkers that join two parts of the scaffold, written as units in the drawing's reading order.
# Reversing the list gives the same linker read from the other side.
LINKERS = {
    "CO": ["C(=O)"],
    "CH2": ["C"],
    "S": ["S"],
    "SO": ["S(=O)"],
    "SO2": ["S(=O)(=O)"],
    "CH=CH": ["C=C"],                    # E/Z not specified in the table
    "(CO)-(CO)": ["C(=O)", "C(=O)"],
    "(COH)-(COH)": ["C(O)", "C(O)"],     # stereocenters not specified in the table
    "CONH": ["C(=O)", "N"],
    "NHCO": ["N", "C(=O)"],
}

def ring_dummy_isotopes(smiles):
    """Slot numbers of the wildcards that sit inside a ring (variable ring atoms)."""
    mol = Chem.MolFromSmiles(smiles)
    return sorted(a.GetIsotope() for a in mol.GetAtoms() if a.GetAtomicNum() == 0 and a.IsInRing())


def fill_ring_atoms(scaffold_smiles, replacements):
    """Swap each variable ring atom [k*] for an aromatic token like c, n, or c(Cl)."""
    for k, token in replacements.items():
        scaffold_smiles = scaffold_smiles.replace(f"[{k}*]", token, 1)
    return scaffold_smiles

def normalize(text):
    """Unify dash characters and remove spaces, so dictionary lookups are reliable."""
    return text.translate(DASHES).replace(" ", "")


def read_text_table(pdf_path, page_no, id_pattern, values_per_row=2, y_tol=3.0, n_parts=None):
    """Return rows of (compound, r_text, parts) from a text table, handling side by side halves.

    n_parts, if given, is the exact number of text columns a row must have, which keeps
    stray sentences in the body text from being read as table rows.
    """
    words = pymupdf.open(pdf_path)[page_no - 1].get_text("words")
    # Group words into visual lines: sort by vertical center, start a new line at a gap > y_tol
    entries = sorted(((y0 + y1) / 2, x0, w) for x0, y0, x1, y1, w, *_ in words)
    lines, current, last_y = [], [], None
    for y, x, w in entries:
        if last_y is not None and y - last_y > y_tol:
            lines.append(current)
            current = []
        current.append((x, w))
        last_y = y
    if current:
        lines.append(current)

    id_re, rows = re.compile(id_pattern), []
    for items in lines:
        tokens = [normalize(w) for _, w in sorted(items)]
        i = 0
        while i < len(tokens):
            if not id_re.match(tokens[i]):
                i += 1
                continue
            compound, i, r_parts = tokens[i], i + 1, []
            while i < len(tokens) and not NUMBER.match(tokens[i]):
                r_parts.append(tokens[i])
                i += 1
            values = tokens[i:i + values_per_row]
            if r_parts and len(values) == values_per_row and (n_parts is None or len(r_parts) == n_parts):
                rows.append({"compound": compound, "r_text": "".join(r_parts), "parts": "|".join(r_parts)})
            i += values_per_row
    return pd.DataFrame(rows)

def attach_single(scaffold_smiles, fragment_smiles):
    """Join a one-slot scaffold and a one-slot fragment at their attachment points."""
    if fragment_smiles in ("*[H]", "[H]*"):
        mol = Chem.RWMol(Chem.MolFromSmiles(scaffold_smiles))
        mol.RemoveAtom([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0][0])
        out = mol.GetMol()
    else:
        out = Chem.molzip(with_map_numbers(scaffold_smiles, 1), with_map_numbers(fragment_smiles, 1))
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out)


def assemble_text_table(results_dir, pdf_dir, paper_id, page, scaffold_crop, id_pattern):
    results_dir = Path(results_dir)
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    scaffold = validated.loc[validated.crop_id == scaffold_crop, "canonical_smiles"].iloc[0]

    table = read_text_table(Path(pdf_dir) / f"{paper_id}.pdf", page, id_pattern)
    records = []
    for r in table.itertuples():
        fragment = SUBSTITUENTS.get(r.r_text)
        smiles, status = None, "unmapped"
        if fragment:
            try:
                smiles, status = attach_single(scaffold, fragment), "assembled"
            except Exception as exc:
                status = f"error: {exc}"
        records.append({"paper_id": str(paper_id), "page": page, "compound": r.compound,
                        "r_text": r.r_text, "fragment": fragment, "assembled_smiles": smiles,
                        "ambiguous": r.r_text in AMBIGUOUS, "scaffold_crop": scaffold_crop,
                        "status": status})

    out = pd.DataFrame(records)
    path = results_dir / f"assembled_text_paper{paper_id}_page{page}.csv"
    out.to_csv(path, index=False)
    print(f"Read {len(out)} rows; assembled {(out.status == 'assembled').sum()} -> {path}")
    unmapped = out[out.status == "unmapped"]
    if len(unmapped):
        print("Unmapped R text (add to SUBSTITUENTS if correct):", unmapped.r_text.tolist())
    return out

def assemble_multi_slot(results_dir, pdf_dir, paper_id, page, scaffold_crop, id_pattern, columns):
    """Assemble rows with several text columns.

    columns has one entry per text column: "ring" for a variable ring atom (CH or N),
    or the scaffold slot number that column fills. The mapping must be checked per table,
    because the scaffold's slot numbers need not match the table's R labels.
    """
    results_dir = Path(results_dir)
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    scaffold = validated.loc[validated.crop_id == scaffold_crop, "canonical_smiles"].iloc[0]
    table = read_text_table(Path(pdf_dir) / f"{paper_id}.pdf", page, id_pattern, n_parts=len(columns))

    ring_isotopes = ring_dummy_isotopes(scaffold)
    records = []
    for r in table.itertuples():
        by_slot, ring_fill, unmapped = {}, {}, []
        for spec, text in zip(columns, r.parts.split("|")):
            if spec.startswith("ring"):
                # "ring:k" names the ring slot; plain "ring" means the scaffold's only ring slot
                k = int(spec.split(":")[1]) if ":" in spec else ring_isotopes[0]
                if text in RING_ATOMS:
                    ring_fill[k] = RING_ATOMS[text]
                else:
                    unmapped.append(text)
            else:
                fragment = SUBSTITUENTS.get(text)
                if fragment is None:
                    unmapped.append(text)
                by_slot[int(spec)] = fragment
        smiles, status = None, "unmapped: " + ",".join(unmapped) if unmapped else ""
        if not unmapped:
            try:
                filled = fill_ring_atoms(scaffold, ring_fill)
                smiles, status = assemble_row(filled, by_slot, "C"), "assembled"
            except Exception as exc:
                status = f"error: {exc}"
        records.append({"paper_id": str(paper_id), "page": page, "compound": r.compound,
                        "r_text": r.parts, "assembled_smiles": smiles,
                        "scaffold_crop": scaffold_crop, "status": status})
    out = pd.DataFrame(records)
    path = results_dir / f"assembled_text_paper{paper_id}_page{page}.csv"
    out.to_csv(path, index=False)
    print(f"Read {len(out)} rows; assembled {(out.status == 'assembled').sum()} -> {path}")
    return out

def assemble_linker_table(results_dir, pdf_dir, paper_id, page, scaffold_crop, id_pattern, reverse):
    """Insert a two point linker where the scaffold has a single wildcard between two parts.

    reverse=True when the SMILES string meets the linker from the opposite side of the drawing.
    """
    results_dir = Path(results_dir)
    validated = pd.read_csv(results_dir / "validated_predictions.csv", dtype={"paper_id": str})
    scaffold = validated.loc[validated.crop_id == scaffold_crop, "canonical_smiles"].iloc[0]
    if scaffold.count("*") != 1:
        raise ValueError("Linker scaffold must contain exactly one wildcard")
    table = read_text_table(Path(pdf_dir) / f"{paper_id}.pdf", page, id_pattern, n_parts=1)

    records = []
    for r in table.itertuples():
        units = LINKERS.get(r.r_text)
        smiles, status = None, "unmapped"
        if units:
            linker = "".join(reversed(units) if reverse else units)
            mol = Chem.MolFromSmiles(scaffold.replace("*", linker, 1))
            if mol is not None:
                smiles, status = Chem.MolToSmiles(mol), "assembled"
            else:
                status = "invalid after insertion"
        records.append({"paper_id": str(paper_id), "page": page, "compound": r.compound,
                        "r_text": r.r_text, "assembled_smiles": smiles,
                        "scaffold_crop": scaffold_crop, "reversed": reverse, "status": status})

    out = pd.DataFrame(records)
    path = results_dir / f"assembled_text_paper{paper_id}_page{page}.csv"
    out.to_csv(path, index=False)
    print(f"Read {len(out)} rows; assembled {(out.status == 'assembled').sum()} -> {path}")
    return out

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble molecules from a text R-group table.")
    parser.add_argument("--results", default="results/A")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--paper", default="19")
    parser.add_argument("--page", type=int, default=6)
    parser.add_argument("--scaffold-crop", default="19_p6_001")
    parser.add_argument("--id-pattern", default=r"^(1|12[a-z])$")
    parser.add_argument("--columns", help="Comma separated mapping per text column, e.g. ring,3,2,4,5")
    parser.add_argument("--linker", choices=["forward", "reverse"],
                        help="Treat R as a two point linker, inserted in drawing order or reversed")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    if args.linker:
        assemble_linker_table(args.results, cfg["data"]["pdf_dir"], args.paper, args.page,
                              args.scaffold_crop, args.id_pattern, args.linker == "reverse")
    elif args.columns:
        assemble_multi_slot(args.results, cfg["data"]["pdf_dir"], args.paper, args.page,
                            args.scaffold_crop, args.id_pattern, args.columns.split(","))
    else:
        assemble_text_table(args.results, cfg["data"]["pdf_dir"], args.paper, args.page,
                            args.scaffold_crop, args.id_pattern)