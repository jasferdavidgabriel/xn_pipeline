"""Expand text labels that recognition left as wildcards, using MolGlyph's own annotations.

MolGlyph's raw output ends with entries like <a>14:SO2C3H7</a>, meaning atom 14 is a
wildcard standing for the text SO2C3H7. Known abbreviations are replaced with their
structures; generic labels (like Heteroaryl/Aryl) stay as a plain *, because the answer
key stores generic compounds that way. Misread labels are matched to the closest known
label only when they are very similar, and every match is recorded for auditing.
"""
import argparse
import difflib
import re
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from src.validate import FRAGMENT_MAX_HEAVY_ATOMS

RDLogger.DisableLog("rdApp.*")

ABBREVIATIONS = {
    "SO2C3H7": "*S(=O)(=O)CCC",
    "NO2": "*[N+](=O)[O-]",
    "OMe": "*OC",
    "OCH3": "*OC",
    "Ph": "*c1ccccc1",
    "Bn": "*Cc1ccccc1",
    "Boc": "*C(=O)OC(C)(C)C",
    "SEM": "*COCC[Si](C)(C)C",
    "DMB": "*Cc1ccc(OC)cc1OC",
    "OBn": "*OCc1ccccc1",
    "BnO": "*OCc1ccccc1",
}
GENERIC = {"<dum>", "Heteroaryl/Aryl", "Heteroaryl", "Aryl", "Ar"}
FUZZY_CUTOFF = 0.75          # similarity needed to accept a misread label (e.g. SO2H7 -> SO2C3H7)
FUZZY_MIN_LENGTH = 4         # short labels change meaning with one letter (OBn vs Bn), so never fuzzy match them
ANNOTATION = re.compile(r"<a>(\d+):(.*?)</a>")


def resolve_label(label):
    """Return (kind, fragment, matched_label), where kind is abbreviation, generic, or unknown."""
    if label in GENERIC:
        return "generic", None, label
    if label in ABBREVIATIONS:
        return "abbreviation", ABBREVIATIONS[label], label
    candidates = list(ABBREVIATIONS) + sorted(GENERIC - {"<dum>"})
    close = []
    if len(label) >= FUZZY_MIN_LENGTH:
        close = difflib.get_close_matches(label, candidates, n=1, cutoff=FUZZY_CUTOFF)
    if close:
        match = close[0]
        if match in GENERIC:
            return "generic", None, match
        return "abbreviation", ABBREVIATIONS[match], match
    return "unknown", None, label


def expand(raw_output, extra=None, dum_fragment=None):
    """Return (smiles, status, notes) after expanding annotated wildcards.
    extra maps labels to fragments for this call only (e.g. {"C": linker C from a scheme's box});
    dum_fragment, if given, fills an unlabeled <dum> wildcard.
    """
    if not isinstance(raw_output, str) or "<sep>" not in raw_output:
        return None, "no annotations", ""
    base, extension = raw_output.split("<sep>", 1)
    annotations = ANNOTATION.findall(extension)
    if not annotations:
        return None, "no annotations", ""
    mol = Chem.MolFromSmiles(base)
    if mol is None:
        return None, "unparseable", ""

    rw, joins, notes, n_generic = Chem.RWMol(mol), [], [], 0
    for idx, label in annotations:
        idx = int(idx)
        if idx >= rw.GetNumAtoms() or rw.GetAtomWithIdx(idx).GetAtomicNum() != 0:
            notes.append(f"{label}: not on a wildcard")
            continue
        if extra and label in extra:
            kind, fragment, matched = "abbreviation", extra[label], label
        elif label == "<dum>" and dum_fragment:
            kind, fragment, matched = "abbreviation", dum_fragment, "blue placeholder"
        else:
            kind, fragment, matched = resolve_label(label)
        notes.append(f"{label} -> {matched}" if kind != "unknown" else f"{label} -> ?")
        if kind == "unknown":
            return None, "unknown label", "; ".join(notes)
        if kind == "generic":
            n_generic += 1
        else:
            map_no = 100 + len(joins)
            rw.GetAtomWithIdx(idx).SetAtomMapNum(map_no)
            joins.append((map_no, fragment))

    out = rw.GetMol()
    for map_no, fragment in joins:
        piece = Chem.MolFromSmiles(fragment)
        for atom in piece.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(map_no)
        out = Chem.molzip(out, piece)
    Chem.SanitizeMol(out)

    remaining = [a for a in out.GetAtoms() if a.GetAtomicNum() == 0]
    if len(remaining) != n_generic:      # a wildcard without a generic label is still incomplete
        return None, "unresolved wildcard", "; ".join(notes)
    for atom in remaining:
        atom.SetIsotope(0)
        atom.SetAtomMapNum(0)
    if not remaining:
        return Chem.MolToSmiles(out), "expanded_complete", "; ".join(notes)
    # A generic compound is a full-size structure with one generic attachment point.
    # Small pieces (table fragments) and two-point pieces (linkers) are not compounds.
    heavy = sum(1 for a in out.GetAtoms() if a.GetAtomicNum() > 1)
    if len(remaining) == 1 and heavy > FRAGMENT_MAX_HEAVY_ATOMS:
        return Chem.MolToSmiles(out), "expanded_generic", "; ".join(notes)
    return Chem.MolToSmiles(out), "generic_piece", "; ".join(notes)

def expand_all(results_dir):
    """Expand annotated wildcards for every crop and save expanded_predictions.csv."""
    results_dir = Path(results_dir)
    raw = pd.read_csv(results_dir / "raw_predictions.csv", dtype={"paper_id": str})
    rows = []
    for r in raw.itertuples():
        smiles, status, notes = expand(r.raw_output)
        rows.append({"crop_id": r.crop_id, "paper_id": r.paper_id,
                     "expanded_smiles": smiles, "status": status, "labels": notes})
    out = pd.DataFrame(rows)
    out.to_csv(results_dir / "expanded_predictions.csv", index=False)
    print(out.status.value_counts().to_string())
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Expand abbreviation labels left as wildcards.")
    parser.add_argument("--results", default="results/A")
    args = parser.parse_args()
    expand_all(args.results)