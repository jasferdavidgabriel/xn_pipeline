"""Classify every missed molecule by how its closest prediction differs from the reference."""
import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


def formula_counts(mol):
    """Element counts, with hydrogens included."""
    counts = {}
    for atom in Chem.AddHs(mol).GetAtoms():
        counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
    return counts


def classify(ref_smiles, pred_smiles):
    if not isinstance(pred_smiles, str):
        return "no prediction", ""
    ref, pred = Chem.MolFromSmiles(ref_smiles), Chem.MolFromSmiles(pred_smiles)
    if ref is None or pred is None:
        return "unparseable", ""
    if Chem.MolToSmiles(ref, isomericSmiles=False) == Chem.MolToSmiles(pred, isomericSmiles=False):
        return "stereo only", "same atoms and bonds"
    fr, fp = formula_counts(ref), formula_counts(pred)
    diff = {el: fp.get(el, 0) - fr.get(el, 0) for el in set(fr) | set(fp)}
    diff = {el: d for el, d in diff.items() if d}
    text = " ".join(f"{'+' if d > 0 else ''}{d}{el}" for el, d in sorted(diff.items()))
    if not diff:
        return "isomer", "same formula, different connections"
    if diff in ({"C": 1, "H": 2}, {"C": -1, "H": -2}):
        return "chain length (CH2)", text
    if sum(abs(d) for el, d in diff.items() if el != "H") <= 2:
        return "small formula error", text
    return "large difference", text


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Classify missed molecules by error type.")
    parser.add_argument("--results", default="results/A300")
    args = parser.parse_args()
    m = pd.read_csv(Path(args.results) / "final_per_molecule.csv", dtype={"paper_id": str})
    miss = m[~m.exact_match].copy()
    miss[["category", "detail"]] = [classify(r, p) for r, p in zip(miss.true_canonical, miss.best_prediction)]
    miss.to_csv(Path(args.results) / "error_analysis.csv", index=False)
    pd.set_option("display.width", 200)
    print(miss[["paper_id", "molecule_name", "best_tanimoto", "category", "detail"]]
          .sort_values(["category", "paper_id"]).to_string(index=False))
    print()
    print(miss.category.value_counts().to_string())