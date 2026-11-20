"""Stage 4: canonicalize raw SMILES with RDKit and flag each prediction's status."""
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")  # hide RDKit warnings for unparseable SMILES


def classify_smiles(raw_smiles):
    """Return (canonical_smiles, status), where status is valid, invalid, or wildcard."""
    if not isinstance(raw_smiles, str) or not raw_smiles.strip():
        return None, "invalid"                    # empty or missing prediction
    mol = Chem.MolFromSmiles(raw_smiles)
    if mol is None:
        return None, "invalid"                       # Blank 1: RDKit couldn't parse it
    canonical = Chem.MolToSmiles(mol)
    # Wildcard atoms ("*") stand for unresolved R-groups
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):   # Blank 2
        return canonical, "wildcard"
    return canonical, "valid"


def validate_predictions(raw_csv, output_csv):
    """Read raw_predictions.csv, add canonical SMILES and status, save validated_predictions.csv."""
    raw = pd.read_csv(raw_csv)
    results = raw["raw_smiles"].apply(classify_smiles)
    validated = raw[["crop_id", "paper_id"]].copy()
    validated["canonical_smiles"] = [canonical for canonical, _ in results]
    validated["status"] = [status for _, status in results]
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    validated.to_csv(output_csv, index=False)
    print(validated["status"].value_counts())
    return validated