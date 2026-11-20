"""Stage 4: canonicalize raw SMILES with RDKit and flag each prediction's status."""
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")  # hide RDKit warnings for unparseable SMILES
# Substituent pieces are small; above this size, a single-star structure is treated as a scaffold
FRAGMENT_MAX_HEAVY_ATOMS = 20


def classify_smiles(raw_smiles):
    """Return (canonical_smiles, status).

    Status is one of:
      valid     complete molecule
      invalid   empty or unparseable
      scaffold  core with R-group slots: numbered stars, one star on a large molecule,
                or several stars on an aromatic core
      fragment  small substituent piece with a single attachment point, like *N1CCOCC1
      linker    chain with an attachment point at each end and no aromatic rings
    """
    if not isinstance(raw_smiles, str) or not raw_smiles.strip():
        return None, "invalid"
    mol = Chem.MolFromSmiles(raw_smiles)
    if mol is None:
        return None, "invalid"
    canonical = Chem.MolToSmiles(mol)

    # Wildcard atoms have atomic number 0; R-group numbers are stored as the atom's isotope
    dummies = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if not dummies:
        return canonical, "valid"

    real_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    has_aromatic = any(atom.GetIsAromatic() for atom in real_atoms)
    numbered = any(atom.GetIsotope() > 0 for atom in dummies)

    if numbered:
        return canonical, "scaffold"
    if len(dummies) == 1:
        if len(real_atoms) <= FRAGMENT_MAX_HEAVY_ATOMS:
            return canonical, "fragment"
        return canonical, "scaffold"          # large molecule with a single R
    if len(dummies) == 2 and not has_aromatic:
        return canonical, "linker"
    return canonical, "scaffold"

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