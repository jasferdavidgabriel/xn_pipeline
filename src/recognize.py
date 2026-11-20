"""Stage 3: recognize molecule images with MolGlyph and convert its output to SMILES.

Mirrors the Phase 1 baseline: MolGlyph from BioMiner's MolScribe package, plus BioMiner's
conversion utilities. Paths and device come from config.yaml, with command line overrides
so the same code runs on a Mac (cpu) or Colab (cuda).
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import pandas as pd
import yaml


def load_molglyph(biominer_dir, weights_path, device):
    """Load the MolGlyph model and BioMiner's SMILES conversion module."""
    import torch

    biominer_dir = Path(biominer_dir)
    sys.path.insert(0, str(biominer_dir / "BioMiner" / "MolScribe"))
    from molscribe import MolGlyph

    model = MolGlyph(str(weights_path), torch.device(device))

    utils_path = biominer_dir / "BioMiner" / "commons" / "molglyph_utils.py"
    spec = importlib.util.spec_from_file_location("molglyph_conversion", utils_path)
    conversion = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = conversion
    spec.loader.exec_module(conversion)
    return model, conversion


def convert_prediction(conversion, raw_output):
    """Turn MolGlyph's raw output into SMILES; None when it flags unresolved parts (<r> or <c>)."""
    if not isinstance(raw_output, str) or "<sep>" not in raw_output:
        raise ValueError("Prediction is missing the <sep> separator.")
    base_smiles, extension = raw_output.split("<sep>", 1)
    if "<r>" in extension or "<c>" in extension:
        return None
    refactored = conversion.get_refactor(conversion.Translator.refactor(raw_output))
    return conversion.get_new_smiles(refactored, base_smiles, extension)


def recognize_images(table, id_col, path_col, model, conversion):
    """Run recognition on every image listed in table; never stop on a single failure."""
    records = []
    for number, row in enumerate(table.to_dict("records"), start=1):
        started = time.perf_counter()
        raw, smiles, error = None, None, ""
        try:
            raw = model.predict_image_file(str(row[path_col]), return_confidence=False)
            smiles = convert_prediction(conversion, raw)
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
        records.append({
            "crop_id": row[id_col],
            "paper_id": str(row["paper_id"]),
            "raw_smiles": smiles,
            "raw_output": raw,
            "inference_time": round(time.perf_counter() - started, 3),
            "error": error,
        })
        print(f"{number}/{len(table)} {row[id_col]}: {smiles if smiles else error or 'unresolved'}", flush=True)
    return pd.DataFrame(records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recognize molecule images with MolGlyph.")
    parser.add_argument("--input", required=True, help="CSV listing images to recognize")
    parser.add_argument("--output", required=True, help="Where to write raw predictions")
    parser.add_argument("--id-col", default="crop_id", help="Column with each image's ID")
    parser.add_argument("--path-col", default="crop_path", help="Column with each image's path")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--biominer-dir", help="Override config biominer.repo_dir")
    parser.add_argument("--weights", help="Override config biominer.molglyph_weights")
    parser.add_argument("--device", help="cpu or cuda; overrides config")
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config)).get("biominer", {})
    biominer_dir = args.biominer_dir or cfg.get("repo_dir", "../BioMiner-main")
    weights = args.weights or cfg.get("molglyph_weights", f"{biominer_dir}/molglyph_large.pt")
    device = args.device or cfg.get("device", "cpu")

    print(f"Loading MolGlyph on {device}...", flush=True)
    model, conversion = load_molglyph(biominer_dir, weights, device)
    table = pd.read_csv(args.input, dtype={"paper_id": str})
    results = recognize_images(table, args.id_col, args.path_col, model, conversion)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, index=False)
    print(f"Saved {len(results)} predictions to {args.output}")