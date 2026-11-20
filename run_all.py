"""Rebuild every result from the committed inputs, in pipeline order. Runs on a laptop CPU.

GPU stages (MolGlyph recognition of missed table cells, targeted OCR) run on Colab; their
outputs are committed under results/A and used here as inputs. See colab_setup.md.

Usage:
    python run_all.py
"""
import argparse
from pathlib import Path

import yaml

from src.assemble import assemble_table
from src.expand_abbreviations import expand_all
from src.score import score
from src.text_tables import assemble_linker_table, assemble_multi_slot, assemble_text_table
from src.validate import validate_predictions


def main(config_path, results_dir):
    cfg = yaml.safe_load(open(config_path))
    pdf_dir, results = cfg["data"]["pdf_dir"], Path(results_dir)

    # Remove old assembly outputs, so a table dropped from the config can't linger in the score
    for old in results.glob("assembled_*.csv"):
        old.unlink()

    print("== 1. Validate recognition output")
    validate_predictions(results / "raw_predictions.csv", results / "validated_predictions.csv")

    print("\n== 2. Assemble drawn R-group tables")
    for job in cfg["tables"].get("drawn", []):
        assemble_table(results, pdf_dir, job["paper"], job["page"])

    print("\n== 3. Assemble text R-group tables")
    for job in cfg["tables"].get("text", []):
        args = (results, pdf_dir, job["paper"], job["page"], job["scaffold_crop"], job["id_pattern"])
        if "linker" in job:
            assemble_linker_table(*args, job["linker"] == "reverse")
        elif "columns" in job:
            assemble_multi_slot(*args, [str(c) for c in job["columns"]])
        else:
            assemble_text_table(*args)

    print("\n== 4. Expand abbreviation labels")
    expand_all(results)

    print("\n== 5. Score")
    score(results, cfg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rebuild all pipeline results.")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--results", default="results/A")
    args = parser.parse_args()
    main(args.config, args.results)