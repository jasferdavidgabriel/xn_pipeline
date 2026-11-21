"""Rebuild every result from the committed inputs, in pipeline order. Runs on a laptop CPU.

GPU stages (MolGlyph recognition of missed table cells, targeted OCR) run on Colab; their
outputs are committed under results/A and used here as inputs. See colab_setup.md.

Usage:
    python run_all.py
"""
import argparse
from pathlib import Path

import yaml
import pandas as pd

from src.assemble import assemble_table
from src.expand_abbreviations import expand_all
from src.score import score
from src.text_tables import assemble_linker_table, assemble_multi_slot, assemble_text_table
from src.validate import validate_predictions
from src.pieces import linker_table, warhead_table
from src.image_tables import assemble_image_table
from src.grid_tables import assemble_grid_table
from src.scheme_placeholders import resolve_scheme


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

    print("\n== 3b. Assemble modular pieces (warheads, linkers, ligands)")
    smiles_of = lambda crop: pd.read_csv(results / "validated_predictions.csv").set_index("crop_id").canonical_smiles[crop]
    for job in cfg["tables"].get("pieces_warhead", []):
        t = warhead_table(results, pdf_dir, job["paper"], job["page"], job["piece_crop"], job["name_pattern"])
        t.to_csv(results / f"assembled_pieces_paper{job['paper']}_page{job['page']}.csv", index=False)
        print(f"Warhead table page {job['page']}: {len(t)} compounds")
    for job in cfg["tables"].get("pieces_linker", []):
        t = linker_table(results, pdf_dir, job["paper"], job["page"], smiles_of(job["warhead_crop"]),
                         smiles_of(job["ligand_piece_crop"]), job["name_pattern"])
        t.to_csv(results / f"assembled_pieces_paper{job['paper']}_page{job['page']}.csv", index=False)
        print(f"Linker table page {job['page']}: {len(t)} compounds")

    print("\n== 3c. Assemble OCR-read image tables")
    for job in cfg["tables"].get("image", []):
        assemble_image_table(results, job["paper"], job["page"], job["scaffolds"],
                             job["ring_slot"], job["r_slot"])

    print("\n== 3d. Assemble grid-detected tables")
    for job in cfg["tables"].get("grid", []):
        assemble_grid_table(results, job["paper"], job["page"], job["scaffold_crop"])

    print("\n== 3e. Resolve scheme placeholders")
    for job in cfg["tables"].get("scheme", []):
        resolve_scheme(results, pdf_dir, job["paper"], job["page"], job["box_crops"],
                       letters=job["letters"], blue_letter=job["blue_letter"])

    print("\n== 4. Expand abbreviation labels")
    expand_all(results)

    print("\n== 5. Score")
    score(results, cfg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rebuild all pipeline results.")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--results", default="results/A300")
    args = parser.parse_args()
    main(args.config, args.results)