# XN Pipeline: Extracting Molecules from Scientific PDFs

Phase 2 pipeline of XN Project Group 5 (Computer Vision, Northeastern University). It turns molecule drawings and R-group tables in medicinal chemistry papers into validated SMILES, and scores them against a curated answer key of 220 molecules from 6 papers.

The pipeline builds on the team's Phase 1 benchmark of BioMiner's components (MolDetv2 detection and MolGlyph recognition). Phase 2 keeps those models and adds the stages that were missing: validation, R-group assembly, table reading, and abbreviation expansion.

## Results

| Paper | Molecules | Phase 1 baseline | Phase 2 pipeline |
|---|---|---|---|
| 1 | 17 | 100.0% | 100.0% |
| 2 | 34 | 26.5% | **88.2%** |
| 4 | 60 | 23.3% | 28.3% |
| 6 | 1 | 0.0% | 0.0% |
| 16 | 4 | 25.0% | 25.0% |
| 19 | 104 | 18.3% | **64.4%** |
| **All** | **220** | **27.3%** | **60.0%** |
| Precision | | 27.2% | **84.9%** |

Metric: exact structure recovery against 220 unique reference structures (duplicates within a paper removed), after RDKit canonicalization. Tanimoto similarity uses Morgan fingerprints (radius 2, 2048 bits) with chirality enabled.

### Where the gains came from

| Stage added | Overall recall | Precision |
|---|---|---|
| Baseline (detection + recognition) | 27.3% | 27.2% |
| + Paper 2 table assembly, filled cells, targeted OCR | 36.8% | ~83% |
| + Paper 19 text tables (Tables 1, 3, 6, 9) | 56.8% | 87.3% |
| + Abbreviation expansion | 60.0% | 84.9% |

The key finding: when recognition sees a complete molecule, it is right about 73% of the time, matching its clean-image benchmark (72%). Most of the end-to-end loss came from molecules that exist only as **scaffold + substituent table**, which detection and recognition alone cannot resolve.

## Quick start

```bash
conda create -n biobench -c conda-forge --override-channels python=3.10 -y
conda activate biobench
pip install -r requirements.txt
python run_all.py
```

`run_all.py` rebuilds every result from the committed inputs on a laptop CPU, in about a minute. Rebuilding is deterministic: the output files are byte for byte identical to the committed ones.

**Data:** place the PDFs in `data/pdfs/` (named `<paperID>.pdf`) and the answer key at `data/ground_truth.csv`. Both are excluded from Git.

**GPU stages** (recognizing missed table cells with MolGlyph, targeted OCR) run on Google Colab. Their outputs are committed under `results/A/`, so `run_all.py` does not need a GPU. To rerun them, see [colab_setup.md](colab_setup.md).

## Pipeline stages

```
PDF pages
  └─ Detection (MolDetv2) ─> crops.csv
       └─ Recognition (MolGlyph) ─> raw_predictions.csv
            └─ Validation (RDKit) ─> validated_predictions.csv
                 ├─ valid: complete molecules ──────────────────────────────┐
                 ├─ scaffold + fragment: drawn R-group tables ─> assembly ──┤
                 ├─ scaffold + text table: text R-group tables ─> assembly ─┤
                 └─ annotated wildcards ─> abbreviation expansion ──────────┤
                                                                            └─> Scoring
```

| Stage | File | What it does |
|---|---|---|
| Import baseline | `src/import_baseline.py` | Converts the Phase 1 run into the shared CSV interfaces |
| Rebuild crops | `src/rebuild_crops.py` | Re-creates crop images from PDFs and saved boxes |
| Recognition | `src/recognize.py` | MolGlyph on any list of images (Colab) |
| Validation | `src/validate.py` | Canonical SMILES plus a status: valid, invalid, scaffold, fragment, linker |
| Drawn table assembly | `src/assemble.py` | Groups fragment crops into table rows and columns by position, then joins them to the scaffold |
| Missed cells | `src/fill_cells.py` | Crops table cells that detection missed, for recognition |
| Table text | `src/read_table_text.py` | Targeted OCR of table text in images; dash detection by ink height (Colab) |
| Text table assembly | `src/text_tables.py` | Reads text tables from the PDF text layer and maps R text to structures |
| Abbreviations | `src/expand_abbreviations.py` | Expands labels MolGlyph left as wildcards (e.g. SO2C3H7), using its own annotations |
| Scoring | `src/score.py`, `src/evaluate.py` | Final predictions, per paper and per molecule results |
| Review | `src/review_crops.py` | Visual review sheets of crops by paper and status |
| All stages | `run_all.py` | Rebuilds every result in order |

Paper specific settings (table pages, scaffold crops, slot mappings, linker direction) live in `config.yaml`, not in code.

## Design decisions and known limitations

1. **Slot mappings are verified per table.** MolGlyph's wildcard numbers do not always match a table's R labels (in paper 19 Table 9, Ra maps to slot 3). Every mapping in `config.yaml` was confirmed by a self test (rebuilding vemurafenib) and by name level scoring.
2. **Linker direction was tested, not assumed.** For paper 19 Table 3, both directions were run; the answer key confirmed the reversed direction for the two directional linkers (CONH, NHCO).
3. **Generic compounds are kept as wildcards** (e.g. "Heteroaryl/Aryl"), because the answer key stores them that way. Only full size structures with one generic attachment point count; table fragments and linkers do not.
4. **Fuzzy label matching is limited and audited.** Misread labels (SO2H7, Bryl) are matched to known labels only at 75% similarity or more, and only for labels of 4 or more characters. Every match is recorded in `expanded_predictions.csv`.
5. **Paper 2 compound numbers are inferred from row order**, confirmed by OCR reading compounds 10 to 21 in the matching rows.
6. **Open decision: compound 10e** (paper 19). The table writes `CH=CH` without E/Z geometry; the answer key uses E. Left as a known miss pending a documented default.
7. **Remaining recognition errors:** B9 (propyl read as ethyl) and the Scheme 2 product (an extra imine) are recognition mistakes no assembly step can fix.
8. **Paper 4** (PROTACs built from interchangeable modules) and **paper 19 Tables 7 and 11** (image tables, merged cells, positional isomers) are not yet handled.

## Repository layout

```
config.yaml          data paths, BioMiner paths, table settings
run_all.py           rebuild all results
src/                 pipeline stages
results/A/           Config A outputs (full page detection baseline + Phase 2 stages)
data/                PDFs and answer key (not in Git)
colab_setup.md       GPU stages on Colab
CONTRIBUTIONS.md     who did what
```

## Credits

See [CONTRIBUTIONS.md](CONTRIBUTIONS.md). Models: MolDetv2 (UniParser) and MolGlyph (BioMiner, Yan et al., 2025); chemistry with RDKit; OCR with EasyOCR.
