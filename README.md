# XN Pipeline: Extracting Molecules from Scientific PDFs

Phase 2 pipeline of XN Project Group 5 (Computer Vision, Northeastern University). It turns molecule drawings and R-group tables in medicinal chemistry papers into validated SMILES, and scores them against a curated answer key of 220 molecules from 6 papers.

The pipeline builds on the team's Phase 1 benchmark of BioMiner's components (MolDetv2 detection and MolGlyph recognition). Phase 2 keeps those models and adds the stages that were missing: validation, R-group assembly, table reading, abbreviation expansion, and assembly from separately drawn pieces.

## Results

| Paper | Molecules | Phase 1 baseline | Phase 2 pipeline |
|---|---|---|---|
| 1 | 17 | 100.0% | 100.0% |
| 2 | 34 | 26.5% | **88.2%** |
| 4 | 60 | 23.3% | **53.3%** |
| 6 | 1 | 0.0% | 0.0% |
| 16 | 4 | 25.0% | 25.0% |
| 19 | 104 | 18.3% | **65.4%** |
| **All** | **220** | **27.3%** | **67.3%** |
| Precision | | 27.2% | **88.0%** |

Metric: exact structure recovery against 220 unique reference structures (duplicates within a paper removed), after RDKit canonicalization. Tanimoto similarity uses Morgan fingerprints (radius 2, 2048 bits) with chirality enabled.

### Where the gains came from

| Stage added | Overall recall | Precision |
|---|---|---|
| Baseline (detection + recognition, 200 DPI) | 27.3% | 27.2% |
| + Paper 2 table assembly, filled cells, targeted OCR | 36.8% | ~83% |
| + Paper 19 text tables (Tables 1, 3, 6, 9) | 56.8% | 87.3% |
| + Abbreviation expansion | 60.0% | 84.9% |
| + Recognition at 300 DPI | 60.5% | 86.8% |
| + Paper 4 assembly from drawn pieces (Tables 1, 3) | 67.3% | 88.0% |

### Key findings

1. **Recognition is not the main bottleneck.** When recognition sees a complete molecule, it is right about 73% of the time, matching its clean-image benchmark (72%). Most end-to-end loss came from molecules defined as **scaffold + substituent table**, or as **separately drawn pieces**, which detection and recognition alone cannot resolve.
2. **Assembling from pieces can beat reading the whole molecule.** CCT400028's full drawing was read with one carbon missing from a long chain. The same chain, drawn on its own in another table, was read correctly, and assembling from pieces recovered the exact molecule.
3. **A Tanimoto of 1.0 does not guarantee the same molecule.** Radius 2 Morgan fingerprints cannot tell a chain of four CH2 groups from a chain of five, so a molecule missing one chain carbon scored 1.0. This is why exact match is the headline metric.
4. **Resolution is not the main lever for recognition errors.** A controlled experiment (same crop boxes and model, 200 vs 300 DPI) gained 1 molecule and lost none, but did not fix miscounted chains or misplaced groups, and made some readings worse. 300 DPI is used as the default because it was strictly better on exact match.

## Quick start

```bash
conda create -n biobench -c conda-forge --override-channels python=3.10 -y
conda activate biobench
pip install -r requirements.txt
python run_all.py
```

`run_all.py` rebuilds every result from the committed inputs on a laptop CPU, in about a minute. Rebuilding is deterministic: the output files are byte for byte identical to the committed ones. Results go to `results/A300/` by default.

**Data:** place the PDFs in `data/pdfs/` (named `<paperID>.pdf`) and the answer key at `data/ground_truth.csv`. Both are excluded from Git.

**GPU stages** (MolGlyph recognition, targeted OCR) run on Google Colab. Their outputs are committed, so `run_all.py` does not need a GPU. To rerun them, see [colab_setup.md](colab_setup.md).

## Pipeline stages

```
PDF pages
  └─ Detection (MolDetv2) ─> crops.csv
       └─ Recognition (MolGlyph, 300 DPI) ─> raw_predictions.csv
            └─ Validation (RDKit) ─> validated_predictions.csv
                 ├─ valid: complete molecules ─────────────────────────────────┐
                 ├─ scaffold + drawn R-group table ─> assembly ────────────────┤
                 ├─ scaffold + text R-group table ─> assembly ─────────────────┤
                 ├─ warheads + linkers + ligand pieces ─> modular assembly ────┤
                 └─ annotated wildcards ─> abbreviation expansion ─────────────┤
                                                                               └─> Scoring
```

| Stage | File | What it does |
|---|---|---|
| Import baseline | `src/import_baseline.py` | Converts the Phase 1 run into the shared CSV interfaces |
| Rebuild crops | `src/rebuild_crops.py` | Re-creates crop images from PDFs and saved boxes, at any DPI |
| Recognition | `src/recognize.py` | MolGlyph on any list of images (Colab) |
| Validation | `src/validate.py` | Canonical SMILES plus a status: valid, invalid, scaffold, fragment, linker |
| Drawn table assembly | `src/assemble.py` | Groups fragment crops into table rows and columns by position, then joins them to the scaffold |
| Missed cells | `src/fill_cells.py` | Crops table cells that detection missed, for recognition |
| Table text | `src/read_table_text.py` | Targeted OCR of table text in images; dash detection by ink height (Colab) |
| Text table assembly | `src/text_tables.py` | Reads text tables from the PDF text layer; maps substituents, ring atoms, and linkers to structures |
| Modular assembly | `src/pieces.py` | Joins separately drawn warheads, linkers, and ligand pieces |
| Abbreviations | `src/expand_abbreviations.py` | Expands labels MolGlyph left as wildcards (e.g. SO2C3H7), using its own annotations |
| Scoring | `src/score.py`, `src/evaluate.py` | Final predictions, per paper and per molecule results |
| Review | `src/review_crops.py` | Visual review sheets of crops by paper and status |
| All stages | `run_all.py` | Rebuilds every result in order |

Paper specific settings (table pages, scaffold and piece crops, slot mappings, linker direction, name patterns) live in `config.yaml`, not in code.

## Design decisions and known limitations

1. **Slot mappings are verified per table.** MolGlyph's wildcard numbers do not always match a table's R labels (in paper 19 Table 9, Ra maps to slot 3). Every mapping in `config.yaml` was confirmed by a self test (rebuilding vemurafenib) and by name level scoring.
2. **Linker direction was tested, not assumed.** In paper 19 Table 3, both directions were run; the answer key confirmed the reversed direction for the two directional linkers (CONH, NHCO). In paper 4, two-ended linkers are oriented by a chemical rule: the end bonded to N joins the ligand's carbonyl. All 7 linkers confirmed it.
3. **Pieces may be reused across tables.** Paper 4 Table 1's scaffold failed recognition, so its warhead and ligand pieces come from Table 3's drawings of the same structures. Compound 44, which appears in both tables, was built identically by both routes.
4. **Some names come from data tables, not figure labels.** Where figure labels are images (paper 4 Table 3), names are read from the data table below, whose rows follow the figure's left to right order. All 10 compounds matched by name, verifying the order.
5. **Generic compounds are kept as wildcards** (e.g. "Heteroaryl/Aryl"), because the answer key stores them that way. Only full size structures with one generic attachment point count; table fragments and linkers do not.
6. **Fuzzy label matching is limited and audited.** Misread labels (SO2H7, Bryl) are matched to known labels only at 75% similarity or more, and only for labels of 4 or more characters. Every match is recorded in `expanded_predictions.csv`.
7. **Paper 2 compound numbers are inferred from row order**, confirmed by OCR reading compounds 10 to 21 in the matching rows.
8. **Open decision: compound 10e** (paper 19). The table writes `CH=CH` without E/Z geometry; the answer key uses E. Left as a known miss pending a documented default.
9. **Remaining recognition errors** (a miscounted chain carbon, a misplaced methyl, an extra methyl, propyl read as ethyl) cannot be fixed by assembly or resolution.
10. **Not yet handled:** paper 4 Scheme 1 intermediates (blue "C" placeholders), paper 4's named literature PROTACs, and paper 19 Tables 7 and 11 (image tables, merged cells, positional isomers).

## Repository layout

```
config.yaml          data paths, BioMiner paths, table settings
run_all.py           rebuild all results
src/                 pipeline stages
results/A/           200 DPI run (Phase 1 baseline crops + Phase 2 stages)
results/A300/        300 DPI run (current default)
data/                PDFs and answer key (not in Git)
colab_setup.md       GPU stages on Colab
CONTRIBUTIONS.md     who did what
```

## Credits

See [CONTRIBUTIONS.md](CONTRIBUTIONS.md). Models: MolDetv2 (UniParser) and MolGlyph (BioMiner, Yan et al., 2025); chemistry with RDKit; OCR with EasyOCR.
