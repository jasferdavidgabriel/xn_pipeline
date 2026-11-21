# XN Pipeline: Extracting Molecules from Scientific PDFs

Phase 2 pipeline of XN Project Group 5 (Computer Vision, Northeastern University). It turns molecule drawings, R-group tables, and synthesis schemes in medicinal chemistry papers into validated SMILES, and scores them against a curated answer key of 220 molecules from 6 papers.

The pipeline builds on the team's Phase 1 benchmark of BioMiner's components (MolDetv2 detection and MolGlyph recognition). Phase 2 keeps those models and adds the stages that were missing: validation, R-group assembly, table reading (text, image, and grid), assembly from separately drawn pieces, abbreviation expansion, and scheme placeholder resolution.

## Results

| Paper | Molecules | Phase 1 baseline | Phase 2 pipeline |
|---|---|---|---|
| 1 | 17 | 100.0% | 100.0% |
| 2 | 34 | 26.5% | **88.2%** |
| 4 | 60 | 23.3% | **80.0%** |
| 6 | 1 | 0.0% | 0.0% |
| 16 | 4 | 25.0% | 25.0% |
| 19 | 104 | 18.3% | **95.2%** |
| **All** | **220** | **27.3%** | **88.6%** |
| Precision | | 27.2% | **90.6%** |

Metric: exact structure recovery against 220 unique reference structures (duplicates within a paper removed), after RDKit canonicalization. Tanimoto similarity uses Morgan fingerprints (radius 2, 2048 bits) with chirality enabled.

### Where the gains came from

| Stage added | Overall recall | Precision |
|---|---|---|
| Baseline (detection + recognition, 200 DPI) | 27.3% | 27.2% |
| + Paper 2 drawn R-group table, filled cells, targeted OCR | 36.8% | ~83% |
| + Paper 19 text tables (Tables 1, 3, 6, 9) | 56.8% | 87.3% |
| + Abbreviation expansion | 60.0% | 84.9% |
| + Recognition at 300 DPI | 60.5% | 86.8% |
| + Paper 4 assembly from drawn pieces (Tables 1, 3) | 67.3% | 88.0% |
| + Paper 19 Table 7: OCR combined with table logic | 70.5% | 88.4% |
| + Paper 19 Table 11: grid detection and merged cells | 81.4% | 89.8% |
| + Paper 4 Scheme 1: placeholders from the definition box and color | 85.0% | 90.3% |
| + Paper 4 Scheme 1: compounds defined by text labels | 88.6% | 90.6% |

Precision rose as recall climbed: each stage added real molecules, not noise.

### Key findings

1. **Recognition is not the main bottleneck.** When recognition sees a complete molecule, it is right about 73% of the time, matching its clean-image benchmark (72%). Most end-to-end loss came from molecules defined indirectly: scaffold + substituent table, separately drawn pieces, or placeholders in synthesis schemes.
2. **Assembling from pieces can beat reading the whole molecule.** CCT400028's full drawing was read with one carbon missing from a long chain. The same chain, drawn on its own in another table, was read correctly, and assembling from pieces recovered the exact molecule.
3. **A Tanimoto of 1.0 does not guarantee the same molecule.** Radius 2 Morgan fingerprints cannot tell a chain of four CH2 groups from a chain of five. This is why exact match is the headline metric.
4. **OCR alone is unreliable on tiny table cells; OCR plus the table's own logic is not.** In paper 19 Table 7, raw OCR swapped or missed the scaffold letters and added stray digits to entry numbers. Combining it with the table's structure (consecutive entries, "a dash means scaffold B") gave 7 of 7 correct. Snapping garbled text to a small vocabulary did the same for Table 11's R1 column.
5. **A positional descriptor is meaningless without its numbering convention.** Table 11's "2,4-diF" and "2,6-diF" number the ring from the sulfonamide side (as in vemurafenib's IUPAC name). Numbering from the carbonyl swapped every pair of compounds; the answer key exposed this immediately.
6. **Color carries information that recognizers discard.** In paper 4 Scheme 1, a blue letter C is a placeholder for a linker, not a carbon atom. Detecting blue pixels identified placeholders that MolGlyph had read as ordinary wildcards.
7. **Resolution is not the main lever for recognition errors.** A controlled experiment (same crop boxes and model, 200 vs 300 DPI) gained 1 molecule and lost none, but did not fix miscounted chains or misplaced groups. 300 DPI is the default because it was strictly better on exact match.

## Quick start

```bash
conda create -n biobench -c conda-forge --override-channels python=3.10 -y
conda activate biobench
pip install -r requirements.txt
python run_all.py
```

`run_all.py` rebuilds every result from the committed inputs on a laptop CPU, in about a minute. Rebuilding is deterministic: the output files are byte for byte identical to the committed ones. Results go to `results/A300/` by default.

**Data:** place the PDFs in `data/pdfs/` (named `<paperID>.pdf`) and the answer key at `data/ground_truth.csv`. Both are excluded from Git.

**GPU stages** run on Google Colab, and their outputs are committed, so `run_all.py` does not need a GPU: MolGlyph recognition (`src/recognize.py`) and the OCR tools (`src/read_table_text.py`, `src/read_columns.py`, `src/grid_ocr.py`, `src/page_ocr.py`). See [colab_setup.md](colab_setup.md).

## Pipeline stages

```
PDF pages
  └─ Detection (MolDetv2) ─> crops.csv
       └─ Recognition (MolGlyph, 300 DPI) ─> raw_predictions.csv
            └─ Validation (RDKit) ─> validated_predictions.csv
                 ├─ valid: complete molecules ──────────────────────────────────────┐
                 ├─ scaffold + drawn R-group table ─> assembly ─────────────────────┤
                 ├─ scaffold + text R-group table ─> assembly ──────────────────────┤
                 ├─ scaffold + image table (OCR + table logic) ─> assembly ─────────┤
                 ├─ scaffold + grid table (lines, merged cells, OCR) ─> assembly ───┤
                 ├─ warheads + linkers + ligand pieces ─> modular assembly ─────────┤
                 ├─ scheme placeholders + definition box ─> assembly ───────────────┤
                 ├─ scheme text labels (page OCR) ─> assembly ──────────────────────┤
                 └─ annotated wildcards ─> abbreviation expansion ──────────────────┤
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
| Image table columns | `src/read_columns.py` | Targeted OCR of text columns in image tables, one cell at a time (Colab) |
| Image table assembly | `src/image_tables.py` | Combines column OCR with table logic (consecutive entries, scaffold rules) |
| Table grid | `src/table_grid.py` | Detects a table's ruled lines, rows, columns, and merged cells from the image |
| Grid OCR | `src/grid_ocr.py` | OCR of grid cells, reading each merged cell once (Colab) |
| Grid table assembly | `src/grid_tables.py` | Vocabulary snapping, positional patterns, R2 assigned by merged cell, scaffold reused with bond cutting |
| Modular assembly | `src/pieces.py` | Joins separately drawn warheads, linkers, and ligand pieces |
| Scheme placeholders | `src/scheme_placeholders.py` | Resolves placeholder labels and blue placeholders with the scheme's definition box |
| Page OCR | `src/page_ocr.py` | Full-page OCR with positions, for text inside images (Colab) |
| Scheme labels | `src/scheme_labels.py` | Assembles compounds defined by labels like "19 R = Br; R1 = H" |
| Abbreviations | `src/expand_abbreviations.py` | Expands labels MolGlyph left as wildcards (e.g. SO2C3H7), using its own annotations |
| Scoring | `src/score.py`, `src/evaluate.py` | Final predictions, per paper and per molecule results |
| Review | `src/review_crops.py` | Visual review sheets of crops by paper and status |
| All stages | `run_all.py` | Rebuilds every result in order |

Paper specific settings (table pages, scaffold and piece crops, slot mappings, linker direction, name patterns, scheme boxes) live in `config.yaml`, not in code.

## How results were verified

Every assembly stage was checked against the answer key in at least one of these ways:

1. **By name:** where compound numbers are readable, each built molecule is compared with the reference of the same name.
2. **By structure:** where names are unreadable, each built molecule is matched against the whole paper's answer key, which reveals the name it corresponds to.
3. **Self tests:** known compounds rebuilt by a new route must match exactly (vemurafenib from three different scaffolds; 11a, 11b, and 16 rebuilt by Table 11's route).
4. **Cross route checks:** compound 44 in paper 4 was built from two different tables, and both routes produced the identical molecule.

## Design decisions and known limitations

1. **Slot mappings are verified per table.** MolGlyph's wildcard numbers do not always match a table's R labels. Every mapping in `config.yaml` was confirmed by a self test and by name or structure.
2. **Directions and orientations were tested, not assumed.** Paper 19 Table 3's linker direction was settled by running both directions; paper 4's two-ended linkers are oriented by a chemical rule (the end on N joins the ligand's carbonyl), confirmed on all 7.
3. **Pieces and scaffolds may be reused across tables** when a table's own drawing failed recognition (paper 4 Table 1, paper 19 Table 11), always with a self test.
4. **Inferred values are flagged, never hidden:** an entry number filled from a consecutive run (Table 7's compound 13), an R2 taken from a neighboring merged cell (Table 11's compound 38), and a blue placeholder resolved as linker C (Scheme 1's compounds 27 and 30). All were confirmed by the answer key.
5. **The pipeline refuses rather than guesses.** A compound whose scaffold has an unresolved label is skipped with the reason recorded (Scheme 1's compounds 36 and 37), rather than built wrong.
6. **Generic compounds are kept as wildcards** (e.g. "Heteroaryl/Aryl"), because the answer key stores them that way. Only full size structures with one generic attachment point count.
7. **Fuzzy label matching is limited and audited.** Misread labels are matched to known ones only at 75% similarity or more, and only for labels of 4 or more characters. Every match is recorded.
8. **Open decision: compound 10e** (paper 19). The table writes `CH=CH` without E/Z geometry; the answer key uses E. Left as a known miss pending a documented default.

### Remaining misses (25)

| Paper | Missed | Cause |
|---|---|---|
| 4 | 12 | Blue C read as an ordinary carbon (15, 16, 28, 32); repeated units in brackets (36, 37); a salt (18); named literature PROTACs (JB300, SK2188, dAURK-4, PROTAC-D); 50 |
| 19 | 5 | E/Z geometry (10e), propyl read as ethyl (B9), and a few others |
| 2 | 4 | Molecules outside the main table |
| 16 | 3 | Recognition errors (e.g. an extra methyl) |
| 6 | 1 | A misplaced methyl (NR-7h) |

Most remaining misses are recognition errors, which no assembly step can repair.

## Repository layout

```
config.yaml          data paths, BioMiner paths, table and scheme settings
run_all.py           rebuild all results
src/                 pipeline stages
results/A/           200 DPI run (Phase 1 baseline crops + early Phase 2 stages)
results/A300/        300 DPI run (current default, all stages)
data/                PDFs and answer key (not in Git)
colab_setup.md       GPU stages on Colab
CONTRIBUTIONS.md     who did what
```

## Credits

See [CONTRIBUTIONS.md](CONTRIBUTIONS.md). Models: MolDetv2 (UniParser) and MolGlyph (BioMiner, Yan et al., 2025); chemistry with RDKit; OCR with EasyOCR.
