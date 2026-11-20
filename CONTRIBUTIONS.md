# Contributions

XN Project Group 5: Sanika Mahesh Potphode, Yihang Xu, Nishtha Patel, Jasfer David Gabriel.

> **To confirm with the team:** please check each line below and correct anything that's inaccurate or missing before submission.

## Phase 1: Benchmarking

| Member | Contribution |
|---|---|
| Sanika Mahesh Potphode | End-to-end benchmark of MolDetv2 + MolGlyph on the six course PDFs, run on CPU. Its crops and predictions are the **Config A baseline** this pipeline builds on. Group report: aggregated results and comparative analysis. |
| Yihang Xu | MolGlyph recognition benchmark on BioVista images (72% exact match), which this pipeline uses to show recognition on complete molecules matches its clean-image accuracy. Group report: conclusions and sponsor questions. |
| Nishtha Patel | Literature review for the group report. |
| Jasfer David Gabriel | Evaluation pipeline (canonicalization, best match scoring, sanity check). Group report: pipeline architecture and Phase 2 plan. |

## Phase 2: Pipeline development

| Member | Contribution |
|---|---|
| Jasfer David Gabriel | Backend pipeline: shared CSV interfaces, validation statuses, crop review, drawn and text R-group table assembly, missed cell recovery, targeted OCR, abbreviation expansion, scoring, `run_all.py`, documentation. |
| *(add)* | *(e.g. reviewing crops, labeling bounding boxes, checking table mappings, report sections)* |

## Development log

| Milestone | Overall recall | Precision |
|---|---|---|
| Reproduced Phase 1 baseline (27.3%), matching Sanika's results to four decimal places | 27.3% | 27.2% |
| Paper 2 scaffold + R-group table assembly | 35.5% | 82.5% |
| Recognized missed table cells (Colab) | 36.4% | |
| Targeted OCR and dash detection: all 21 paper 2 table compounds correct | 36.8% | |
| Paper 19 text tables: Tables 6, 9, 1, 3 (45 of 46 correct by name) | 56.8% | 87.3% |
| Abbreviation expansion from MolGlyph annotations | 60.0% | 84.9% |
| `run_all.py`: deterministic rebuild of all results | 60.0% | 84.9% |
