# Running the GPU stages on Google Colab

These stages use models that are slow on a laptop CPU, so they run on a free Colab T4 GPU. Their outputs are committed under `results/A300/`, so `run_all.py` works without a GPU. Follow this guide only to regenerate them.

| Stage | Tool | Output |
|---|---|---|
| Render crops at 300 DPI | `src/rebuild_crops.py` | `results/A300/crops/` |
| MolGlyph recognition | `src/recognize.py` | `raw_predictions.csv`, `cell_predictions_*.csv` |
| Table text OCR (paper 2) | `src/read_table_text.py` | `table_text_paper2_page5.csv` |
| Image table columns (paper 19 Table 7) | `src/read_columns.py` | `table_cells_paper19_page7.csv` |
| Grid OCR (paper 19 Table 11) | `src/grid_ocr.py` | `grid_cells_paper19_page9.csv` |
| Page OCR (paper 4 Scheme 1) | `src/page_ocr.py` | `page_ocr_paper4_page3.csv` |

## One time setup

1. **Hugging Face access.** The MolGlyph weights are gated. Create a free account at huggingface.co, open `jiaxianustc/MolGlyph`, and accept the conditions. Then create a **Read** access token (profile picture, Settings, Access Tokens).
2. **Colab runtime.** Runtime, Change runtime type, T4 GPU.
3. **Colab secret.** In the key panel on the left, add a secret named `HF_TOKEN` with your token, and switch on Notebook access.

## Every session

```python
!git clone https://github.com/jasferdavidgabriel/xn_pipeline.git
%cd /content/xn_pipeline
```

The PDFs are not in Git, so upload the ones you need into `data/pdfs/`:

```python
import os, shutil
from google.colab import files
os.makedirs("data/pdfs", exist_ok=True)
for name in files.upload():
    shutil.move(name, f"data/pdfs/{name}")
```

If a session restarts, rerun the variable cell at the top of Part A before any MolGlyph command.

## Part A: MolGlyph recognition

MolGlyph needs torch 1.13, which does not support Colab's default Python, so it runs in a separate Python 3.10 environment created with micromamba.

```python
# A1: BioMiner code and a Python 3.10 environment
!git clone https://github.com/jiaxianyan/BioMiner.git /content/BioMiner-main
!curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
!./bin/micromamba create -y -q -p /content/py310 -c conda-forge python=3.10
```

```python
# A2: pinned packages (same versions as the Phase 1 baseline)
PY = "/content/py310/bin/python"
!{PY} -m pip install -q torch==1.13.1 torchvision==0.14.1 numpy==1.26.4 opencv-python==4.11.0.86 \
    SmilesPE==0.0.3 OpenNMT-py==2.2.0 albumentations==1.1.0 timm==0.4.12
!{PY} -m pip install -q /content/BioMiner-main/BioMiner/MolScribe
!{PY} -m pip install -q rdkit pandas pyyaml pymupdf huggingface_hub
!{PY} -c "import torch; print('GPU available:', torch.cuda.is_available())"
```

```python
# A3: variables (rerun after any restart) and weights
import os
from google.colab import userdata
PY = "/content/py310/bin/python"
os.environ["HF_TOKEN"] = userdata.get("HF_TOKEN")
os.environ["MPLBACKEND"] = "Agg"   # Colab's default plotting backend is missing in the py310 env
!{PY} -c "from huggingface_hub import hf_hub_download; hf_hub_download('jiaxianustc/MolGlyph', 'molglyph_large.pt', local_dir='/content/BioMiner-main')"
```

```python
# A4: render all crops at 300 DPI, then recognize them
!pip install -q pymupdf pyyaml
!python -m src.rebuild_crops --results results/A300 --dpi 300
!{PY} -m src.recognize --input results/A300/crops.csv --output results/A300/raw_predictions.csv \
    --biominer-dir /content/BioMiner-main --device cuda
```

```python
# A5: recognize table cells that detection missed (paper 2)
!{PY} -m src.recognize --input results/A300/missing_cells_paper2_page5.csv \
    --output results/A300/cell_predictions_paper2_page5.csv \
    --id-col cell_id --biominer-dir /content/BioMiner-main --device cuda
```

## Part B: OCR tools

EasyOCR works with Colab's default Python, so no special environment is needed.

```python
!pip install -q easyocr rdkit pymupdf pyyaml
```

```python
# B1: paper 2 table text (compound numbers, X atoms, dash detection)
!python -m src.read_table_text --results results/A300 --paper 2 --page 5

# B2: paper 19 Table 7 text columns (entry, scaffold, ring atom)
!python -m src.read_columns --paper 19 --page 7 \
  --columns '{"entry": [395, 460, "0123456789"], "scaffold": [450, 520, "AB"], "A": [530, 640, "NCH-"]}'

# B3: paper 19 Table 11 grid cells (needs grid_paper19_page9.json from src/table_grid.py)
!python -m src.grid_ocr --paper 19 --page 9 \
  --columns '{"R1": "nCHbenzyl0123456789()-F", "A": "0123456789,-diF", "X": "CHN"}'

# B4: paper 4 Scheme 1 labels
!python -m src.page_ocr --paper 4 --page 3
```

## Bringing results back

Download each output, move it into `results/A300/` on your machine, and run `python run_all.py`:

```python
from google.colab import files
files.download("results/A300/raw_predictions.csv")   # and any other outputs you regenerated
```

## Notes from development

1. **Full page OCR drops single characters** (digits 1 to 9, the letters C and N). Recognizing small cropped regions one at a time, with a restricted character list, solves this.
2. **Thin table text renders as light gray at 200 DPI**, which OCR misreads. The column and grid tools render at 400 DPI and stretch contrast.
3. **OCR on tiny cells is still unreliable alone.** The assembly stages combine it with each table's own logic and with small vocabularies, and flag every inferred value.
4. **Dashes are found by ink height, not OCR.** A dash cell's dark pixels span about 2 to 4% of the cell height, versus 24% or more for any letter.
5. **Docling needs cell matching on** when full page OCR is used, or table cells come back empty. It was explored for paper 2's image table but not used in the final pipeline.
