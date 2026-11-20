# Running the GPU stages on Google Colab

Two stages use models that are slow on a laptop CPU, so they run on a free Colab T4 GPU:

1. **MolGlyph recognition** of table cells that detection missed (`src/recognize.py`)
2. **Targeted OCR** of table text inside images (`src/read_table_text.py`)

Their outputs are committed under `results/A/`, so `run_all.py` works without a GPU. Follow this guide only to regenerate them.

## One time setup

1. **Hugging Face access.** The MolGlyph weights are gated. Create a free account at huggingface.co, open `jiaxianustc/MolGlyph`, and accept the conditions. Then create a **Read** access token (profile picture, Settings, Access Tokens).
2. **Colab runtime.** Runtime, Change runtime type, T4 GPU.
3. **Colab secret.** In the key panel on the left, add a secret named `HF_TOKEN` with your token, and switch on Notebook access.

## Part A: MolGlyph recognition

MolGlyph needs torch 1.13, which does not support Colab's default Python, so it runs in a separate Python 3.10 environment created with micromamba.

```python
# Cell 1: code
!git clone https://github.com/jasferdavidgabriel/xn_pipeline.git
!git clone https://github.com/jiaxianyan/BioMiner.git /content/BioMiner-main
```

```python
# Cell 2: Python 3.10 environment
!curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
!./bin/micromamba create -y -q -p /content/py310 -c conda-forge python=3.10
PY = "/content/py310/bin/python"
```

```python
# Cell 3: pinned packages (same versions as the Phase 1 baseline)
!{PY} -m pip install -q torch==1.13.1 torchvision==0.14.1 numpy==1.26.4 opencv-python==4.11.0.86 \
    SmilesPE==0.0.3 OpenNMT-py==2.2.0 albumentations==1.1.0 timm==0.4.12
!{PY} -m pip install -q /content/BioMiner-main/BioMiner/MolScribe
!{PY} -m pip install -q rdkit pandas pyyaml pymupdf huggingface_hub
!{PY} -c "import torch; print('GPU available:', torch.cuda.is_available())"
```

```python
# Cell 4: weights, plus a matplotlib fix
import os
from google.colab import userdata
os.environ["HF_TOKEN"] = userdata.get("HF_TOKEN")
os.environ["MPLBACKEND"] = "Agg"   # Colab's default plotting backend is missing in the py310 env
!{PY} -c "from huggingface_hub import hf_hub_download; hf_hub_download('jiaxianustc/MolGlyph', 'molglyph_large.pt', local_dir='/content/BioMiner-main')"
```

```python
# Cell 5: recognize the missed table cells
%cd /content/xn_pipeline
!{PY} -m src.recognize \
    --input results/A/missing_cells_paper2_page5.csv \
    --output results/A/cell_predictions_paper2_page5.csv \
    --id-col cell_id --biominer-dir /content/BioMiner-main --device cuda
```

## Part B: Targeted OCR

EasyOCR works with Colab's default Python, so no special environment is needed.

```python
!pip install -q easyocr rdkit pymupdf pyyaml
```

Upload the paper's PDF to `data/pdfs/` (PDFs are not in Git), then:

```python
%cd /content/xn_pipeline
!python -m src.read_table_text --paper 2 --page 5
```

This writes `table_text_paper2_page5.csv` (compound numbers and X atoms) and adds dash detection to `missing_cells_paper2_page5.csv`.

## Bringing results back

```python
from google.colab import files
files.download("results/A/cell_predictions_paper2_page5.csv")
files.download("results/A/table_text_paper2_page5.csv")
files.download("results/A/missing_cells_paper2_page5.csv")
```

Move the files into `results/A/` on your machine, then run `python run_all.py`.

## Notes from development

1. **Full page OCR drops single characters.** Docling (with EasyOCR) found the paper 2 table's structure (21 rows, 7 columns) and read multi character text, but skipped the digits 1 to 9 and the letters C and N. Recognizing small cropped regions one at a time, with a restricted character list, solved this.
2. **Docling needs cell matching on** when full page OCR is used, or table cells come back empty.
3. **Dashes are found by ink height, not OCR.** A dash cell's dark pixels span about 2% of the cell height, versus 46% or more for any drawn substituent.
