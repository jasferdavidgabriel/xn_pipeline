"""Config B, stage 1: extract figure and table images from the PDFs with Docling.

Settings follow the course benchmark's Methods (Gulluoglu et al., 2026): Heron layout model,
page images rendered at images_scale 3.0, the picture detection threshold lowered to 0.4, and
images under 5 KB dropped as noise. Their extra pass cropping individual table cells is not
repeated here, because stage 2 runs a molecule detector inside every table image anyway.
"""
import argparse
import io
from pathlib import Path

import pandas as pd
import yaml

POINTS_TO_PX = 200 / 72   # PDF points to 200 DPI pixels, the units of Config A's crop boxes


def lower_picture_threshold(value):
    """Lower Docling's confidence filter for the Picture label, if this version exposes it."""
    try:
        from docling.utils.layout_postprocessor import LayoutPostprocessor
        from docling_core.types.doc import DocItemLabel
        old = LayoutPostprocessor.CONFIDENCE_THRESHOLDS[DocItemLabel.PICTURE]
        LayoutPostprocessor.CONFIDENCE_THRESHOLDS[DocItemLabel.PICTURE] = value
        print(f"Picture confidence threshold: {old} -> {value}")
    except (ImportError, AttributeError, KeyError) as exc:
        print(f"Picture threshold not adjustable in this version ({exc.__class__.__name__}); using defaults")


def build_converter(images_scale):
    from docling.datamodel.base_models import InputFormat
    from docling.datamodel.layout_model_specs import DOCLING_LAYOUT_HERON
    from docling.datamodel.pipeline_options import PdfPipelineOptions
    from docling.document_converter import DocumentConverter, PdfFormatOption
    opts = PdfPipelineOptions()
    opts.images_scale = images_scale
    opts.generate_picture_images = True
    if hasattr(opts, "generate_table_images"):
        opts.generate_table_images = True
    opts.layout_options.model_spec = DOCLING_LAYOUT_HERON
    return DocumentConverter(format_options={InputFormat.PDF: PdfFormatOption(pipeline_options=opts)})


def extract_figures(pdf_dir, results_dir, paper_ids, images_scale=3.0, threshold=0.4, min_bytes=5 * 1024):
    results_dir = Path(results_dir)
    image_dir = results_dir / "figures"
    image_dir.mkdir(parents=True, exist_ok=True)
    lower_picture_threshold(threshold)
    converter = build_converter(images_scale)

    rows = []
    for pid in paper_ids:
        doc = converter.convert(str(Path(pdf_dir) / f"{pid}.pdf")).document
        counts = {"picture": 0, "table": 0, "dropped": 0}
        for kind, items in (("picture", doc.pictures), ("table", doc.tables)):
            for item in items:
                img = item.get_image(doc)
                if img is None or not item.prov:
                    continue
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                if buf.tell()