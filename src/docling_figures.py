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
    from docling.datamodel.pipeline_options import PdfPipelineOptions
    from docling.document_converter import DocumentConverter, PdfFormatOption
    opts = PdfPipelineOptions()
    opts.images_scale = images_scale
    opts.generate_picture_images = True
    opts.generate_page_images = True      # table images are cropped from page images in this version
    # Heron is the default layout model in this version; print it so the run records what was used
    print("Layout model:", opts.layout_options.model_spec.name)
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
                if buf.tell() < min_bytes:                       # tiny images are noise
                    counts["dropped"] += 1
                    continue
                prov = item.prov[0]
                page_height = doc.pages[prov.page_no].size.height
                box = prov.bbox.to_top_left_origin(page_height=page_height)
                figure_id = f"{pid}_p{prov.page_no}_{kind[0]}{counts[kind]:03d}"
                path = image_dir / f"{figure_id}.png"
                path.write_bytes(buf.getvalue())
                counts[kind] += 1
                rows.append({"paper_id": str(pid), "page": prov.page_no, "figure_id": figure_id, "kind": kind,
                             "x1": round(box.l * POINTS_TO_PX), "y1": round(box.t * POINTS_TO_PX),
                             "x2": round(box.r * POINTS_TO_PX), "y2": round(box.b * POINTS_TO_PX),
                             "images_scale": images_scale, "image_path": str(path)})
        print(f"paper {pid}: {counts['picture']} pictures, {counts['table']} tables, {counts['dropped']} dropped as noise")

    df = pd.DataFrame(rows)
    df.to_csv(results_dir / "figures.csv", index=False)
    print(f"{len(df)} figure and table images -> {results_dir / 'figures.csv'}")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Config B stage 1: Docling figure extraction.")
    parser.add_argument("--results", default="results/B")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--papers", default="1,2,4,6,16,19")
    args = parser.parse_args()
    cfg = yaml.safe_load(open(args.config))
    extract_figures(cfg["data"]["pdf_dir"], args.results, args.papers.split(","))