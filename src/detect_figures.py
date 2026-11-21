"""Config B, stage 2: run MolDetv2 inside each Docling figure and map the boxes back to the page.

Detection settings match the Phase 1 baseline (confidence 0.5, image size 960, 10 px margin),
so Config A and Config B differ only in how pages are split into images.
"""
import argparse
from pathlib import Path

import pandas as pd
from PIL import Image

PX200_PER_POINT = 200 / 72


def detect(results_dir, weights, conf=0.5, imgsz=960, margin=10, device=0):
    from ultralytics import YOLO
    results_dir = Path(results_dir)
    figures = pd.read_csv(results_dir / "figures.csv", dtype={"paper_id": str})
    crop_dir = results_dir / "crops"
    crop_dir.mkdir(parents=True, exist_ok=True)
    model = YOLO(weights)

    rows, per_page = [], {}
    for f in figures.itertuples():
        img = Image.open(f.image_path).convert("RGB")
        result = model.predict(img, conf=conf, imgsz=imgsz, device=device, verbose=False)[0]
        scale = PX200_PER_POINT / f.images_scale           # figure pixels -> 200 DPI page pixels
        for (bx1, by1, bx2, by2), score in zip(result.boxes.xyxy.tolist(), result.boxes.conf.tolist()):
            box = (max(0, int(bx1) - margin), max(0, int(by1) - margin),
                   min(img.width, int(bx2) + margin), min(img.height, int(by2) + margin))
            key = (f.paper_id, f.page)
            per_page[key] = per_page.get(key, 0) + 1
            crop_id = f"{f.paper_id}_p{f.page}_{per_page[key]:03d}"   # same ID format as Config A
            path = crop_dir / f"{crop_id}.png"
            img.crop(box).save(path)
            rows.append({"crop_id": crop_id, "figure_id": f.figure_id, "paper_id": f.paper_id, "page": f.page,
                         "x1": round(f.x1 + box[0] * scale), "y1": round(f.y1 + box[1] * scale),
                         "x2": round(f.x1 + box[2] * scale), "y2": round(f.y1 + box[3] * scale),
                         "confidence": round(score, 3), "crop_path": str(path)})

    df = pd.DataFrame(rows)
    df.to_csv(results_dir / "crops.csv", index=False)
    print(df.groupby("paper_id").size().rename("crops").to_string())
    print(f"{len(df)} crops from {len(figures)} figures -> {results_dir / 'crops.csv'}")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Config B stage 2: MolDetv2 inside Docling figures.")
    parser.add_argument("--results", default="results/B")
    parser.add_argument("--weights", required=True)
    parser.add_argument("--cpu", action="store_true")
    args = parser.parse_args()
    detect(args.results, args.weights, device="cpu" if args.cpu else 0)