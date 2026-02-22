#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


def make_placeholder(path: Path, title: str, subtitle: str) -> None:
    width, height = 1600, 1000
    img = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    # Simple border
    draw.rectangle([40, 40, width - 40, height - 40], outline=(40, 40, 40), width=4)

    # Fonts (fall back to default if unavailable)
    try:
        font_title = ImageFont.truetype("Arial.ttf", 56)
        font_sub = ImageFont.truetype("Arial.ttf", 32)
    except Exception:
        font_title = ImageFont.load_default()
        font_sub = ImageFont.load_default()

    # Centered text
    title_bbox = draw.textbbox((0, 0), title, font=font_title)
    subtitle_bbox = draw.textbbox((0, 0), subtitle, font=font_sub)
    title_w = title_bbox[2] - title_bbox[0]
    title_h = title_bbox[3] - title_bbox[1]
    sub_w = subtitle_bbox[2] - subtitle_bbox[0]
    sub_h = subtitle_bbox[3] - subtitle_bbox[1]

    x_title = (width - title_w) // 2
    y_title = (height - (title_h + 20 + sub_h)) // 2
    x_sub = (width - sub_w) // 2
    y_sub = y_title + title_h + 20

    draw.text((x_title, y_title), title, fill=(20, 20, 20), font=font_title)
    draw.text((x_sub, y_sub), subtitle, fill=(60, 60, 60), font=font_sub)

    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path, format="PNG")


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    out_dir = repo_root / "plots" / "publication" / "png"

    make_placeholder(
        out_dir / "figure1_placeholder.png",
        "Figure 1 (Placeholder)",
        "Study overview and benchmark design",
    )
    make_placeholder(
        out_dir / "figure2_placeholder.png",
        "Figure 2 (Placeholder)",
        "Domain-quality benchmark summary",
    )
    make_placeholder(
        out_dir / "figure3_placeholder.png",
        "Figure 3 (Placeholder)",
        "Stability and sensitivity summaries",
    )
    make_placeholder(
        out_dir / "figure4_placeholder.png",
        "Figure 4 (Placeholder)",
        "Compute feasibility and guidance",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

