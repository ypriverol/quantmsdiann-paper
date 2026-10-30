"""Supplementary Fig. S1: (a) per-step runtime on the EBI cluster next to (b) the MDC cluster run.

Inputs: analysis/figures/performance/runtime_per_step.svg (stage performance_trace) and
mdc_cluster_runtime.svg (stage mdc_cluster_runtime). Output: FigS1_runtime_per_step.png/.svg
in the same folder (the SVG wraps the PNG so Word shows both panels).

Usage: python -m scripts.compose_figs1
"""
from __future__ import annotations

import base64
import io
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

PERF = Path(__file__).resolve().parents[1] / "analysis/figures/performance"
HEIGHT = 2400  # px per panel


def render(svg: Path) -> Image.Image:
    png = subprocess.run(["rsvg-convert", "-h", str(HEIGHT), "-b", "white", str(svg)],
                         check=True, capture_output=True).stdout
    return Image.open(io.BytesIO(png)).convert("RGB")


def main() -> int:
    a, b = render(PERF / "runtime_per_step.svg"), render(PERF / "mdc_cluster_runtime.svg")
    pad, top = 120, 160
    out = Image.new("RGB", (a.width + b.width + 3 * pad, HEIGHT + top), "white")
    out.paste(a, (pad, top))
    out.paste(b, (2 * pad + a.width, top))
    draw = ImageDraw.Draw(out)
    try:
        font = ImageFont.truetype("Arial Bold.ttf", 110)
    except OSError:
        font = ImageFont.load_default()
    draw.text((pad // 2, 20), "(a)", fill="black", font=font)
    draw.text((pad + a.width + pad // 2, 20), "(b)", fill="black", font=font)
    png = PERF / "FigS1_runtime_per_step.png"
    out.save(png, dpi=(600, 600))
    w, h = out.size
    data = base64.b64encode(png.read_bytes()).decode()
    (PERF / "FigS1_runtime_per_step.svg").write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{w / 6:.1f}" height="{h / 6:.1f}" viewBox="0 0 {w} {h}">'
        f'<image width="{w}" height="{h}" xlink:href="data:image/png;base64,{data}"/></svg>\n')
    print(png, w, h)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
