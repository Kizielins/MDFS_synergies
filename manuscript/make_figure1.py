#!/usr/bin/env python3
"""Assemble Figure 1 from its panels: a) schematic, b) Figure1b.png (all directed pairs), c) Figure1c.png (additivity).
The schematic is taken from Figure1a_schematic.png if present (export the slide from Figure1a_schematic.pptx at high
resolution); otherwise from the panel-a area of the previously embedded figure (low resolution, draft quality)."""
import os, sys
from PIL import Image, ImageDraw, ImageFont
HERE = os.path.dirname(os.path.abspath(__file__)); WIDTH = 3200
def fit(im): return im.resize((WIDTH, round(im.height * WIDTH / im.width)), Image.LANCZOS)
src = os.path.join(HERE, "Figure1a_schematic.png")
a = Image.open(src).convert("RGB") if os.path.exists(src) else None
if a is None:
    prev = sys.argv[1] if len(sys.argv) > 1 else None
    im = Image.open(prev).convert("RGB"); a = im.crop((0, 0, im.width, round(im.height * 0.615)))
panels = [fit(a)] + [fit(Image.open(os.path.join(HERE, f)).convert("RGB")) for f in ("Figure1b.png", "Figure1c.png")]
gap = 90; H = sum(p.height for p in panels) + gap * 2
out = Image.new("RGB", (WIDTH, H), "white"); d = ImageDraw.Draw(out)
try: font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 64)
except Exception: font = ImageFont.load_default()
y = 0
for i, p in enumerate(panels):
    out.paste(p, (0, y))
    if i > 0: d.text((20, y - 70), "bc"[i - 1] + ")", fill="black", font=font)
    y += p.height + gap
out.save(os.path.join(HERE, "Figure1.png"), dpi=(300, 300)); print("Figure1.png", out.size)
