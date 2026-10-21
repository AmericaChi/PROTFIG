"""
make_gallery.py — gallery of the example figures for the documentation.

  python3 make_gallery.py [figures_dir] [docs_dir]      (from examples/; or with pymol -cq ... --)

Writes <docs>/gallery.jpg (all figures in a grid, labeled with style/name) and one
thumbnail per figure in <docs>/presets/<style>__<name>.jpg for the README.
Requires Pillow (pip install pillow).
"""
import glob
import os
import sys

from PIL import Image, ImageDraw, ImageFont

FIG = sys.argv[1] if len(sys.argv) > 1 else "figures"
DOCS = sys.argv[2] if len(sys.argv) > 2 else os.path.join("..", "docs")
CELL, COLS, PAD = 420, 5, 10
os.makedirs(os.path.join(DOCS, "presets"), exist_ok=True)
try:
    font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 18)
except OSError:
    font = ImageFont.load_default()


def flat(path):
    im = Image.open(path).convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    return Image.alpha_composite(bg, im).convert("RGB")


pngs = sorted(p for p in glob.glob(os.path.join(FIG, "*", "*", "*.png"))
              if os.path.basename(p)[:-4] == os.path.basename(os.path.dirname(p)))
tiles = []
for p in pngs:
    preset, name = p.split(os.sep)[-3], p.split(os.sep)[-2]
    im = flat(p)
    thumb = im.copy()
    thumb.thumbnail((900, 900))
    thumb.save(os.path.join(DOCS, "presets", f"{preset}__{name}.jpg"), quality=85)
    im.thumbnail((CELL, CELL - 28))
    cell = Image.new("RGB", (CELL, CELL), "white")
    cell.paste(im, ((CELL - im.width) // 2, 28 + (CELL - 28 - im.height) // 2))
    ImageDraw.Draw(cell).text((6, 4), f"{preset} / {name}", fill=(40, 40, 40), font=font)
    tiles.append(cell)
rows = (len(tiles) + COLS - 1) // COLS
grid = Image.new("RGB", (COLS * CELL + (COLS + 1) * PAD, rows * CELL + (rows + 1) * PAD), (225, 225, 225))
for i, t in enumerate(tiles):
    grid.paste(t, (PAD + (i % COLS) * (CELL + PAD), PAD + (i // COLS) * (CELL + PAD)))
grid.save(os.path.join(DOCS, "gallery.jpg"), quality=85)
print(f"{len(tiles)} figures -> {os.path.join(DOCS, 'gallery.jpg')}")
