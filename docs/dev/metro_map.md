# Metro map

The pipeline overview is generated from `assets/metro_map.mmd` with [nf-metro](https://github.com/pinin4fjords/nf-metro). Update the source whenever a public route or reuse boundary changes, then regenerate both committed images from the project development shell (`nix-shell` or direnv):

```bash
python3 -m venv .venv-nf-metro
.venv-nf-metro/bin/pip install 'nf-metro==1.1.0'

.venv-nf-metro/bin/nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-gwas_metro_map.svg \
  --theme light --mode light --x-spacing 90 \
  --responsive --validate --compact-offsets --no-chrome-css \
  --logo docs/images/nf-core-gwas_logo_light.png

cairosvg docs/images/nf-core-gwas_metro_map.svg \
  -o docs/images/nf-core-gwas_metro_map.png --output-width 2200
```

Open the SVG and PNG after rendering and check their text, line labels, contrast and cropping. The README uses the SVG; the PNG is a static presentation fallback.

The SVG and PNG must retain transparent backgrounds. Check both the raw alpha channel and white/dark composites rather than judging transparency from a single image viewer.
