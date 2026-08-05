# Metro map

The pipeline overview is generated from `assets/metro_map.mmd` with [nf-metro](https://github.com/pinin4fjords/nf-metro). Update the source whenever a public route or reuse boundary changes, then regenerate both committed images:

```bash
python3 -m venv .venv-nf-metro
.venv-nf-metro/bin/pip install 'nf-metro==1.1.0'

.venv-nf-metro/bin/nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-gwas_metro_map.svg \
  --theme nfcore-light --mode light --x-spacing 60 --y-spacing 40 \
  --responsive --validate --compact-offsets \
  --logo docs/images/nf-core-gwas_logo_light.png

magick -background white docs/images/nf-core-gwas_metro_map.svg \
  -resize 2200x -depth 8 -strip docs/images/nf-core-gwas_metro_map.png
```

Open the SVG and PNG after rendering and check their text, line labels, contrast and cropping. The README uses the SVG; the PNG is a static presentation fallback.
