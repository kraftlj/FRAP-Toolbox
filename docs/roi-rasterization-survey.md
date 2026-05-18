# ROI Rasterization Survey

This note records the follow-up search for a circular ROI rasterizer that can
replace the current `skimage.draw.disk` implementation while reproducing the
MATLAB guide exports.

## Parity Target

The diffusion guide fixture `Venus_Cytoplasm_1.lsm` uses a circular bleach ROI
with center `(256, 23)` and radius `9` in MATLAB's 1-based UI coordinates. The
exported Raw FRAP vector implies a 252-pixel bleach mask. Supplying that
inferred mask to the Python loader reproduces the MATLAB Time, Raw FRAP,
Normalized FRAP, Distance, and Post-bleach Profile vectors to export precision.

That means the remaining strict live-loader parity gap is specifically the
circle-to-mask rasterization step, not BioIO/LSM loading, timestamp extraction,
profile construction, normalization, or photodecay correction.

## Tools Checked

- OpenCV drawing APIs:
  - `cv2.circle(..., thickness=cv2.FILLED)`
  - `cv2.circle(..., shift=...)` with fixed-point subpixel centers/radii
  - `cv2.fillPoly(...)` with MATLAB's 41-point circle polygon
  - `cv2.ellipse2Poly(...)` followed by `cv2.fillPoly(...)`
- scikit-image:
  - `skimage.draw.disk`
  - `skimage.draw.polygon`
  - `skimage.draw.polygon2mask`
- Pillow:
  - `ImageDraw.ellipse`
  - `ImageDraw.polygon`
- matplotlib diagnostic:
  - `Path.contains_points(...)` with candidate pixel-center offsets

Official API references:

- OpenCV drawing functions: https://docs.opencv.org/3.4/d6/d6e/group__imgproc__draw.html
- scikit-image draw API: https://scikit-image.org/docs/stable/api/skimage.draw.html
- Pillow ImageDraw API: https://pillow.readthedocs.io/en/stable/reference/ImageDraw.html
- MATLAB `poly2mask`: https://www.mathworks.com/help/images/ref/poly2mask.html

## Best Local Matches

All counts are compared against the inferred 252-pixel MATLAB guide mask.

| Tool / method | Mask count | Missing target pixels | Extra pixels | Total differing pixels |
| --- | ---: | ---: | ---: | ---: |
| Inferred MATLAB guide mask | 252 | 0 | 0 | 0 |
| Tuned `matplotlib.path.Path.contains_points` diagnostic | 254 | 0 | 2 | 2 |
| Tuned `matplotlib.path.Path.contains_points` diagnostic | 250 | 2 | 0 | 2 |
| OpenCV subpixel `cv2.circle` scan | 257 | 6 | 11 | 17 |
| OpenCV `cv2.circle`, center `(255, 22)`, radius `9` | 253 | 11 | 12 | 23 |
| Current `skimage.draw.disk` | 249 | 13 | 10 | 23 |
| `skimage.draw.polygon` with MATLAB circle polygon, offset `-1` | 253 | 11 | 12 | 23 |
| OpenCV `ellipse2Poly` plus `fillPoly` | 241 | 19 | 8 | 27 |
| Pillow `ImageDraw.ellipse` best bounding box | 277 | 2 | 27 | 29 |
| OpenCV `fillPoly` with fixed-point MATLAB circle polygon | 277 | 4 | 29 | 33 |

## Conclusion

OpenCV is not a drop-in fix for this parity gap. Its filled circle and polygon
scan-conversion behavior is useful for modern image processing, but the tested
variants do not reproduce the MATLAB `poly2mask` mask implied by the guide
exports.

The next production-quality path is to implement and test a MATLAB-compatible
polygon scan conversion mode for circular diffusion and FRAP-FRET ROIs. The
test should target the 252-pixel inferred guide mask first, then expand to
small synthetic polygons/circles where MATLAB or Octave reference masks can be
generated. Until that is complete, the inferred mask test should remain as the
diagnostic proof that the rest of the loader is already at parity.
