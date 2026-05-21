# Distribution-Aware Ratio Indices 

A small collection of experimental ratio indices for R, where the core math is fixed but the output is reshaped to a chosen probability mass
function (PMF). Not arbitrary sample distribution! The idea: the same band ratio (e.g. in remote sensing Vegetation Indices) can be delivered raw, scaled,
rank-uniform, or Gaussianized depending on what you intend to do with it (mapping, ML, or parametric statistics).

The usage is NOT limited to remote sensing, but any other field were ratios are used.

All functions take two numeric vectors `p` and `q` (e.g. two reflectance bands)
and are fully vectorized. Edge cases (`p == q`, `0/0`, all-`NA`, single value)
are handled internally.

## Indices

### FXVI family (Functional Vegetation Index) is meant to behave like a bijective function.
## Imagined use: Test if it mitigates NDVI saturation.
A modified index combining an NDVI-style term with a Canberra-distance log term:

```
((p - q) * |p - q|) / (p + q)  *  (1 + 0.5 * log( (K^2 - p*q) / (p - q)^2 ))
```

`K` is the max possible value (bit-depth dependent): 255 (8-bit), 1023 (10-bit),
65535 (16-bit). If unknown, the `scaled` variant auto-detects it from the data.

The five variants share this math and differ only in post-processing:

| Function       | Output PMF / range           | Use when…                                               |
|----------------|------------------------------|---------------------------------------------------------|
| `FXVI`         | Raw math (unbounded)         | Debugging the pure theory of the equation               |
| `FXVIminmax`   | Min-max to [-1, 1]           | Plotting next to NDVI / standard maps                   |
| `FXVIspaced`   | Rank → uniform [-1, 1]       | ML, rank-based stats, maximum visual contrast           |
| `FXVInormal`   | Rank → standard normal       | Parametric tests (ANOVA, t-test, linear regression)     |
| `FXVIscaled`   | Min-max [-1, 1], dynamic `K` | You don't know the bit-depth (auto-detects `K`)         | <- Standard

### VIBE (Vegetation Index for Biomass Estimation)
## Imagined use: Biomass estimation
Symmetrically scaled to [-1, 1], strictly preserving directional sign.
Output is bimodal. `V(p, p) = 0`; `V(0, 0)` is returned as `NA` (undefined).

### URFI (Uniform Relative Frequency Index)
## Imagined use: Edge Detection, Higher contrast (in UAV images)
A rank-preserving analytical compression of the normalized ratio
`x = (p - q)/(p + q)` via `(2x)/(1 + |x|)`, rounded to 2 decimals. Produces a
near-uniform spread. Based on Kaiblinger & Wabnig (2019): Modified NDVI formula
for uniformly distributed values (unpublished manuscript).


## Notes
- `K` must match your data's bit-depth for `FXVI`, `FXVIminmax`, `FXVIspaced`,
  and `FXVInormal` (default `255`). `FXVIscaled` and `VIBE` infer it if `K = NULL`.
- These are exploratory / experimental indices, not validated remote-sensing
  products — use accordingly.


Copyright (c) [2026] [Thomas Wabnig]
Wabnig, T. (2026). Distribution-Aware Ratio Indices. https://github.com/wabnigt-byte/Distribution-Aware-Ratio-Indices
