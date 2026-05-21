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


### R-Functions

# If you do not know the bit-depth or max value: Just run FXVI_scaled(p, q). The function will find the highest peak in your data and use it for $K$.

FXVIscaled <- function(p, q, K = NULL) {
  # 1. Handle K dynamically if not provided
  if (is.null(K)) {
    # Find the absolute maximum value across both input vectors
    # suppressWarnings prevents a warning if both vectors happen to be entirely NA
    K <- suppressWarnings(max(c(p, q), na.rm = TRUE))
    
    # Fallback safety: if the data is entirely empty/NA, set K to 1
    if (is.infinite(K)) K <- 1 
  }
  
  # 2. Calculate raw mathematical components
  leading_fraction <- ((p - q) * abs(p - q)) / (p + q)
  log_arg <- (K^2 - p * q) / (p - q)^2
  
  # Suppress expected log() warnings when p and q are exactly equal
  raw_values <- suppressWarnings(leading_fraction * (1 + 0.5 * log(log_arg)))
  
  # 3. Find the true boundaries for scaling
  finite_vals <- raw_values[is.finite(raw_values)]
  
  # Safety check: If the user passes only a single valid value
  if (length(finite_vals) <= 1) {
    return(rep(0, length(raw_values)))
  }
  
  min_val <- min(finite_vals, na.rm = TRUE)
  max_val <- max(finite_vals, na.rm = TRUE)
  
  # Handle the edge case where the math produces absolutely no variance
  if (min_val == max_val) {
    return(rep(0, length(raw_values)))
  }
  
  # 4. Min-Max scale to exactly [-1, 1]
  scaled_values <- 2 * ((raw_values - min_val) / (max_val - min_val)) - 1
  
  # 5. Handle mathematical singularities (where p == q)
  # By standard index convention, perfectly identical inputs equal 0
  scaled_values[!is.finite(scaled_values)] <- 0
  
  return(scaled_values)
}


# "I am trying to debug the pure mathematical theory of my equation."

FXVI <- function(p, q, K = 255) {
  leading_fraction <- ((p - q) * abs(p - q)) / (p + q)
  log_arg <- (K^2 - p * q) / (p - q)^2
  
  raw_values <- suppressWarnings(leading_fraction * (1 + 0.5 * log(log_arg)))
  
  # Safely handle the p=q edge cases
  raw_values[!is.finite(raw_values)] <- 0
  
  return(raw_values)
}

# "I just want to plot it next to NDVI or make a standard map."

FXVIminmax <- function(p, q, K = 255) {
  leading_fraction <- ((p - q) * abs(p - q)) / (p + q)
  log_arg <- (K^2 - p * q) / (p - q)^2
  
  raw_values <- suppressWarnings(leading_fraction * (1 + 0.5 * log(log_arg)))
  
  finite_vals <- raw_values[is.finite(raw_values)]
  
  # Safety check for single-value inputs
  if (length(finite_vals) <= 1) return(rep(0, length(raw_values)))
  
  min_val <- min(finite_vals, na.rm = TRUE)
  max_val <- max(finite_vals, na.rm = TRUE)
  
  if (min_val == max_val) return(rep(0, length(raw_values)))
  
  # Min-Max scale to [-1, 1]
  scaled <- 2 * ((raw_values - min_val) / (max_val - min_val)) - 1
  
  # Safely handle the p=q edge cases
  scaled[!is.finite(scaled)] <- 0
  
  return(scaled)
}

# "I am doing machine learning, rank-based stats, or need maximum visual contrast."

FXVIspaced <- function(p, q, K = 255) {
  leading_fraction <- ((p - q) * abs(p - q)) / (p + q)
  log_arg <- (K^2 - p * q) / (p - q)^2
  
  raw_values <- suppressWarnings(leading_fraction * (1 + 0.5 * log(log_arg)))
  
  n_obs <- length(raw_values)
  
  # Safety check for single-value inputs
  if (n_obs <= 1) return(rep(0, n_obs))
  
  # Rank the values and stretch evenly between -1 and 1
  spaced <- (rank(raw_values, ties.method = "average", na.last = "keep") - 1) / (n_obs - 1) * 2 - 1
  
  # Safely handle the p=q edge cases
  spaced[!is.finite(spaced)] <- 0
  
  return(spaced)
}

# "I need to run a parametric statistical test (like an ANOVA, t-test, or Linear Regression)."

FXVInormal <- function(p, q, K = 255) {
  leading_fraction <- ((p - q) * abs(p - q)) / (p + q)
  log_arg <- (K^2 - p * q) / (p - q)^2
  
  raw_values <- suppressWarnings(leading_fraction * (1 + 0.5 * log(log_arg)))
  
  n_obs <- length(raw_values)
  
  # Safety check for single-value inputs
  if (n_obs <= 1) return(rep(0, n_obs))
  
  # Calculate percentiles and map to a standard normal distribution
  percentile <- (rank(raw_values, ties.method = "average", na.last = "keep") - 0.5) / n_obs
  normed <- suppressWarnings(qnorm(percentile))
  
  # Safely handle the p=q edge cases
  normed[!is.finite(normed)] <- 0
  
  return(normed)
}


# VIBE Index: Symmetrically scaled to [-1, 1] to strictly preserve directional sign.

VIBE <- function(p, q, K = NULL) {
  # 1. Auto-detect K if not provided
  if (is.null(K)) {
    K <- suppressWarnings(max(c(p, q), na.rm = TRUE))
    if (is.infinite(K)) K <- 1 
  }
  
  # 2. Core Math
  denominator <- p + q
  
  # suppressWarnings safely catches the 0/0 edge case when p=0 and q=0
  fraction <- suppressWarnings((p - q) / denominator)
  
  # sign() natively returns 1, -1, or 0
  raw_values <- denominator * sign(fraction)
  
  # 3. Symmetrical Scaling to [-1, 1]
  # The absolute theoretical maximum of (p + q) is 2 * K
  scaled_values <- raw_values / (2 * K)
  
  # 4. Enforce strict mathematical conventions
  # V(p,p) = 0 is naturally handled by the math because sign(0) = 0
  # V(0,0) is explicitly marked as undefined (NA) per the formula rules
  scaled_values[p == 0 & q == 0] <- NA
  
  return(scaled_values)
}

# I want more contrast, or do edge-detection.
# The standard, URFI Index: Applies rank-preserving analytical compression to the normalized ratio.

URFI <- function(p, q) {
  # 1. Calculate the base normalized ratio (x)
  denominator <- p + q
  
  # suppressWarnings catches the 0/0 edge case when p=0 and q=0
  x <- suppressWarnings((p - q) / denominator)
  
  # Handle the undefined 0/0 case (standard practice is NA)
  x[p == 0 & q == 0] <- NA
  
  # 2. Core Math: Nonlinear Compression. Kaiblinger & Wabnig (2019): Modified NDVI formula for uniformly distributed values. Unpublished Mansucript.
  raw_values <- (2 * x) / (1 + abs(x))
  
  # 3. Enforce the strict 2-decimal rounding rule
  rounded_values <- round(raw_values, 2)
  
  return(rounded_values)
}

Copyright (c) [2026] [Thomas Wabnig]
Wabnig, T. (2026). Distribution-Aware Ratio Indices. https://github.com/wabnigt-byte/Distribution-Aware-Ratio-Indices
