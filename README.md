
<!-- badges: start -->

[![R-CMD-check](https://github.com/djghosh1123/ivgls/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/djghosh1123/ivgls/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# ivgls

**ivgls** implements network-aware instrumental variable (IV) regression
with a graph-fused Lasso penalty for causal variable selection in
high-dimensional, graph-structured settings.

## Estimators

| Function     | Graph penalty | Invalid-IV robust |
|--------------|:-------------:|:-----------------:|
| `iv_lasso()` |      No       |        No         |
| `ivgl()`     |      Yes      |        No         |
| `ivgl_s()`   |      Yes      |        Yes        |

## Installation

`glmgraph` is required but not on CRAN — install it first:

``` r
devtools::install_github("cran/glmgraph")
install.packages("ivgls")
```

## Quick Example

``` r
library(ivgls)

set.seed(1)
A    <- make_graph(p = 20, type = "chain")
L    <- get_laplacian(A)
bobj <- generate_beta(A, s2 = 4, signal = 2)
dat  <- generate_data(n = 120, p = 20, q = 60,
                      s_alpha = 5, alpha_strength = 3,
                      beta_true = bobj$beta_true)

fit <- ivgl_s(dat$Y, dat$X, dat$Z, L)
get_mcc(bobj$active_set, which(abs(fit$beta) > 1e-4), p = 20)
```

## Citation

> Pal, S. & Ghosh, D. (2026). Network-aware IV regression for causal
> node discovery and estimation.
