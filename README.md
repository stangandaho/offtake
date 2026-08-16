# offtake

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/stangandaho/offtake/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stangandaho/offtake/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/stangandaho/offtake/graph/badge.svg)](https://app.codecov.io/gh/stangandaho/offtake)
<!-- badges: end -->

**offtake** is a tidyverse-friendly R toolkit for assessing the sustainability
of wildlife exploitation (bushmeat hunting, recreational hunting and
fisheries). It implements the main index-based and model-based in the systematic 
review of **Adounkè et al. (2026)** and the methodological
review of **Weinbaum et al. (2013)**.

Every function takes a data frame and returns a classed tibble with an explicit `sustainable`
verdict, so results drop straight into a `dplyr`/`ggplot2` workflow.

## Installation

```r
# install.packages("pak")
pak::pkg_install("stangandaho/offtake")
```

## Methods implemented

All functions use the `ot_` prefix (offtake) for a distinct,
recognisable API.

| Function      | Family | Method | Source |
|---------------|--------|--------|--------|
| `ot_asc()`    | index  | Age structure comparison (hunted vs reference; lower juvenile share = unsustainable) | Adounkè et al. 2026; Weinbaum et al. 2013 |
| `ot_hyco()`   | index  | Hunting yield / CPUE comparison (lower yield at more-hunted site = unsustainable) | Adounkè et al. 2026; Weinbaum et al. 2013 |
| `ot_pdc()`    | index  | Population density comparison (lower density at hunted site = unsustainable) | Adounkè et al. 2026; Weinbaum et al. 2013 |
| `ot_pro()`    | model  | Production model, `P = 0.6·K·(λmax − 1)·F` | Robinson & Redford 1991 |
| `ot_pbr()`    | model  | Potential biological removal, `PBR = Nmin·½·Rmax·FR` | Wade 1998 |
| `ot_biode()`  | model  | Spatial source–sink biodemographic model | Levi et al. 2009 |
| `ot_msy()`    | model  | Logistic maximum sustainable yield, `MSY = rK/4` *(optional add-on)* | Schaefer 1954 |
| `ot_samse()`  | model  | Stochastic sustainable-mortality limit *(optional add-on)* | Manlik et al. 2022 |

Helper `ot_lambda_max()` solves Cole's (1954) equation for the maximum finite
rate of increase used by `ot_pro()`.

## Quick start

```r
library(offtake)

# Index: is density lower at the hunted site?
ot_pdc(bushmeat_sites, density = density, group = site_type,
       reference = "reference")

# Model: does offtake exceed the Robinson & Redford production?
ot_pro(duiker_demography, k = density_k, harvest = annual_take,
       b = b, a = a, w = w, longevity = lifespan) # λmax from Cole's equation

# Spatial model: map hunting-induced depletion
res <- ot_biode(manu_settlements, x = x_km, y = y_km, humans = hunters,
                k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6)
res
ot_biode_surface(res) # per-cell density surface, ready for ggplot2::geom_raster()
```

## References

- Adounkè, G. R. M. *et al.* (2026) Systematic review of sustainability
  assessment approaches for wildlife exploitation. *Biological Conservation*
  313, 111606. <https://doi.org/10.1016/j.biocon.2025.111606>
- Weinbaum, K. Z. *et al.* (2013) Searching for sustainability: are assessments
  of wildlife harvests behind the times? *Ecology Letters* 16, 99–111.
  <https://doi.org/10.1111/ele.12008>
- Robinson, J. G. & Redford, K. H. (1991) Sustainable harvest of Neotropical
  forest mammals.
- Wade, P. R. (1998) *Marine Mammal Science* 14, 1–37.
- Levi, T. *et al.* (2009) *Journal of Applied Ecology* 46, 804–814.
- Manlik, O. *et al.* (2022) *Conservation Biology* 36, e13897.
- Cole, L. C. (1954) *The Quarterly Review of Biology* 29, 103–137.
