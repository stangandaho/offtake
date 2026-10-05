# Potential biological removal (PBR)

Estimates the maximum number of animals that may be removed from a
population per year while allowing it to reach or stay near its optimal
size, and compares it with the observed human-caused removal. Originally
developed for marine mammal management (Wade 1998) and increasingly
applied to terrestrial harvest (Adounke et al. 2026).

## Usage

``` r
ot_pbr(
  data,
  rmax,
  removal,
  nmin = NULL,
  n = NULL,
  cv = NULL,
  fr = 0.5,
  percentile = 0.2,
  uncertainty = NULL,
  n_sim = 1000,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  A data frame with **one row per population/stock**. Required columns:
  the maximum growth rate (`rmax`) and the observed annual removal
  (`removal`). You must also supply the population size, either as a
  minimum estimate (`nmin`) or as a point estimate (`n`) with an
  optional coefficient of variation (`cv`) from which `nmin` is derived.
  May be grouped with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- rmax:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Maximum annual net productivity rate \\R\_{max}\\, i.e.
  \\\lambda\_{max} - 1\\ (e.g. `0.04` for many large mammals, `0.12` for
  fast-breeding species). An instantaneous rate `r` can be converted
  with `rmax = exp(r) - 1`.

- removal:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual human-caused removal (offtake), on the same basis as
  the population size.

- nmin:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Minimum population estimate. Optional if `n` (and optionally `cv`) are
  supplied.

- n:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Point abundance estimate, used with `cv` to derive `nmin`. Use the
  current abundance, not the carrying capacity.

- cv:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Coefficient of variation of `n` (e.g. `0.3` for a 30% CV). Treated as
  0 (so `nmin = n`) when omitted.

- fr:

  Recovery factor \\F_R\\ in `[0.1, 1]` (default `0.5`). May be a single
  value or an embraced column. Use 0.1 for endangered or poorly known
  populations.

- percentile:

  Lower percentile used to convert `n`/`cv` to `nmin` (default `0.20`,
  giving `z = 0.842`).

- uncertainty:

  Optional named list of coefficients of variation (numbers or columns
  of `data`). Allowed names: `n` or `nmin` (whichever is used), `rmax`
  and `removal`. If `cv` already describes the uncertainty of `n`, do
  not add `n` here as well.

- n_sim:

  Number of draws used when `uncertainty` is given (default `1000`).

- level:

  Width of the interval `limit_lo` to `limit_hi` (default `0.95`).

- seed:

  Optional integer seed, for reproducible draws.

## Value

An offtake tibble with one row per input row and the columns:

- nmin:

  Minimum population estimate used (supplied, or derived from `n` and
  `cv`).

- rmax:

  Maximum annual net productivity rate used.

- fr:

  Recovery factor applied.

- limit:

  Potential biological removal, the maximum sustainable annual removal.

- observed:

  The observed annual removal (`removal`).

- ratio:

  Exploitation ratio, `observed / limit`.

- sustainable:

  Logical verdict: `TRUE` when `observed <= limit`.

- limit_lo, limit_hi, p_unsustainable:

  Only with `uncertainty`; see the Uncertainty section.

Group columns come first when `data` is grouped.

## Details

\$\$PBR = N\_{min}\\ \tfrac{1}{2} R\_{max}\\ F_R\$\$ where `N_min` is a
conservative (minimum) population estimate, \\R\_{max}\\ is the maximum
annual net productivity rate (the discrete rate \\\lambda\_{max} - 1\\)
and \\F_R\\ is a recovery factor between 0.1 and

1.  The term \\\tfrac{1}{2} R\_{max}\\ is the per-capita growth rate of
    a logistic population at half its carrying capacity, where it is
    most productive. **A removal exceeding PBR indicates
    unsustainability.**

`N_min` may be supplied directly through `nmin`, or computed from an
abundance estimate `n` and its coefficient of variation `cv` as the
lower tail of a log-normal distribution, \\N\_{min} = N /
\exp\\\big(z\sqrt{\ln(1+CV^2)}\big)\\, with `z` the standard-normal
quantile for the chosen percentile (`z = 0.842` for the 20th percentile,
as recommended by Wade 1998).

Grouped data frames (from
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))
are accepted; the group columns are kept in the output.

## Uncertainty

Model inputs are rarely known precisely: \\\lambda\_{max}\\ from Cole's
equation and from the Caughley & Krebs (1983) equation can differ by up
to four times for the same species, and densities depend on the survey
method (Adounke et al. 2026). With `uncertainty`, each named input is
drawn `n_sim` times from a lognormal distribution with the given value
as its mean and the given coefficient of variation (CV), the safe limit
is recomputed for every draw, and the output gains three columns:
`limit_lo` and `limit_hi` (the central `level` interval of the limit)
and `p_unsustainable` (the share of draws in which the observed offtake
exceeds the limit). This turns the TRUE/FALSE verdict into a
probability. For `lambda`, the CV applies to the annual surplus
\\\lambda\_{max} - 1\\, so that drawn growth rates stay above 1.

## References

Wade, P. R. (1998) Calculating limits to the allowable human-caused
mortality of cetaceans and pinnipeds. *Marine Mammal Science* 14, 1-37.
[doi:10.1111/j.1748-7692.1998.tb00688.x](https://doi.org/10.1111/j.1748-7692.1998.tb00688.x)

## See also

[`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md),
[`ot_samse()`](https://stangandaho.github.io/offtake/reference/ot_samse.md),
[`ot_msy()`](https://stangandaho.github.io/offtake/reference/ot_msy.md)

## Examples

``` r
# One row per stock. Columns:
#   abund = point abundance estimate (animals)
#   cv = coefficient of variation of that estimate
#   rmax = maximum annual net productivity rate
#   take = observed annual removal (animals)
d <- data.frame(
  stock = c("A", "B"),
  abund = c(1200, 800),
  cv = c(0.3, 0.2),
  rmax = c(0.04, 0.12),
  take = c(15, 40)
)
ot_pbr(d, rmax = rmax, removal = take, n = abund, cv = cv, fr = 0.5)
#> model-based: PBR (potential biological removal)
#> 
#> # A tibble: 2 × 7
#>    nmin  rmax    fr limit observed ratio sustainable
#> * <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl> <lgl>      
#> 1  937.  0.04   0.5  9.37       15  1.60 FALSE      
#> 2  677.  0.12   0.5 20.3        40  1.97 FALSE      
```
