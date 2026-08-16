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
  percentile = 0.2
)
```

## Arguments

- data:

  A data frame with **one row per population/stock**. Required columns:
  the maximum growth rate (`rmax`) and the observed annual removal
  (`removal`). You must also supply the population size, either as a
  minimum estimate (`nmin`) or as a point estimate (`n`) with an
  optional coefficient of variation (`cv`) from which `nmin` is derived.

- rmax:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Maximum per-capita annual rate of increase \\R\_{max}\\ (e.g. `0.04`
  for many large mammals, `0.12` for fast-breeding species).

- removal:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual human-caused removal (offtake), in numbers of animals.

- nmin:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Minimum population estimate. Optional if `n` (and optionally `cv`) are
  supplied.

- n:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Point abundance estimate, used with `cv` to derive `nmin`.

- cv:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Coefficient of variation of `n` (e.g. `0.3` for a 30% CV). Treated as
  0 (so `nmin = n`) when omitted.

- fr:

  Recovery factor \\F_R\\ in `[0.1, 1]` (default `0.5`). May be a single
  value or an embraced column.

- percentile:

  Lower percentile used to convert `n`/`cv` to `nmin` (default `0.20`,
  giving `z = 0.842`).

## Value

An offtake tibble with **one row per input row** and the columns:

- nmin:

  Minimum population estimate used (supplied, or derived from `n` and
  `cv`).

- rmax:

  Maximum per-capita rate of increase used.

- fr:

  Recovery factor applied.

- pbr:

  Potential biological removal – the maximum sustainable annual removal,
  in numbers of animals.

- removal:

  The observed annual removal, echoed back for comparison.

- sustainable:

  Logical verdict: `TRUE` when `removal <= pbr`, `FALSE` when the
  observed removal exceeds PBR.

## Details

\$\$PBR = N\_{min}\\ \tfrac{1}{2} R\_{max}\\ F_R\$\$ where `N_min` is a
conservative (minimum) population estimate, \\R\_{max}\\ is the maximum
per-capita annual rate of increase and \\F_R\\ is a recovery factor
between 0.1 and 1. **A removal exceeding PBR indicates
unsustainability.**

`N_min` may be supplied directly through `nmin`, or computed from an
abundance estimate `n` and its coefficient of variation `cv` as the
lower tail of a log-normal distribution, \\N\_{min} = N /
\exp\\\big(z\sqrt{\ln(1+CV^2)}\big)\\, with `z` the standard-normal
quantile for the chosen percentile (`z = 0.842` for the 20th percentile,
as recommended by Wade 1998).

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
#   rmax = maximum per-capita rate of increase (per year)
#   take = observed annual removal (animals)
d <- data.frame(
  stock = c("A", "B"),
  abund = c(1200, 800),
  cv = c(0.3, 0.2),
  rmax = c(0.04, 0.12),
  take = c(15, 40)
)
ot_pbr(d, rmax = rmax, removal = take, n = abund, cv = cv, fr = 0.5)
#> <offtake: PBR (potential biological removal)>  (model-based)
#> Reference: Wade (1998) 
#> 
#> # A tibble: 2 × 6
#>    nmin  rmax    fr   pbr removal sustainable
#> * <dbl> <dbl> <dbl> <dbl>   <dbl> <lgl>      
#> 1  937.  0.04   0.5  9.37      15 FALSE      
#> 2  677.  0.12   0.5 20.3       40 FALSE      
```
