# Robinson & Redford production model (Pro)

Estimates the maximum sustainable harvest (production) of a population
from its carrying capacity, maximum rate of increase and a
longevity-based mortality factor, and compares it with the observed
annual offtake. This is the most widely used model in bushmeat
sustainability assessments (Adounke et al. 2026).

## Usage

``` r
ot_pro(
  data,
  k,
  harvest,
  lambda = NULL,
  f = NULL,
  longevity = NULL,
  b = NULL,
  a = NULL,
  w = NULL,
  uncertainty = NULL,
  n_sim = 1000,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  A data frame with **one row per species or population unit**. Required
  columns: carrying capacity (`k`) and observed annual offtake
  (`harvest`), on the *same basis* (both per km^2, or both absolute
  counts). You must also supply the growth rate – either a `lambda`
  column or the life-history columns `b`, `a`, `w` – and the mortality
  factor – either an `f` column or a `longevity` column. May be grouped
  with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- k:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Carrying capacity `K` (density per km^2, or absolute population size),
  strictly positive.

- harvest:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual offtake, on the same basis as `k`.

- lambda:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Maximum finite rate of increase \\\lambda\_{max}\\ (dimensionless, \>
  1). Optional if `b`, `a`, `w` are given.

- f:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Mortality factor `F` (0.2, 0.4 or 0.6). Optional if `longevity` is
  given.

- longevity:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Age at last reproduction / lifespan (years), used to derive `F` when
  `f` is not given.

- b, a, w:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Life-history columns passed to
  [`ot_lambda_max()`](https://stangandaho.github.io/offtake/reference/ot_lambda_max.md)
  when `lambda` is not given: annual female offspring per female (`b`),
  age at first reproduction (`a`) and age at last reproduction (`w`),
  all in years.

- uncertainty:

  Optional named list of coefficients of variation, one per uncertain
  input; each value can be a number or a column of `data`. Here the
  allowed names are `k`, `lambda` and `harvest`, e.g.
  `list(k = 0.3, lambda = 0.2)`. See the Uncertainty section.

- n_sim:

  Number of draws used when `uncertainty` is given (default `1000`).

- level:

  Width of the interval `limit_lo` to `limit_hi` (default `0.95`).

- seed:

  Optional integer seed, for reproducible draws.

## Value

An offtake tibble with **one row per input row** and the columns:

- lambda_max:

  Maximum finite rate of increase used, whether supplied or computed
  from `b`, `a`, `w`.

- f:

  Mortality factor `F` used (0.2, 0.4 or 0.6).

- limit:

  Estimated production `P`, the maximum sustainable harvest, on the same
  basis as `k` and `harvest` (e.g. individuals per km^2 per year).

- observed:

  The observed annual offtake (`harvest`).

- ratio:

  Exploitation ratio, `observed / limit`. Above 1 means the offtake
  exceeds the production; 3 means three times too much.

- sustainable:

  Logical verdict: `TRUE` when `observed <= limit`.

- limit_lo, limit_hi, p_unsustainable:

  Only with `uncertainty`: the interval of the limit and the probability
  that the offtake exceeds it.

Group columns come first when `data` is grouped.

## Details

The production (maximum sustainable number that can be taken per year)
is \$\$P = 0.6\\K\\(\lambda\_{max} - 1)\\F\$\$ where `K` is
carrying-capacity density (or population size), \\\lambda\_{max}\\ is
the maximum finite rate of increase and `F` is the mortality factor set
by longevity: `F = 0.6` for short-lived species (lifespan ~ 5 years),
`0.4` for medium-lived (5-10 years) and `0.2` for long-lived species (\>
10 years) (Robinson & Redford 1991). The factor `0.6 K` is the density
at which production is assumed maximal. **If the observed harvest
exceeds `P` the harvest is considered unsustainable.**

Supply \\\lambda\_{max}\\ through `lambda` directly, or leave it `NULL`
and provide the life-history columns `b`, `a`, `w` so it is computed
with
[`ot_lambda_max()`](https://stangandaho.github.io/offtake/reference/ot_lambda_max.md)
(Cole's equation). Provide the mortality factor either through `f`
directly or through a `longevity` column.

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

Robinson, J. G. & Redford, K. H. (1991) Sustainable harvest of
Neotropical forest mammals. In *Neotropical Wildlife Use and
Conservation* (eds Robinson & Redford), 415-429. University of Chicago
Press.

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.
[doi:10.1016/j.biocon.2025.111606](https://doi.org/10.1016/j.biocon.2025.111606)

## See also

[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md),
[`ot_msy()`](https://stangandaho.github.io/offtake/reference/ot_msy.md),
[`ot_lambda_max()`](https://stangandaho.github.io/offtake/reference/ot_lambda_max.md)

## Examples

``` r
# One row per species. Columns:
#   dens = carrying-capacity density (ind / km^2)
#   offtake = observed annual harvest (ind / km^2 / year)
#   lifespan = age at last reproduction (years) -> sets F
#   lam = maximum finite rate of increase (lambda_max)
d <- data.frame(
  species = c("red duiker", "blue duiker"),
  dens = c(10, 25),
  offtake = c(3.0, 6.0),
  lifespan = c(9, 6),
  lam = c(1.35, 1.55)
)
ot_pro(d, k = dens, harvest = offtake, lambda = lam, longevity = lifespan)
#> model-based: Pro (Robinson & Redford production model)
#> 
#> # A tibble: 2 × 6
#>   lambda_max     f limit observed ratio sustainable
#> *      <dbl> <dbl> <dbl>    <dbl> <dbl> <lgl>      
#> 1       1.35   0.4  0.84        3  3.57 FALSE      
#> 2       1.55   0.4  3.3         6  1.82 FALSE      

# With 30% uncertainty on K and 20% on the growth surplus
ot_pro(d, k = dens, harvest = offtake, lambda = lam, longevity = lifespan,
       uncertainty = list(k = 0.3, lambda = 0.2), seed = 1)
#> model-based: Pro (Robinson & Redford production model)
#> 
#> # A tibble: 2 × 9
#>   lambda_max     f limit observed ratio sustainable limit_lo limit_hi
#> *      <dbl> <dbl> <dbl>    <dbl> <dbl> <lgl>          <dbl>    <dbl>
#> 1       1.35   0.4  0.84        3  3.57 FALSE          0.372     1.60
#> 2       1.55   0.4  3.3         6  1.82 FALSE          1.49      6.23
#> # ℹ 1 more variable: p_unsustainable <dbl>
```
