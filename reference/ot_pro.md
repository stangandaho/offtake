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
  w = NULL
)
```

## Arguments

- data:

  A data frame with **one row per species or population unit**. Required
  columns: carrying capacity (`k`) and observed annual offtake
  (`harvest`), on the *same basis* (both per km^2, or both absolute
  counts). You must also supply the growth rate – either a `lambda`
  column or the life-history columns `b`, `a`, `w` – and the mortality
  factor – either an `f` column or a `longevity` column.

- k:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Carrying capacity `K` (density per km^2, or absolute population size).

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

## Value

An offtake tibble with **one row per input row** and the columns:

- lambda_max:

  Maximum finite rate of increase used, whether supplied or computed
  from `b`, `a`, `w`.

- f:

  Mortality factor `F` used (0.2, 0.4 or 0.6).

- production:

  Estimated maximum sustainable harvest `P`, on the same basis as `k`
  and `harvest` (e.g. individuals per km^2 per year).

- harvest:

  The observed annual offtake, echoed back for comparison.

- sustainable:

  Logical verdict: `TRUE` when `harvest <= production`, `FALSE` when the
  observed offtake exceeds the estimated production.

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
#> <offtake: Pro (Robinson & Redford production model)>  (model-based)
#> Reference: Robinson & Redford (1991) 
#> 
#> # A tibble: 2 × 5
#>   lambda_max     f production harvest sustainable
#> *      <dbl> <dbl>      <dbl>   <dbl> <lgl>      
#> 1       1.35   0.4       0.84       3 FALSE      
#> 2       1.55   0.4       3.3        6 FALSE      
```
