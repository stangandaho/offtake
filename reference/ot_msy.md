# Maximum sustainable yield of the logistic model (MSY)

Estimates the maximum sustainable yield from the logistic (surplus
production) model and compares it with the observed annual harvest. This
is the stock-recruitment / surplus-production benchmark used across
fisheries and, increasingly, wildlife harvest (Adounke et al. 2026).

## Usage

``` r
ot_msy(
  data,
  r,
  k,
  harvest,
  n = NULL,
  uncertainty = NULL,
  n_sim = 1000,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- data:

  A data frame with **one row per population/stock**. Required columns:
  intrinsic growth rate (`r`), carrying capacity (`k`) and observed
  annual harvest (`harvest`); optionally current abundance (`n`). May be
  grouped with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- r:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Intrinsic rate of natural increase (per year), strictly positive.

- k:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Carrying capacity `K` (population size or density), strictly positive.

- harvest:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual harvest, on the same basis as `k`.

- n:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional current abundance, used to compute `surplus_at_n`.

- uncertainty:

  Optional named list of coefficients of variation (numbers or columns
  of `data`). Allowed names: `r`, `k` and `harvest`.

- n_sim:

  Number of draws used when `uncertainty` is given (default `1000`).

- level:

  Width of the interval `limit_lo` to `limit_hi` (default `0.95`).

- seed:

  Optional integer seed, for reproducible draws.

## Value

An offtake tibble with one row per input row and the columns:

- r:

  Intrinsic growth rate used.

- k:

  Carrying capacity used.

- limit:

  Maximum sustainable yield, `(r * k) / 4`.

- observed:

  The observed annual harvest (`harvest`).

- ratio:

  Exploitation ratio, `observed / limit`.

- sustainable:

  Logical verdict: `TRUE` when `observed <= limit`.

- limit_lo, limit_hi, p_unsustainable:

  Only with `uncertainty`; see the Uncertainty section.

- surplus_at_n:

  Only present when `n` is supplied: the sustainable yield at the
  *current* abundance, `(r * n) * (1 - n / k)`. Equals `limit` when
  `n = k / 2`.

Group columns come first when `data` is grouped.

## Details

The logistic model of population growth is \$\$\frac{dN}{dt} = r N
\left(1 - \frac{N}{K}\right)\$\$ whose surplus production is maximised
at \\N = K/2\\, giving \$\$MSY = \frac{rK}{4}.\$\$ **A harvest above MSY
is considered unsustainable.** When a current abundance `n` is supplied,
the instantaneous surplus production at that abundance, \\rN(1 - N/K)\\,
is also returned as `surplus_at_n`, which is the sustainable yield at
the *current* (not optimal) population size. A population already below
\\K/2\\ cannot sustain MSY.

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

Schaefer, M. B. (1954) Some aspects of the dynamics of populations
important to the management of the commercial marine fisheries.
*Bulletin of the Inter-American Tropical Tuna Commission* 1, 27-56.

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.
[doi:10.1016/j.biocon.2025.111606](https://doi.org/10.1016/j.biocon.2025.111606)

## See also

[`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md),
[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md)

## Examples

``` r
# One row per stock. Columns:
#   r = intrinsic growth rate (per year)
#   K = carrying capacity
#   take = observed annual harvest
#   now = current abundance (optional, for surplus_at_n)
d <- data.frame(stock = "A", r = 0.4, K = 1000, take = 80, now = 600)
ot_msy(d, r = r, k = K, harvest = take, n = now)
#> model-based: MSY (logistic maximum sustainable yield)
#> 
#> # A tibble: 1 × 7
#>       r     k limit observed ratio sustainable surplus_at_n
#> * <dbl> <dbl> <dbl>    <dbl> <dbl> <lgl>              <dbl>
#> 1   0.4  1000   100       80   0.8 TRUE                  96
```
