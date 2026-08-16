# Maximum sustainable yield of the logistic model (MSY)

Estimates the maximum sustainable yield from the logistic (surplus
production) model and compares it with the observed annual harvest. This
is the stock-recruitment / surplus-production benchmark used across
fisheries and, increasingly, wildlife harvest (Adounke et al. 2026).

## Usage

``` r
ot_msy(data, r, k, harvest, n = NULL)
```

## Arguments

- data:

  A data frame with **one row per population/stock**. Required columns:
  intrinsic growth rate (`r`), carrying capacity (`k`) and observed
  annual harvest (`harvest`); optionally current abundance (`n`).

- r:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Intrinsic rate of natural increase (per year).

- k:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Carrying capacity `K` (population size or density).

- harvest:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual harvest, on the same basis as `k`.

- n:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional current abundance, used to compute `surplus_at_n`.

## Value

An offtake tibble with **one row per input row** and the columns:

- r:

  Intrinsic growth rate used.

- k:

  Carrying capacity used.

- msy:

  Maximum sustainable yield, `(r * k) / 4`.

- harvest:

  The observed annual harvest, echoed back for comparison.

- sustainable:

  Logical verdict: `TRUE` when `harvest <= msy`.

- surplus_at_n:

  Only present when `n` is supplied: the sustainable yield at the
  *current* abundance, `(r * n) * (1 - n / k)`. Equals `msy` when
  `n = k / 2`.

## Details

The logistic model of population growth is \$\$\frac{dN}{dt} = r N
\left(1 - \frac{N}{K}\right)\$\$ whose surplus production is maximised
at \\N = K/2\\, giving \$\$MSY = \frac{rK}{4}.\$\$ **A harvest above MSY
is considered unsustainable.** When a current abundance `n` is supplied,
the instantaneous surplus production at that abundance, \\rN(1 - N/K)\\,
is also returned as `surplus_at_n`, which is the sustainable yield at
the *current* (not optimal) population size.

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
#> <offtake: MSY (logistic maximum sustainable yield)>  (model-based)
#> Reference: Schaefer (1954) 
#> 
#> # A tibble: 1 × 6
#>       r     k   msy harvest sustainable surplus_at_n
#> * <dbl> <dbl> <dbl>   <dbl> <lgl>              <dbl>
#> 1   0.4  1000   100      80 TRUE                  96
```
