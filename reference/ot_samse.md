# Sustainable anthropogenic mortality in stochastic environments (SAMSE)

Estimates, by Monte-Carlo simulation, the largest constant annual
removal that does **not** produce a negative stochastic growth rate from
the current population size, given environmental variability in the
growth rate and density dependence, and compares it with the observed
removal. SAMSE is a stochastic successor to
[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md)
proposed by Manlik et al. (2022) and highlighted by Adounke et al.
(2026).

## Usage

``` r
ot_samse(
  data,
  rmax,
  sd_env,
  removal,
  n0,
  k,
  theta = 1,
  years = 50,
  nsims = 500,
  max_extinction = 0.05,
  tol = 0.001,
  seed = NULL
)
```

## Arguments

- data:

  A data frame with **one row per population/stock**. Required columns:
  mean maximum growth rate (`rmax`), its environmental standard
  deviation (`sd_env`), the observed annual removal (`removal`), the
  current population size (`n0`) and the carrying capacity (`k`). May be
  grouped with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- rmax:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Maximum annual growth rate \\r\_{max}\\ on the log scale (growth of a
  small population is about \\e^{r\_{max}}\\ per year; e.g. `0.10` for
  ~10% growth).

- sd_env:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Environmental standard deviation of the annual growth rate
  \\\sigma_e\\ (larger = more year-to-year variability, lower
  sustainable removal).

- removal:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual human-caused removal, on the same basis as `n0` and
  `k`.

- n0:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Current population size (or density), the starting point of the
  projections.

- k:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Carrying capacity, on the same basis as `n0`.

- theta:

  Shape of density dependence \\\theta\\ (default `1`, logistic; values
  above 1 keep growth high until the population is close to `K`).

- years:

  Projection horizon in years (default `50`).

- nsims:

  Number of Monte-Carlo trajectories per candidate removal (default
  `500`).

- max_extinction:

  Largest acceptable probability that a trajectory goes extinct within
  `years` (default `0.05`).

- tol:

  Convergence tolerance of the bisection, as a fraction of `n0` (default
  `0.001`).

- seed:

  Optional integer seed for reproducibility (common random numbers are
  used across candidate removals).

## Value

An offtake tibble with one row per input row and the columns:

- n0:

  Current population size used.

- k:

  Carrying capacity used.

- rmax:

  Maximum growth rate used.

- sd_env:

  Environmental standard deviation used.

- limit:

  Estimated SAMSE limit: the largest constant annual removal keeping the
  stochastic growth rate non-negative and the extinction probability at
  or below `max_extinction`.

- observed:

  The observed annual removal (`removal`).

- ratio:

  Exploitation ratio, `observed / limit`.

- sustainable:

  Logical verdict: `TRUE` when `observed <= limit`.

- p_extinct:

  Probability of extinction within `years` if the observed removal
  continues.

Group columns come first when `data` is grouped.

## Model

The population is projected with a theta-logistic (Ricker type) model
with environmental noise and a constant annual take `H`: \$\$N\_{t+1} =
\max\\\Big(0,\\ N_t \exp\\\big\[r\_{max}\big(1 - (N_t/K)^{\theta}\big) +
\varepsilon_t\big\] - H\Big), \quad \varepsilon_t \sim
\mathrm{Normal}(0, \sigma_e).\$\$ Growth slows as the population
approaches its carrying capacity `K`, so a population near `K` has
little surplus to give. For each candidate `H`, `nsims` trajectories of
`years` years are run from `n0`. The stochastic growth rate is the mean
of the annual \\\log(N\_{t+1}/N_t)\\ over all trajectories and years in
which the population is still present, and the extinction probability is
the share of trajectories that reach zero. The SAMSE limit is the
largest `H` for which the stochastic growth rate is not negative **and**
the extinction probability does not exceed `max_extinction`, found by
bisection with common random numbers. **An observed removal above the
SAMSE limit indicates unsustainability.**

Because the criterion is "no decline from the current size", the limit
depends on `n0`: it is close to the surplus the population produces at
that size, reduced by environmental variability. Use the current
abundance for `n0`; a population at its carrying capacity has (almost)
no surplus.

## Implementation note

Manlik et al. (2022) obtain the SAMSE limit with an individual-based
population viability analysis in the *Vortex* software. This function
reproduces the principle (largest removal without a negative stochastic
growth rate) with a transparent count-based projection. It does not
include age structure, demographic stochasticity, inbreeding or
catastrophes, and the extinction constraint (`max_extinction`) is an
addition of this package.

Grouped data frames (from
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))
are accepted; the group columns are kept in the output.

## References

Manlik, O., Lacy, R. C., Sherwin, W. B., Finn, H., Loneragan, N. R. &
Allen, S. J. (2022) A stochastic model for estimating sustainable limits
to wildlife mortality in a changing world. *Conservation Biology* 36,
e13897. [doi:10.1111/cobi.13897](https://doi.org/10.1111/cobi.13897)

## See also

[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md)

## Examples

``` r
# One row per stock. Columns:
#   rmax = maximum growth rate (log scale, per year)
#   sd = environmental SD of the annual growth rate
#   n0 = current population size
#   K = carrying capacity
#   take = observed annual removal (animals)
d <- data.frame(stock = "A", rmax = 0.10, sd = 0.25, n0 = 500, K = 1000, take = 20)
ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
         nsims = 200, years = 40, seed = 1)
#> model-based: SAMSE (sustainable anthropogenic mortality, stochastic)
#> Note: Count-based re-implementation of the SAMSE principle with density dependence; not the original Vortex-based procedure. 
#> 
#> # A tibble: 1 × 9
#>      n0     k  rmax sd_env limit observed ratio sustainable p_extinct
#> * <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl> <dbl> <lgl>           <dbl>
#> 1   500  1000   0.1   0.25  7.81       20  2.56 FALSE           0.305
```
