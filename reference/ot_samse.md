# Sustainable anthropogenic mortality in stochastic environments (SAMSE)

Estimates, by Monte-Carlo simulation, the largest constant annual
removal that does **not** drive a negative long-run stochastic growth
rate, given environmental variability in the population growth rate, and
compares it with the observed removal. SAMSE is a stochastic successor
to
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
  years = 50,
  nsims = 500,
  tol = 0.001,
  seed = NULL
)
```

## Arguments

- data:

  A data frame with **one row per population/stock**. Required columns:
  mean maximum growth rate (`rmax`), its environmental standard
  deviation (`sd_env`), the observed annual removal (`removal`) and the
  starting population size (`n0`).

- rmax:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Mean maximum annual growth rate \\r\_{max}\\ on the log scale (so
  \\\lambda = e^{r\_{max}}\\; e.g. `0.10` for ~10% growth).

- sd_env:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Environmental standard deviation of the annual growth rate
  \\\sigma_e\\ (larger = more year-to-year variability, lower
  sustainable removal).

- removal:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Observed annual human-caused removal (animals per year).

- n0:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Starting (ideally minimum, `N_min`-style) population size.

- years:

  Projection horizon in years (default `50`).

- nsims:

  Number of Monte-Carlo trajectories per candidate removal (default
  `500`).

- tol:

  Convergence tolerance of the bisection, as a fraction of `n0` (default
  `0.001`).

- seed:

  Optional integer seed for reproducibility (uses common random numbers
  across candidate removals).

## Value

An offtake tibble with **one row per input row** and the columns:

- n0:

  Starting population size used.

- rmax:

  Mean maximum growth rate used.

- sd_env:

  Environmental standard deviation used.

- samse_limit:

  Estimated SAMSE limit – the largest constant annual removal (animals
  per year) keeping the stochastic growth rate non-negative.

- removal:

  The observed annual removal, echoed back for comparison.

- sustainable:

  Logical verdict: `TRUE` when `removal <= samse_limit`.

## Implementation note

The original SAMSE limit of Manlik et al. (2022) is obtained with an
individual-based population viability analysis run in the *Vortex*
software, iterating removal levels until the forecast stochastic growth
rate is no longer negative. This function reproduces that **principle**
with a transparent, self-contained stochastic projection; it is **not**
a re-implementation of Vortex and does not include age structure,
inbreeding or catastrophes. Use it as a precautionary screening tool and
cite Manlik et al. (2022) for the concept.

The projection is a stochastic exponential model with a constant annual
take `H`: \$\$N\_{t+1} = \max\\\big(0,\\ N_t\\ e^{r_t} - H\big), \quad
r_t \sim \mathrm{Normal}(r\_{max},\\ \sigma_e).\$\$ The stochastic
growth rate is \\\rho(H) = \mathrm{mean}\\ \log(N\_{t}/N\_{0})/T\\. The
SAMSE limit is the largest `H` with \\\rho(H) \ge 0\\, found by
bisection. **An observed removal above the SAMSE limit indicates
unsustainability.**

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
#   rmax = mean max growth rate (log scale, per year)
#   sd = environmental SD of the annual growth rate
#   n0 = starting population size
#   take = observed annual removal (animals)
d <- data.frame(stock = "A", rmax = 0.10, sd = 0.25, n0 = 500, take = 20)
ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0,
         nsims = 200, years = 40, seed = 1)
#> <offtake: SAMSE (sustainable anthropogenic mortality, stochastic)>  (model-based)
#> Reference: Manlik et al. (2022) 
#> Note: Simplified stochastic re-implementation of the SAMSE principle; not the original Vortex-based procedure. 
#> 
#> # A tibble: 1 × 6
#>      n0  rmax sd_env samse_limit removal sustainable
#> * <dbl> <dbl>  <dbl>       <dbl>   <dbl> <lgl>      
#> 1   500   0.1   0.25        17.7      20 FALSE      
```
