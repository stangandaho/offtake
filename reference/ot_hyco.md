# Hunting yield comparison (HYCo)

Compares harvested biomass (or catch-per-unit-effort, CPUE) between
more- and less-hunted sites. Lower yields at the more-hunted site(s) are
interpreted as unsustainable (Adounke et al. 2026; Weinbaum et al.
2013). If an `effort` column is supplied the metric compared is CPUE
(`yield / effort`); otherwise raw yield is compared.

## Usage

``` r
ot_hyco(data, yield, group, reference, effort = NULL, alpha = 0.05)
```

## Arguments

- data:

  A data frame in *long* format with one row per harvest record. It must
  contain at least a yield column (`yield`) and a site column (`group`);
  optionally a hunting-effort column (`effort`).

- yield:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column of harvested biomass (e.g. kg) or number of animals taken.

- group:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column identifying the site / hunting intensity.

- reference:

  Character scalar naming the level of `group` treated as the
  less-hunted reference (e.g. `"low"`).

- effort:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional column of hunting effort (e.g. hunter-days). When supplied,
  `yield / effort` (CPUE) is compared instead of raw yield.

- alpha:

  Significance level for the one-sided test (default `0.05`).

## Value

An offtake tibble with **one row per hunted site**. The columns are the
same as for
[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md),
except that `hunted_value` and `reference_value` hold the mean yield (or
mean CPUE when `effort` is supplied) at each site rather than density;
`sustainable` is `FALSE` when the more-hunted site yields significantly
less.

## References

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.
[doi:10.1016/j.biocon.2025.111606](https://doi.org/10.1016/j.biocon.2025.111606)

Weinbaum, K. Z., Brashares, J. S., Golden, C. D. & Getz, W. M. (2013)
Searching for sustainability: are assessments of wildlife harvests
behind the times? *Ecology Letters* 16, 99-111.
[doi:10.1111/ele.12008](https://doi.org/10.1111/ele.12008)

## See also

[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md),
[`ot_asc()`](https://stangandaho.github.io/offtake/reference/ot_asc.md)

## Examples

``` r
# One row per harvest record. `site` = intensity, `kg` = biomass taken,
# `days` = hunter-days of effort.
d <- data.frame(
  site = rep(c("low", "high"), each = 3),  # less- vs more-hunted
  kg = c(40, 45, 38, 20, 25, 18), # harvested biomass (kg)
  days = c(10, 11, 9, 10, 12, 9) # hunting effort (hunter-days)
)
# Compare CPUE (kg per hunter-day):
ot_hyco(d, yield = kg, group = site, reference = "low", effort = days)
#> <offtake: HYCo (hunting yield comparison, CPUE)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> 
#> # A tibble: 1 × 7
#>   hunted_site reference_site hunted_value reference_value pct_change   p_value
#> * <chr>       <chr>                 <dbl>           <dbl>      <dbl>     <dbl>
#> 1 high        low                    2.03            4.10      -50.6 0.0000864
#> # ℹ 1 more variable: sustainable <lgl>
```
