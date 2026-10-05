# Hunting yield comparison (HYCo)

Compares harvested biomass (or catch-per-unit-effort, CPUE) between
more- and less-hunted sites. Lower yields at the more-hunted site(s) are
interpreted as unsustainable (Adounke et al. 2026; Weinbaum et al.
2013). If an `effort` column is supplied the metric compared is CPUE
(`yield / effort`, one value per record); otherwise raw yield is
compared. The statistics, the outcome and the caveats are those of
[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md).

## Usage

``` r
ot_hyco(
  data,
  yield,
  group,
  reference,
  effort = NULL,
  alpha = 0.05,
  p_adjust = "holm",
  decline = 0.2
)
```

## Arguments

- data:

  A data frame in *long* format with one row per harvest record. It must
  contain at least a yield column (`yield`) and a site column (`group`);
  optionally a hunting-effort column (`effort`). May be grouped with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

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
  Optional column of hunting effort (e.g. hunter-days), strictly
  positive. When supplied, `yield / effort` (CPUE) is compared instead
  of raw yield.

- alpha:

  Significance level of the one-sided test; the confidence intervals
  have level `1 - alpha` (default `0.05`).

- p_adjust:

  Method used to adjust p-values when several hunted sites are compared
  with the same reference, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) (default
  `"holm"`; `"none"` for no adjustment).

- decline:

  Smallest relative decline considered substantial, between 0 and 1
  (default `0.2`, i.e. 20%). Used to tell "no substantial decline" from
  "inconclusive".

## Value

An offtake tibble with one row per hunted site. The columns are the same
as for
[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md),
except that `hunted_value` and `reference_value` hold the mean yield (or
mean CPUE when `effort` is supplied) at each site rather than density.

## Outcome

The verdict has three levels, so that "no difference found" is not
mistaken for evidence of sustainability:

- `"unsustainable"`: the hunted site is significantly lower (adjusted
  p-value below `alpha`); `sustainable = FALSE`.

- `"no substantial decline"`: not significant, and the confidence
  interval of the response ratio rules out a decline of `decline`
  (default 20%) or more; `sustainable = TRUE`.

- `"inconclusive"`: no test was possible (fewer than two replicates per
  site) or the interval is too wide to rule out a substantial decline;
  `sustainable = NA`.

## Pseudoreplication

Replicates taken inside one hunted site and one reference site
(transects, camera stations) are not independent replicates of hunting.
The test then shows that the two *sites* differ, not that *hunting*
causes the difference (Hurlbert 1984). Prefer several hunted and several
reference sites, keep the sites otherwise comparable, and read the
result as a warning sign rather than a proof.

Grouped data frames (from
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html))
are analysed group by group, for example one comparison per species, and
the group columns are kept in the output.

## References

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.
[doi:10.1016/j.biocon.2025.111606](https://doi.org/10.1016/j.biocon.2025.111606)

Weinbaum, K. Z., Brashares, J. S., Golden, C. D. & Getz, W. M. (2013)
Searching for sustainability: are assessments of wildlife harvests
behind the times? *Ecology Letters* 16, 99-111.
[doi:10.1111/ele.12008](https://doi.org/10.1111/ele.12008)

Hedges, L. V., Gurevitch, J. & Curtis, P. S. (1999) The meta-analysis of
response ratios in experimental ecology. *Ecology* 80, 1150-1156.

Hurlbert, S. H. (1984) Pseudoreplication and the design of ecological
field experiments. *Ecological Monographs* 54, 187-211.

## See also

[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md),
[`ot_asc()`](https://stangandaho.github.io/offtake/reference/ot_asc.md)

## Examples

``` r
# One row per harvest record. `site` = intensity, `kg` = biomass taken,
# `days` = hunter-days of effort.
d <- data.frame(
  site = rep(c("low", "high"), each = 3), # less- vs more-hunted
  kg = c(40, 45, 38, 20, 25, 18), # harvested biomass (kg)
  days = c(10, 11, 9, 10, 12, 9) # hunting effort (hunter-days)
)
# Compare CPUE (kg per hunter-day):
ot_hyco(d, yield = kg, group = site, reference = "low", effort = days)
#> index-based: HYCo (hunting yield comparison, CPUE)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference hunted_value reference_value
#> * <chr>       <chr>             <int>       <int>        <dbl>           <dbl>
#> 1 high        low                   3           3         2.03            4.10
#> # ℹ 8 more variables: pct_change <dbl>, lnrr <dbl>, lnrr_lo <dbl>,
#> #   lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>, sustainable <lgl>
```
