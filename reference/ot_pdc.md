# Population density comparison (PDC)

Compares wildlife abundance/density between hunted and reference
(unhunted or lightly hunted) sites. Lower density at the hunted site(s)
is interpreted as unsustainable (Adounke et al. 2026; Weinbaum et al.
2013).

## Usage

``` r
ot_pdc(
  data,
  density,
  group,
  reference,
  alpha = 0.05,
  p_adjust = "holm",
  decline = 0.2
)
```

## Arguments

- data:

  A data frame in *long* format with one row per density estimate. It
  must contain at least: a numeric density column (`density`) and a
  site/treatment column (`group`) whose values include the reference
  level and one or more hunted levels. Extra columns are ignored. May be
  grouped with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- density:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column of abundance or density estimates (e.g. individuals per km^2).
  Higher = more animals.

- group:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column identifying the site or hunting treatment. Must contain the
  `reference` level and at least one hunted level.

- reference:

  Character scalar naming the level of `group` treated as the unhunted /
  lightly hunted reference (e.g. `"reference"`).

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

An offtake tibble with one row per hunted site (and per group, if `data`
is grouped) and the columns:

- hunted_site:

  The hunted level of `group` being assessed.

- reference_site:

  The reference level it is compared against.

- n_hunted, n_reference:

  Number of non-missing values at each site.

- hunted_value:

  Mean density at the hunted site.

- reference_value:

  Mean density at the reference site.

- pct_change:

  Percentage change of the hunted site relative to the reference,
  `(hunted - reference) / reference * 100`. Negative means the hunted
  site is depleted.

- lnrr:

  Log response ratio, `log(hunted / reference)`. A value of `-0.69`
  means the hunted site has half the reference value.

- lnrr_lo, lnrr_hi:

  Confidence interval of `lnrr` (level `1 - alpha`), from the
  delta-method variance of Hedges et al. (1999). `NA` when a site has
  fewer than two replicates.

- p_value:

  P-value of the one-sided Welch *t*-test that the hunted site has a
  *lower* mean than the reference. `NA` when no test is possible.

- p_adj:

  `p_value` adjusted across the hunted sites with `p_adjust`.

- outcome:

  `"unsustainable"`, `"no substantial decline"` or `"inconclusive"` (see
  the Outcome section).

- sustainable:

  `FALSE`, `TRUE` or `NA`, matching `outcome`.

## Details

Provide one row per density estimate (e.g. per line transect or
camera-trap station). For each hunted site the function reports the size
of the difference (percent change and log response ratio with its
confidence interval) and a one-sided Welch *t*-test (hunted \<
reference). When several hunted sites are compared with the same
reference, p-values are adjusted for multiple comparisons (`p_adjust`).

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

[`ot_hyco()`](https://stangandaho.github.io/offtake/reference/ot_hyco.md),
[`ot_asc()`](https://stangandaho.github.io/offtake/reference/ot_asc.md)

## Examples

``` r
# One row per transect. `site` = treatment, `dens` = animals per km^2.
d <- data.frame(
  site = rep(c("control", "hunted"), each = 4), # reference vs hunted
  dens = c(12, 14, 11, 13, 6, 7, 5, 8) # density per transect
)
ot_pdc(d, density = dens, group = site, reference = "control")
#> index-based: PDC (population density comparison)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference hunted_value reference_value
#> * <chr>       <chr>             <int>       <int>        <dbl>           <dbl>
#> 1 hunted      control               4           4          6.5            12.5
#> # ℹ 8 more variables: pct_change <dbl>, lnrr <dbl>, lnrr_lo <dbl>,
#> #   lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>, sustainable <lgl>

# The bundled example data
ot_pdc(bushmeat_sites, density = density, group = site_type,
       reference = "reference")
#> index-based: PDC (population density comparison)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference hunted_value reference_value
#> * <chr>       <chr>             <int>       <int>        <dbl>           <dbl>
#> 1 hunted      reference             6           6         13.5            29.5
#> # ℹ 8 more variables: pct_change <dbl>, lnrr <dbl>, lnrr_lo <dbl>,
#> #   lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>, sustainable <lgl>
```
