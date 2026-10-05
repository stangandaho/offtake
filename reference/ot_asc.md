# Age structure comparison (ASC)

Compares the age/sex structure of a population between hunted and
reference (unhunted or lightly hunted) sites. Following Adounke et al.
(2026), a **lower proportion of the juvenile class at the hunted site is
interpreted as unsustainable** (a signal of recruitment failure).
Optionally the full class-by-site frequency distribution can be compared
with a chi-squared test.

## Usage

``` r
ot_asc(
  data,
  stage,
  group,
  reference,
  juvenile = NULL,
  count = NULL,
  method = c("proportion", "distribution"),
  alpha = 0.05,
  p_adjust = "holm",
  decline = 0.2
)
```

## Arguments

- data:

  A data frame with one row per sampled/harvested individual (or one row
  per class if `count` is supplied). It must contain an age/sex class
  column (`stage`) and a site column (`group`). If your data are already
  tallied, add a `count` column and pass it to `count`. May be grouped
  with
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).

- stage:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column giving the age or sex class of each individual (e.g.
  `"juvenile"`/`"adult"`).

- group:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column identifying the site / hunting treatment.

- reference:

  Character scalar naming the level of `group` treated as the unhunted /
  lightly hunted reference.

- juvenile:

  Character vector of the level(s) of `stage` that represent the
  juvenile (pre-reproductive) class. Required for
  `method = "proportion"`.

- count:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Optional column of counts, used when `data` is already aggregated to
  class \\\times\\ site.

- method:

  `"proportion"` (default) compares the juvenile proportion between the
  two sites; `"distribution"` compares the whole class \\\times\\ site
  table with a chi-squared test.

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
is grouped). Columns depend on `method`.

For `method = "proportion"`:

- hunted_site, reference_site:

  The two levels being compared.

- n_hunted, n_reference:

  Number of individuals at each site.

- juv_prop_hunted:

  Proportion of individuals in the juvenile class at the hunted site
  (0-1).

- juv_prop_reference:

  Juvenile proportion at the reference site.

- pct_change:

  Percentage change in juvenile proportion relative to the reference;
  negative means fewer juveniles at the hunted site.

- lnrr, lnrr_lo, lnrr_hi:

  Log ratio of the two juvenile proportions and its confidence interval
  (level `1 - alpha`).

- p_value:

  One-sided
  [`stats::prop.test()`](https://rdrr.io/r/stats/prop.test.html) p-value
  for a *lower* juvenile proportion at the hunted site.

- p_adj:

  `p_value` adjusted across the hunted sites with `p_adjust`.

- outcome, sustainable:

  Three-level verdict, as in
  [`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md).

For `method = "distribution"`:

- hunted_site, reference_site:

  The two levels being compared.

- statistic:

  Pearson chi-squared statistic for the class \\\times\\ site table.

- df:

  Degrees of freedom of the test.

- p_value, p_adj:

  Chi-squared p-value, raw and adjusted.

- outcome:

  `"unsustainable"` when the two structures differ significantly,
  otherwise `"inconclusive"`.

- sustainable:

  `FALSE` or `NA`, matching `outcome`.

## Details

With `method = "proportion"` the size of the difference is reported as a
log ratio of the two juvenile proportions with its confidence interval,
the test is a one-sided two-sample test of proportions, and the
three-level outcome is that of
[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md).
With `method = "distribution"` a chi-squared test cannot tell the
direction of the shift, nor show that two structures are similar: a
significant difference gives `"unsustainable"` (as in Weinbaum et al.
2013) and anything else `"inconclusive"`.

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
[`ot_hyco()`](https://stangandaho.github.io/offtake/reference/ot_hyco.md)

## Examples

``` r
# Individual-level data: one row per animal, `class` = age class of that
# animal, `site` = where it was sampled.
set.seed(1)
d <- data.frame(
  site = rep(c("control", "hunted"), c(60, 60)),
  class = c(sample(c("juv", "adult"), 60, TRUE, c(0.40, 0.60)),
            sample(c("juv", "adult"), 60, TRUE, c(0.15, 0.85)))
)
ot_asc(d, stage = class, group = site, reference = "control", juvenile = "juv")
#> index-based: ASC (age structure comparison, juvenile proportion)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference juv_prop_hunted
#> * <chr>       <chr>             <dbl>       <dbl>           <dbl>
#> 1 hunted      control              60          60            0.15
#> # ℹ 9 more variables: juv_prop_reference <dbl>, pct_change <dbl>, lnrr <dbl>,
#> #   lnrr_lo <dbl>, lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>,
#> #   sustainable <lgl>
```
