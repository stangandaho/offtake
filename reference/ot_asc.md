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
  alpha = 0.05
)
```

## Arguments

- data:

  A data frame with **one row per sampled/harvested individual** (or one
  row per class if `count` is supplied). It must contain an age/sex
  class column (`stage`) and a site column (`group`). If your data are
  already tallied, add a `count` column and pass it to `count`.

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
  two sites with a two-sample test of proportions; `"distribution"`
  compares the whole class \\\times\\ site table with a chi-squared
  test.

- alpha:

  Significance level (default `0.05`).

## Value

An offtake tibble with **one row per hunted site**. Columns depend on
`method`.

For `method = "proportion"`:

- hunted_site, reference_site:

  The two levels being compared.

- juv_prop_hunted:

  Proportion of individuals in the juvenile class at the hunted site
  (0-1).

- juv_prop_reference:

  Juvenile proportion at the reference site.

- pct_change:

  Percentage change in juvenile proportion relative to the reference;
  negative means fewer juveniles at the hunted site.

- p_value:

  One-sided
  [`stats::prop.test()`](https://rdrr.io/r/stats/prop.test.html) p-value
  for a *lower* juvenile proportion at the hunted site.

- sustainable:

  `FALSE` when the juvenile share is significantly lower at the hunted
  site.

For `method = "distribution"`:

- hunted_site, reference_site:

  The two levels being compared.

- statistic:

  Pearson chi-squared statistic for the class \\\times\\ site table.

- df:

  Degrees of freedom of the test.

- p_value:

  Chi-squared p-value.

- sustainable:

  `FALSE` when the two age structures differ significantly
  (`p_value < alpha`); interpret alongside the direction of the shift
  (Weinbaum et al. 2013).

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
#> <offtake: ASC (age structure comparison, juvenile proportion)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> 
#> # A tibble: 1 × 7
#>   hunted_site reference_site juv_prop_hunted juv_prop_reference pct_change
#> * <chr>       <chr>                    <dbl>              <dbl>      <dbl>
#> 1 hunted      control                   0.15              0.417        -64
#> # ℹ 2 more variables: p_value <dbl>, sustainable <lgl>
```
