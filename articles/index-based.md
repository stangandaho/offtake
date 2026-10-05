# Index-based approaches

## The idea

Index-based methods ask a simple question: is the hunted site worse off
than a site that is not hunted (or only lightly hunted)? You measure the
same thing at both places, for example animal density or hunting yield,
and compare them. If the hunted site is clearly lower, that is a warning
sign that hunting may not be sustainable.

These methods are popular because they need little data and are easy to
explain (Adounke et al. 2026; Weinbaum et al. 2013). Their main limit is
that they only compare sites. They do not tell you the exact level of
hunting at which things tip over, and a difference between two sites can
have causes other than hunting (soil, habitat, or the size of the forest
patch). So treat the verdict as a signal to look closer, not as a final
answer.

The package covers three of these indices:

| Function | Compares | Warning sign |
|----|----|----|
| [`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md) | animal density | lower density at the hunted site |
| [`ot_hyco()`](https://stangandaho.github.io/offtake/reference/ot_hyco.md) | hunting yield (or catch per unit effort) | lower yield at the more hunted site |
| [`ot_asc()`](https://stangandaho.github.io/offtake/reference/ot_asc.md) | share of young animals | fewer juveniles at the hunted site |

All three take a tidy data frame and unquoted column names, and all
three return a small table with an `outcome` and a `sustainable` column.

We use the bundled example data `bushmeat_sites` (density and yield per
transect) and `bushmeat_ages` (one row per animal). Both are made up for
the sake of the examples, not real field data.

``` r

head(bushmeat_sites)
#>   site_type transect density yield_kg effort_days
#> 1 reference      T01      28       41          10
#> 2 reference      T02      31       44          11
#> 3 reference      T03      26       39           9
#> 4 reference      T04      33       45          10
#> 5 reference      T05      29       42          12
#> 6 reference      T06      30       43           9
```

## Population density comparison, `ot_pdc()`

Here you compare animal density at the hunted site with the reference
site. You give one row per density estimate (for example one row per
transect), the density column, the column that marks the site, and which
site is the reference.

``` r

ot_pdc(bushmeat_sites,
       density = density,
       group = site_type,
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

The output answers two questions.

**How big is the gap?** `hunted_value` and `reference_value` are the
mean densities at each site, and `pct_change` is the change of the
hunted site relative to the reference,

\\\text{pct\\change} = \frac{\bar{x}\_{\text{hunted}} -
\bar{x}\_{\text{reference}}}{\bar{x}\_{\text{reference}}} \times 100.\\

The same gap on the log scale is the log response ratio, a standard
effect size in ecology (Hedges et al. 1999):

\\\text{lnrr} =
\ln\frac{\bar{x}\_{\text{hunted}}}{\bar{x}\_{\text{reference}}}, \qquad
\text{SE} = \sqrt{\frac{s_h^2}{n_h\\\bar{x}\_h^2} +
\frac{s_r^2}{n_r\\\bar{x}\_r^2}}.\\

A value of \\-0.69\\ means the hunted site has half the reference value.
`lnrr_lo` and `lnrr_hi` give its 95% confidence interval.

**Is it more than chance?** `p_value` comes from a one-sided Welch
*t*-test asking whether the hunted site is really lower. When several
hunted sites are compared with the same reference, `p_adj` corrects for
the number of comparisons (Holm’s method by default).

## Reading the outcome

The `outcome` column has three levels, so that “no difference found” is
never mistaken for “sustainable”:

| `outcome` | `sustainable` | Meaning |
|----|----|----|
| `"unsustainable"` | `FALSE` | The hunted site is significantly lower. |
| `"no substantial decline"` | `TRUE` | Not significant, and the confidence interval rules out a decline of 20% or more (`decline` argument). |
| `"inconclusive"` | `NA` | Too few replicates, or an interval so wide that a large decline cannot be ruled out. |

Here is what two noisy transects per site give. The hunted site looks
17% lower, but the data cannot tell a small decline from a large one:

``` r

small <- data.frame(site = c("ref", "ref", "hunted", "hunted"),
                    density = c(10, 14, 11, 9))
ot_pdc(small, density = density, group = site, reference = "ref")
#> index-based: PDC (population density comparison)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference hunted_value reference_value
#> * <chr>       <chr>             <int>       <int>        <dbl>           <dbl>
#> 1 hunted      ref                   2           2           10              12
#> # ℹ 8 more variables: pct_change <dbl>, lnrr <dbl>, lnrr_lo <dbl>,
#> #   lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>, sustainable <lgl>
```

## Hunting yield comparison, `ot_hyco()`

The logic is the same, but the metric is what hunters bring back rather
than what is left in the forest. If you also pass an `effort` column,
the function compares catch per unit effort (yield divided by effort)
instead of raw yield. That is often more honest, because more effort
naturally brings more catch.

``` r

ot_hyco(bushmeat_sites,
        yield = yield_kg,
        group = site_type,
        reference = "reference",
        effort = effort_days)
#> index-based: HYCo (hunting yield comparison, CPUE)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference hunted_value reference_value
#> * <chr>       <chr>             <int>       <int>        <dbl>           <dbl>
#> 1 hunted      reference             6           6         2.03            4.20
#> # ℹ 8 more variables: pct_change <dbl>, lnrr <dbl>, lnrr_lo <dbl>,
#> #   lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>, sustainable <lgl>
```

The output columns match
[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md).
Here `hunted_value` and `reference_value` hold the mean catch per
hunter-day at each site.

## Age structure comparison, `ot_asc()`

This one looks at who is being caught. In a healthy population you
expect a fair share of young animals. When hunting is heavy, juveniles
often become scarce, which is a sign that the population is not
replacing itself. So a lower share of juveniles at the hunted site is
the warning sign.

You give one row per animal, the class column, the site column, and
which class counts as juvenile.

``` r

ot_asc(bushmeat_ages,
       stage = age_class,
       group = site_type,
       reference = "reference",
       juvenile = "juvenile")
#> index-based: ASC (age structure comparison, juvenile proportion)
#> 
#> # A tibble: 1 × 14
#>   hunted_site reference_site n_hunted n_reference juv_prop_hunted
#> * <chr>       <chr>             <dbl>       <dbl>           <dbl>
#> 1 hunted      reference           120         120           0.158
#> # ℹ 9 more variables: juv_prop_reference <dbl>, pct_change <dbl>, lnrr <dbl>,
#> #   lnrr_lo <dbl>, lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>,
#> #   sustainable <lgl>
```

Here `juv_prop_hunted` and `juv_prop_reference` are the juvenile shares
(between 0 and 1) at each site, and `lnrr` is the log ratio of the two
shares. The test is a one-sided test of two proportions. If your data
are already tallied, pass the column of counts to `count`.

If you prefer to compare the whole age structure at once, use
`method = "distribution"`. That runs a chi-squared test on the full
class by site table. A chi-squared test cannot say in which direction
the structure shifted, nor show that two structures are alike, so a
significant difference gives `"unsustainable"` and anything else gives
`"inconclusive"`.

``` r

ot_asc(bushmeat_ages,
       stage = age_class,
       group = site_type,
       reference = "reference",
       method = "distribution")
#> index-based: ASC (age structure comparison, distribution)
#> Note: A chi-squared test cannot give the direction of the shift nor show similarity: a significant difference is flagged, anything else is inconclusive. 
#> 
#> # A tibble: 1 × 8
#>   hunted_site reference_site statistic    df p_value   p_adj outcome sustainable
#> * <chr>       <chr>              <dbl> <int>   <dbl>   <dbl> <chr>   <lgl>      
#> 1 hunted      reference           13.3     1 2.63e-4 2.63e-4 unsust… FALSE
```

## Several species at once

All three functions accept data grouped with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
Each group (for example each species) is then analysed on its own, and
the group column is kept in the result.

``` r

two_species <- rbind(
  transform(bushmeat_sites, species = "blue duiker"),
  transform(bushmeat_sites, species = "red duiker",
            density = density * rep(c(1, 1.8), each = 6))
)
ot_pdc(dplyr::group_by(two_species, species),
       density = density,
       group = site_type,
       reference = "reference")
#> index-based: PDC (population density comparison)
#> 
#> # A tibble: 2 × 15
#>   species     hunted_site reference_site n_hunted n_reference hunted_value
#> * <chr>       <chr>       <chr>             <int>       <int>        <dbl>
#> 1 blue duiker hunted      reference             6           6         13.5
#> 2 red duiker  hunted      reference             6           6         24.3
#> # ℹ 9 more variables: reference_value <dbl>, pct_change <dbl>, lnrr <dbl>,
#> #   lnrr_lo <dbl>, lnrr_hi <dbl>, p_value <dbl>, p_adj <dbl>, outcome <chr>,
#> #   sustainable <lgl>
```

## A note of caution

These indices are quick screens. A `FALSE` verdict says the hunted site
looks worse on this one measure, which is worth following up. It does
not prove that hunting is the cause. Two points deserve care.

- **Comparable sites.** The two sites should differ in hunting only, as
  far as possible (habitat, soil, size of the forest patch).
- **Pseudoreplication.** Transects inside one hunted site and one
  reference site are not independent replicates of hunting. The test
  then shows that the two *sites* differ, not that *hunting* causes the
  difference (Hurlbert 1984). Several hunted and several reference sites
  make a much stronger design.

Whenever you can, read the three indices together rather than trusting
any single one.

## References

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.

Hedges, L. V., Gurevitch, J. & Curtis, P. S. (1999) The meta-analysis of
response ratios in experimental ecology. *Ecology* 80, 1150 to 1156.

Hurlbert, S. H. (1984) Pseudoreplication and the design of ecological
field experiments. *Ecological Monographs* 54, 187 to 211.

Weinbaum, K. Z., Brashares, J. S., Golden, C. D. & Getz, W. M. (2013)
Searching for sustainability: are assessments of wildlife harvests
behind the times? *Ecology Letters* 16, 99 to 111.
