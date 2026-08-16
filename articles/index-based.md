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
three return a small table with a `sustainable` column (`TRUE` or
`FALSE`).

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
#> <offtake: PDC (population density comparison)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> 
#> # A tibble: 1 × 7
#>   hunted_site reference_site hunted_value reference_value pct_change     p_value
#> * <chr>       <chr>                 <dbl>           <dbl>      <dbl>       <dbl>
#> 1 hunted      reference              13.5            29.5      -54.2 0.000000150
#> # ℹ 1 more variable: sustainable <lgl>
```

How to read the output:

- `hunted_value` and `reference_value` are the mean densities at each
  site.

- `pct_change` is the change of the hunted site relative to the
  reference,

  \\\text{pct\\change} = \frac{\text{hunted} -
  \text{reference}}{\text{reference}} \times 100.\\

  A negative value means the hunted site holds fewer animals.

- `p_value` comes from a one-sided *t*-test asking whether the hunted
  site is really lower (it is `NA` if a site has fewer than two
  replicates).

- `sustainable` is `FALSE` when the hunted site is significantly lower,
  and `TRUE` otherwise.

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
#> <offtake: HYCo (hunting yield comparison, CPUE)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> 
#> # A tibble: 1 × 7
#>   hunted_site reference_site hunted_value reference_value pct_change   p_value
#> * <chr>       <chr>                 <dbl>           <dbl>      <dbl>     <dbl>
#> 1 hunted      reference              2.03            4.20      -51.7 0.0000219
#> # ℹ 1 more variable: sustainable <lgl>
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
#> <offtake: ASC (age structure comparison, juvenile proportion)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> 
#> # A tibble: 1 × 7
#>   hunted_site reference_site juv_prop_hunted juv_prop_reference pct_change
#> * <chr>       <chr>                    <dbl>              <dbl>      <dbl>
#> 1 hunted      reference                0.158              0.375      -57.8
#> # ℹ 2 more variables: p_value <dbl>, sustainable <lgl>
```

Here `juv_prop_hunted` and `juv_prop_reference` are the juvenile shares
(between 0 and 1) at each site, and `pct_change` is the change relative
to the reference.

If you prefer to compare the whole age structure at once rather than
just the juvenile share, use `method = "distribution"`. That runs a
chi-squared test on the full class by site table and flags a significant
difference. In that case read the result together with the direction of
the shift, since a difference on its own does not say which way things
are going.

``` r

ot_asc(bushmeat_ages,
       stage = age_class,
       group = site_type,
       reference = "reference",
       method = "distribution")
#> <offtake: ASC (age structure comparison, distribution)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> Note: A significant class x site difference is flagged as unsustainable; interpret alongside the direction of the shift (Weinbaum et al. 2013). 
#> 
#> # A tibble: 1 × 6
#>   hunted_site reference_site statistic    df  p_value sustainable
#> * <chr>       <chr>              <dbl> <int>    <dbl> <lgl>      
#> 1 hunted      reference           13.3     1 0.000263 FALSE
```

## A note of caution

These indices are quick screens. A `FALSE` verdict says the hunted site
looks worse on this one measure, which is worth following up. It does
not prove that hunting is the cause. Whenever you can, make sure the two
sites are otherwise comparable, and read the three indices together
rather than trusting any single one.

## References

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.

Weinbaum, K. Z., Brashares, J. S., Golden, C. D. & Getz, W. M. (2013)
Searching for sustainability: are assessments of wildlife harvests
behind the times? *Ecology Letters* 16, 99 to 111.
