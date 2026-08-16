# Population density comparison (PDC)

Compares wildlife abundance/density between hunted and reference
(unhunted or lightly hunted) sites. Lower density at the hunted site(s)
is interpreted as unsustainable (Adounke et al. 2026; Weinbaum et al.
2013).

## Usage

``` r
ot_pdc(data, density, group, reference, alpha = 0.05)
```

## Arguments

- data:

  A data frame in *long* format with one row per density estimate. It
  must contain at least: a numeric density column (`density`) and a
  site/treatment column (`group`) whose values include the reference
  level and one or more hunted levels. Extra columns are ignored.

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

  Significance level for the one-sided test (default `0.05`).

## Value

An offtake tibble with **one row per hunted site** and the columns:

- hunted_site:

  The hunted level of `group` being assessed.

- reference_site:

  The reference level it is compared against.

- hunted_value:

  Mean density at the hunted site.

- reference_value:

  Mean density at the reference site.

- pct_change:

  Percentage change of the hunted site relative to the reference,
  `(hunted - reference) / reference * 100`. Negative means the hunted
  site is depleted.

- p_value:

  P-value of the one-sided Welch *t*-test that the hunted site has a
  *lower* mean than the reference. `NA` when a site has fewer than two
  replicates.

- sustainable:

  Logical verdict. `FALSE` when the hunted site is significantly lower
  than the reference (`p_value < alpha`); otherwise `TRUE`.

## Details

Provide **one row per density estimate** (e.g. per line transect or
camera-trap station). With two or more replicates per site a one-sided
Welch *t*-test is used (hunted \< reference); with a single value per
site the verdict falls back to the sign of the difference and a warning
is issued.

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
#> <offtake: PDC (population density comparison)>  (index-based)
#> Reference: Adounke et al. (2026); Weinbaum et al. (2013) 
#> 
#> # A tibble: 1 × 7
#>   hunted_site reference_site hunted_value reference_value pct_change  p_value
#> * <chr>       <chr>                 <dbl>           <dbl>      <dbl>    <dbl>
#> 1 hunted      control                 6.5            12.5        -48 0.000297
#> # ℹ 1 more variable: sustainable <lgl>
```
