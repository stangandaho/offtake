# Illustrative site-level bushmeat survey data

Synthetic paired-transect data contrasting a lightly hunted *reference*
site with a *hunted* site, for use with the index functions
[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md)
and
[`ot_hyco()`](https://stangandaho.github.io/offtake/reference/ot_hyco.md).
**These are simulated values, not measurements from any real study.**

## Usage

``` r
bushmeat_sites
```

## Format

A data frame with 12 rows and 5 variables:

- site_type:

  `"reference"` or `"hunted"`.

- transect:

  Transect identifier.

- density:

  Wildlife density (individuals per km^2).

- yield_kg:

  Harvested biomass (kg).

- effort_days:

  Hunting effort (hunter-days).

## Source

Simulated.

## See also

[`ot_pdc()`](https://stangandaho.github.io/offtake/reference/ot_pdc.md),
[`ot_hyco()`](https://stangandaho.github.io/offtake/reference/ot_hyco.md)
