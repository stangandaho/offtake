# Illustrative duiker demographic parameters

Synthetic species-level life-history and offtake parameters for the
non-spatial harvest models
[`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md),
[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md),
[`ot_msy()`](https://stangandaho.github.io/offtake/reference/ot_msy.md)
and
[`ot_samse()`](https://stangandaho.github.io/offtake/reference/ot_samse.md).
Ranges are loosely inspired by the duiker literature but the values are
**simulated and illustrative only.**

## Usage

``` r
duiker_demography
```

## Format

A data frame with 3 rows and 8 variables:

- species:

  Common name.

- density_k:

  Carrying-capacity density (individuals per km^2).

- annual_take:

  Observed annual offtake (individuals per km^2 per year).

- lifespan:

  Age at last reproduction (years).

- b:

  Mean annual female offspring per female.

- a:

  Age at first reproduction (years).

- w:

  Age at last reproduction (years).

- rmax:

  Maximum instantaneous rate of increase (per year).

## Source

Simulated.

## See also

[`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md),
[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md),
[`ot_msy()`](https://stangandaho.github.io/offtake/reference/ot_msy.md),
[`ot_samse()`](https://stangandaho.github.io/offtake/reference/ot_samse.md)
