# Extract per-settlement catch per unit effort from a Biode result

Extract per-settlement catch per unit effort from a Biode result

## Usage

``` r
ot_biode_cpue(x)
```

## Arguments

- x:

  An
  [offtake](https://stangandaho.github.io/offtake/reference/ot_biode.md)
  object returned by
  [`ot_biode()`](https://stangandaho.github.io/offtake/reference/ot_biode.md).

## Value

A tibble with **one row per settlement** and the columns:

- settlement:

  Row index of the settlement in the input data.

- x, y:

  Settlement coordinates.

- humans:

  Number of hunters at the settlement.

- local_cpue:

  Expected kills per hunt for that settlement's hunters: the density
  averaged over where they hunt (weighted by their own effort kernel)
  times `er * kill_rate`. Lower values indicate a more depleted hunting
  ground. See the `cpue` column of
  [`ot_biode()`](https://stangandaho.github.io/offtake/reference/ot_biode.md)
  for how this relates to the catch per unit effort of Levi et al.
  (2011).

## See also

[`ot_biode()`](https://stangandaho.github.io/offtake/reference/ot_biode.md)

## Examples

``` r
settlements <- data.frame(xkm = c(10, 30), ykm = c(10, 20), hunters = c(80, 120))
res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
                k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
                extent = c(-10, 50, -10, 40), resolution = 2)
ot_biode_cpue(res)
#> # A tibble: 2 × 5
#>   settlement     x     y humans local_cpue
#>        <int> <dbl> <dbl>  <dbl>      <dbl>
#> 1          1    10    10     80    0.00636
#> 2          2    30    20    120    0.00461
```
