# Extract the predicted density surface from a Biode result

Extract the predicted density surface from a Biode result

## Usage

``` r
ot_biode_surface(x)
```

## Arguments

- x:

  An
  [offtake](https://stangandaho.github.io/offtake/reference/ot_biode.md)
  object returned by
  [`ot_biode()`](https://stangandaho.github.io/offtake/reference/ot_biode.md).

## Value

A tibble with **one row per grid cell** and the columns:

- x:

  Cell-centre x coordinate (coordinate unit of the settlements).

- y:

  Cell-centre y coordinate.

- density:

  Predicted game density in that cell; `NA` outside the study area when
  a `boundary` was given.

This long format is ready for mapping, e.g. with
`ggplot2::geom_raster(ggplot2::aes(x, y, fill = density))`.

## See also

[`ot_biode()`](https://stangandaho.github.io/offtake/reference/ot_biode.md)

## Examples

``` r
settlements <- data.frame(xkm = 10, ykm = 10, hunters = 100)
res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
                k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
                extent = c(-10, 30, -10, 30), resolution = 2)
ot_biode_surface(res)
#> # A tibble: 400 × 3
#>        x     y density
#>    <dbl> <dbl>   <dbl>
#>  1    -9    -9    25.0
#>  2    -9    -7    25.0
#>  3    -9    -5    25.0
#>  4    -9    -3    25.0
#>  5    -9    -1    25.0
#>  6    -9     1    25.0
#>  7    -9     3    24.9
#>  8    -9     5    24.9
#>  9    -9     7    24.9
#> 10    -9     9    24.8
#> # ℹ 390 more rows
```
