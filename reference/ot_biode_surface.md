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

  Predicted steady-state game density in that cell.

This long format is ready for mapping, e.g. with
`ggplot2::geom_raster(ggplot2::aes(x, y, fill = density))`.

## See also

[`ot_biode()`](https://stangandaho.github.io/offtake/reference/ot_biode.md)

## Examples

``` r
settlements <- data.frame(xkm = 10, ykm = 10, hunters = 100)
res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
                k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
                resolution = 2)
ot_biode_surface(res)
#> # A tibble: 361 × 3
#>        x     y density
#>    <dbl> <dbl>   <dbl>
#>  1    -8    -8    25.0
#>  2    -8    -6    25.0
#>  3    -8    -4    25.0
#>  4    -8    -2    25.0
#>  5    -8     0    24.9
#>  6    -8     2    24.9
#>  7    -8     4    24.8
#>  8    -8     6    24.8
#>  9    -8     8    24.7
#> 10    -8    10    24.7
#> # ℹ 351 more rows
```
