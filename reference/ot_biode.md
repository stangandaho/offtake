# Biodemographic spatial model of hunting depletion (Biode)

Predicts the spatial distribution of game density across a landscape
given the location and number of hunters of each settlement, following
the spatially explicit model of Levi et al. (2009), in its
multi-settlement form (Levi et al. 2011). Unlike the non-spatial models
([`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md),
[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md),
[`ot_msy()`](https://stangandaho.github.io/offtake/reference/ot_msy.md)),
Biode shows *where* hunting depletes game: around settlements, while
remote areas stay close to carrying capacity.

## Usage

``` r
ot_biode(
  data,
  x,
  y,
  humans,
  k,
  r,
  hphy,
  kill_rate,
  sigma,
  er = 0.02,
  theta = 1,
  method = c("steady", "dynamic"),
  extent = NULL,
  boundary = NULL,
  resolution = 1,
  margin = 3 * sigma,
  years = 100,
  dispersal = 0.02,
  extirpation_threshold = 1,
  max_extirpated = 0.05
)
```

## Arguments

- data:

  A data frame of settlements, **one row per settlement**. Required
  columns: the settlement coordinates (`x`, `y`) in the same length unit
  as `sigma` (e.g. kilometres), and the number of hunters (`humans`).
  The biological and hunting parameters (`k`, `r`, `hphy`, `kill_rate`,
  `sigma`, ...) are passed as scalar arguments, not columns, because
  they describe the target species and hunting technology, not
  individual settlements.

- x, y:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Settlement coordinate columns, in the same length unit as `sigma`.

- humans:

  \<[`data-masked`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Column giving the number of hunters (human population) at each
  settlement.

- k:

  Carrying-capacity density `K` (individuals per unit area).

- r:

  Intrinsic population growth rate (per year).

- hphy:

  Hunts per hunter per year (`h`).

- kill_rate:

  Per-encounter kill probability `d` (weapon efficiency, 0-1).

- sigma:

  Spatial spread of hunting effort \\\sigma\\ (mean hunting distance
  from a settlement), in the coordinate unit.

- er:

  Coefficient converting game density to encounter rate `e` (default
  `0.02`, as in Levi et al. 2011).

- theta:

  Density-dependence exponent \\\theta\\ (default `1`, logistic).

- method:

  `"steady"` (closed-form steady state, no movement) or `"dynamic"`
  (time-stepped, with movement between cells). See the Two methods
  section.

- extent:

  Optional study area as a rectangle, `c(xmin, xmax, ymin, ymax)` in the
  coordinate unit.

- boundary:

  Optional study area as a polygon: a data frame or matrix whose first
  two columns are the x and y coordinates of its vertices, in order.
  Takes precedence over `extent`. (From an `sf` polygon, use
  `sf::st_coordinates(polygon)[, 1:2]`.)

- resolution:

  Grid cell size in coordinate units (default `1`).

- margin:

  Buffer, in coordinate units, by which the simulation grid extends
  beyond the study area (default `3 * sigma`, rounded up to whole
  cells). When neither `extent` nor `boundary` is given, the study area
  itself is the settlements' bounding box plus `margin`.

- years:

  Number of years simulated by `method = "dynamic"` (default `100`).

- dispersal:

  Share of the animals of a cell moving to each of its four neighbours
  every year, between 0 and 0.25 (default `0.02`, as in Levi et al.
  2011). It is defined per cell, so it depends on `resolution`. Used
  only by `method = "dynamic"`.

- extirpation_threshold:

  Density at or below which a cell is counted as locally extirpated
  (default `1`).

- max_extirpated:

  Maximum fraction of the study area allowed to be extirpated for the
  assessment to be flagged `sustainable` (default `0.05`).

## Value

A one-row offtake tibble summarising the study area, with the columns:

- n_settlements:

  Number of settlements in `data`.

- area:

  Size of the study area, in squared coordinate units (km^2 when
  coordinates are in km): the number of grid cells whose centre lies in
  the study area times the cell area, so it approximates the true area
  at the resolution of the grid.

- mean_density:

  Mean predicted game density over the study area.

- min_density:

  Lowest predicted density (0 where a sink is fully extirpated).

- area_extirpated:

  Area where the density is at or below `extirpation_threshold`, in the
  unit of `area`.

- frac_extirpated:

  `area_extirpated / area`: the share of the study area locally wiped
  out.

- frac_below_half_k:

  Share of the study area with density below half of carrying capacity
  `k`, a broader depletion footprint.

- cpue:

  Expected kills per hunt, averaged over all hunts on the grid: the
  effort-weighted mean density times `er * kill_rate`. This definition
  is specific to this package and differs in scale and weighting from
  the catch per unit effort reported by Levi et al. (2011), so compare
  settlements with each other rather than with values from that study.

- sustainable:

  Logical verdict: `TRUE` when `frac_extirpated <= max_extirpated`.

The full predicted density surface and per-settlement CPUE are attached
as attributes and are most easily retrieved with
[`ot_biode_surface()`](https://stangandaho.github.io/offtake/reference/ot_biode_surface.md)
and
[`ot_biode_cpue()`](https://stangandaho.github.io/offtake/reference/ot_biode_cpue.md).

## Hunting pressure

Each settlement `c` spreads hunting effort around it with the kernel
\$\$\varphi(\delta) = \frac{\exp\\\big(-\delta^2 / 2\sigma^2\big)}
{2\pi\\\delta + 1},\$\$ where \\\delta\\ is the distance to the
settlement and \\\sigma\\ the spread of hunting effort. The annual
per-capita kill rate in a cell is \\e\\d\\h \sum_c
P_c\\\varphi(\delta_c)\\, with `P_c` the number of hunters, `h` (`hphy`)
hunts per hunter per year, `d` (`kill_rate`) the probability that an
encounter ends in a kill and `e` (`er`) the coefficient converting game
density into encounters (Levi et al. 2009, 2011).

## Two methods

- `method = "steady"` (default) uses the closed-form steady state of the
  multi-settlement model (Levi et al. 2011). In each cell, logistic
  growth balances hunting: \$\$N(x,y) = \Big\[\\K^{\theta}\Big(1 -
  \frac{e\\ d\\ h}{r}\\ \sum\_{c} P_c\\
  \varphi(\delta_c)\Big)\Big\]\_{+}^{1/\theta}.\$\$ Animals do **not**
  move between cells in this form: it maps local depletion but not the
  resupply of hunted areas by immigration.

- `method = "dynamic"` runs the time-stepped model of Levi et al.
  (2011), which adds animal movement: each year a fraction `dispersal`
  of the animals in a cell moves to each of its four neighbours, then
  the population grows logistically and is hunted, for `years` years
  starting from carrying capacity everywhere. Remote, lightly hunted
  cells then act as *sources* that resupply the hunted *sinks* around
  settlements. Edges are reflecting (no animals leave the grid), a
  choice of this package; the implementation published with Levi et
  al. (2011) handles the grid edges differently, which only matters
  close to them.

## Study area

Areas and fractions are computed inside the study area, so the verdict
does not depend on an arbitrary map size. Give the study area as
`extent` (a rectangle) or `boundary` (a polygon). The model itself runs
on a grid that extends `margin` beyond the study area, so that the
hunting grounds of settlements near its edge, and animal movement across
it, are not cut off. Without a study area, the grid is the bounding box
of the settlements plus `margin`, a warning is issued, and
`frac_extirpated` then changes with `margin`.

## References

Levi, T., Shepard, G. H., Ohl-Schacherer, J., Peres, C. A. & Yu, D. W.
(2009) Modelling the long-term sustainability of indigenous hunting in
Manu National Park, Peru: landscape-scale management implications for
Amazonia. *Journal of Applied Ecology* 46, 804-814.

Levi, T., Shepard, G. H., Ohl-Schacherer, J., Wilmers, C. C., Peres, C.
A. & Yu, D. W. (2011) Spatial tools for modeling the sustainability of
subsistence hunting in tropical forests. *Ecological Applications* 21,
1802-1818.

Adounke, G. R. M. et al. (2026) Systematic review of sustainability
assessment approaches for wildlife exploitation. *Biological
Conservation* 313, 111606.
[doi:10.1016/j.biocon.2025.111606](https://doi.org/10.1016/j.biocon.2025.111606)

## See also

[`ot_biode_surface()`](https://stangandaho.github.io/offtake/reference/ot_biode_surface.md),
[`ot_biode_cpue()`](https://stangandaho.github.io/offtake/reference/ot_biode_cpue.md)

## Examples

``` r
# One row per settlement. `xkm`/`ykm` = location (km), `hunters` = people.
settlements <- data.frame(
  village = c("A", "B"),
  xkm = c(10, 30),
  ykm = c(15, 20),
  hunters = c(80, 120)
)
# Study area: a 50 x 40 km rectangle
res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
                k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
                extent = c(-5, 45, -5, 35), resolution = 2)
res
#> model-based: Biode (Levi et al. spatial model, steady)
#> 
#> # A tibble: 1 × 9
#>   n_settlements  area mean_density min_density area_extirpated frac_extirpated
#> *         <int> <dbl>        <dbl>       <dbl>           <dbl>           <dbl>
#> 1             2  2000         16.3           0             392           0.196
#> # ℹ 3 more variables: frac_below_half_k <dbl>, cpue <dbl>, sustainable <lgl>
head(ot_biode_surface(res))
#> # A tibble: 6 × 3
#>       x     y density
#>   <dbl> <dbl>   <dbl>
#> 1    -4    -4    25.0
#> 2    -4    -2    25.0
#> 3    -4     0    24.9
#> 4    -4     2    24.9
#> 5    -4     4    24.8
#> 6    -4     6    24.5

# With animal movement between cells
ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
         k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
         extent = c(-5, 45, -5, 35), resolution = 2, method = "dynamic")
#> model-based: Biode (Levi et al. spatial model, dynamic)
#> 
#> # A tibble: 1 × 9
#>   n_settlements  area mean_density min_density area_extirpated frac_extirpated
#> *         <int> <dbl>        <dbl>       <dbl>           <dbl>           <dbl>
#> 1             2  2000         16.6    0.000235             192           0.096
#> # ℹ 3 more variables: frac_below_half_k <dbl>, cpue <dbl>, sustainable <lgl>
```
