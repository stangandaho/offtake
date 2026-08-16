# Biodemographic spatial source-sink model (Biode)

Predicts the steady-state spatial distribution of game density across a
landscape given the location and human population of hunting
settlements, following the spatially explicit reaction-diffusion
source-sink model of Levi et al. (2009). Unlike the non-spatial models
([`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md),
[`ot_pbr()`](https://stangandaho.github.io/offtake/reference/ot_pbr.md),
[`ot_msy()`](https://stangandaho.github.io/offtake/reference/ot_msy.md)),
Biode maps where hunting creates population *sinks* around settlements
and *sources* in the remote matrix.

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
  resolution = 1,
  margin = 3 * sigma,
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

  Carrying-capacity density `K` (individuals per cell area).

- r:

  Intrinsic population growth rate.

- hphy:

  Hunts per hunter per year (`h`).

- kill_rate:

  Per-encounter kill probability `d` (weapon efficiency, 0-1).

- sigma:

  Spatial spread of hunting effort \\\sigma\\ (mean hunting distance
  from a settlement), in the coordinate unit.

- er:

  Coefficient converting game density to encounter rate `e` (default
  `0.02`, the value used by Levi et al. 2009).

- theta:

  Density-dependence exponent \\\theta\\ (default `1`, logistic).

- resolution:

  Grid cell size in coordinate units (default `1`).

- margin:

  Landscape margin added around the settlement bounding box, in
  coordinate units (default `3 * sigma`).

- extirpation_threshold:

  Density at or below which a cell is counted as locally extirpated
  (default `1`).

- max_extirpated:

  Maximum fraction of the landscape allowed to be extirpated for the
  assessment to be flagged `sustainable` (default `0.05`).

## Value

A one-row offtake tibble summarising the landscape, with the columns:

- n_settlements:

  Number of settlements in `data`.

- mean_density:

  Mean predicted game density across all grid cells.

- min_density:

  Lowest predicted density (0 where a sink is fully extirpated).

- frac_extirpated:

  Fraction of grid cells at or below `extirpation_threshold` – the share
  of the landscape locally wiped out.

- frac_below_half_k:

  Fraction of cells with density below half of carrying capacity `k` – a
  broader depletion footprint.

- cpue:

  Effort-weighted mean catch per unit effort predicted across the
  landscape.

- sustainable:

  Logical verdict: `TRUE` when `frac_extirpated <= max_extirpated`.

The full predicted density surface and per-settlement local CPUE are
attached as attributes and are most easily retrieved with
[`ot_biode_surface()`](https://stangandaho.github.io/offtake/reference/ot_biode_surface.md)
and
[`ot_biode_cpue()`](https://stangandaho.github.io/offtake/reference/ot_biode_cpue.md).

## Details

At steady state the model reduces to a closed form for the game density
at each point of the landscape: \$\$N(x,y) = \Big\[\\K^{\theta}\Big(1 -
\frac{e\\ d\\ h}{r}\\ \sum\_{c} P_c\\
\varphi(\delta_c)\Big)\Big\]\_{+}^{1/\theta},\$\$ with the Gaussian
hunting-effort kernel \$\$\varphi(\delta) = \frac{\exp\\\big(-\delta^2 /
2\sigma^2\big)} {2\pi\\\delta + 1},\$\$ where \\\delta_c\\ is the
distance from a point to settlement `c`, `P_c` its hunter population,
`K` carrying capacity, `r` the intrinsic growth rate, `h` (`hphy`) hunts
per hunter per year, `d` (`kill_rate`) the per-encounter kill
probability (weapon efficiency), `e` (`er`) the density-to-encounter
coefficient, \\\theta\\ the density-dependence exponent, and
\\\[\cdot\]\_+\\ denotes truncation at zero. The form is taken directly
from Levi et al. (2009) and their published solver code.

## References

Levi, T., Shepard, G. H., Ohl-Schacherer, J., Peres, C. A. & Yu, D. W.
(2009) Modelling the long-term sustainability of indigenous hunting in
Manu National Park, Peru: landscape-scale management implications for
Amazonia. *Journal of Applied Ecology* 46, 804-814.

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
res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
                k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
                resolution = 2)
res
#> <offtake: Biode (Levi et al. spatial source-sink model)>  (model-based)
#> Reference: Levi et al. (2009) 
#> 
#> # A tibble: 1 × 7
#>   n_settlements mean_density min_density frac_extirpated frac_below_half_k
#> *         <int>        <dbl>       <dbl>           <dbl>             <dbl>
#> 1             2         17.8           0           0.163             0.258
#> # ℹ 2 more variables: cpue <dbl>, sustainable <lgl>
head(ot_biode_surface(res))
#> # A tibble: 6 × 3
#>       x     y density
#>   <dbl> <dbl>   <dbl>
#> 1    -8    -3    25.0
#> 2    -8    -1    25.0
#> 3    -8     1    25.0
#> 4    -8     3    25.0
#> 5    -8     5    25.0
#> 6    -8     7    24.9
```
