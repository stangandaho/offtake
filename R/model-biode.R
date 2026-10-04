#' Biodemographic spatial model of hunting depletion (Biode)
#'
#' Predicts the spatial distribution of game density across a landscape given
#' the location and number of hunters of each settlement, following the
#' spatially explicit model of Levi et al. (2009), in its multi-settlement form
#' (Levi et al. 2011). Unlike the non-spatial models
#' ([ot_pro()], [ot_pbr()], [ot_msy()]), Biode shows *where* hunting depletes
#' game: around settlements, while remote areas stay close to carrying
#' capacity.
#'
#' @section Hunting pressure:
#' Each settlement `c` spreads hunting effort around it with the kernel
#' \deqn{\varphi(\delta) = \frac{\exp\!\big(-\delta^2 / 2\sigma^2\big)}
#'        {2\pi\,\delta + 1},}
#' where \eqn{\delta} is the distance to the settlement and \eqn{\sigma} the
#' spread of hunting effort. The annual per-capita kill rate in a cell is
#' \eqn{e\,d\,h \sum_c P_c\,\varphi(\delta_c)}, with `P_c` the number of hunters,
#' `h` (`hphy`) hunts per hunter per year, `d` (`kill_rate`) the probability
#' that an encounter ends in a kill and `e` (`er`) the coefficient converting
#' game density into encounters (Levi et al. 2009, 2011).
#'
#' @section Two methods:
#' * `method = "steady"` (default) uses the closed-form steady state of the
#'   multi-settlement model (Levi et al. 2011). In each cell, logistic growth
#'   balances hunting:
#'   \deqn{N(x,y) = \Big[\,K^{\theta}\Big(1 - \frac{e\, d\, h}{r}\,
#'         \sum_{c} P_c\, \varphi(\delta_c)\Big)\Big]_{+}^{1/\theta}.}
#'   Animals do **not** move between cells in this form: it maps local
#'   depletion but not the resupply of hunted areas by immigration.
#' * `method = "dynamic"` runs the time-stepped model of Levi et al. (2011),
#'   which adds animal movement: each year a fraction
#'   `dispersal` of the animals in a cell moves to each of its four neighbours,
#'   then the population grows logistically and is hunted, for `years` years
#'   starting from carrying capacity everywhere. Remote, lightly hunted cells
#'   then act as *sources* that resupply the hunted *sinks* around
#'   settlements. Edges are reflecting (no animals leave the grid), a choice
#'   of this package; the implementation published with Levi et al. (2011)
#'   handles the grid edges differently, which only matters close to them.
#'
#' @section Study area:
#' Areas and fractions are computed inside the study area, so the verdict
#' does not depend on an arbitrary map size. Give the study area as `extent`
#' (a rectangle) or `boundary` (a polygon). The model itself runs on a grid
#' that extends `margin` beyond the study area, so that the hunting grounds of
#' settlements near its edge, and animal movement across it, are not cut off.
#' Without a study area, the grid is the bounding box of the settlements plus
#' `margin`, a warning is issued, and `frac_extirpated` then changes with
#' `margin`.
#'
#' @param data A data frame of settlements, **one row per settlement**. Required
#'   columns: the settlement coordinates (`x`, `y`) in the same length unit as
#'   `sigma` (e.g. kilometres), and the number of hunters (`humans`). The
#'   biological and hunting parameters (`k`, `r`, `hphy`, `kill_rate`, `sigma`,
#'   ...) are passed as scalar arguments, not columns, because they describe the
#'   target species and hunting technology, not individual settlements.
#' @param x,y <[`data-masked`][rlang::args_data_masking]> Settlement
#'   coordinate columns, in the same length unit as `sigma`.
#' @param humans <[`data-masked`][rlang::args_data_masking]> Column giving the
#'   number of hunters (human population) at each settlement.
#' @param k Carrying-capacity density `K` (individuals per unit area).
#' @param r Intrinsic population growth rate (per year).
#' @param hphy Hunts per hunter per year (`h`).
#' @param kill_rate Per-encounter kill probability `d` (weapon efficiency, 0-1).
#' @param sigma Spatial spread of hunting effort \eqn{\sigma} (mean hunting
#'   distance from a settlement), in the coordinate unit.
#' @param er Coefficient converting game density to encounter rate `e`
#'   (default `0.02`, as in Levi et al. 2011).
#' @param theta Density-dependence exponent \eqn{\theta} (default `1`, logistic).
#' @param method `"steady"` (closed-form steady state, no movement) or
#'   `"dynamic"` (time-stepped, with movement between cells). See the Two
#'   methods section.
#' @param extent Optional study area as a rectangle,
#'   `c(xmin, xmax, ymin, ymax)` in the coordinate unit.
#' @param boundary Optional study area as a polygon: a data frame or matrix
#'   whose first two columns are the x and y coordinates of its vertices, in
#'   order. Takes precedence over `extent`. (From an `sf` polygon, use
#'   `sf::st_coordinates(polygon)[, 1:2]`.)
#' @param resolution Grid cell size in coordinate units (default `1`).
#' @param margin Buffer, in coordinate units, by which the simulation grid
#'   extends beyond the study area (default `3 * sigma`, rounded up to whole
#'   cells). When neither `extent` nor `boundary` is given, the study area itself
#'   is the settlements' bounding box plus `margin`.
#' @param years Number of years simulated by `method = "dynamic"` (default
#'   `100`).
#' @param dispersal Share of the animals of a cell moving to each of its four
#'   neighbours every year, between 0 and 0.25 (default `0.02`, as in Levi et
#'   al. 2011). It is defined per cell, so it depends on
#'   `resolution`. Used only by `method = "dynamic"`.
#' @param extirpation_threshold Density at or below which a cell is counted as
#'   locally extirpated (default `1`).
#' @param max_extirpated Maximum fraction of the study area allowed to be
#'   extirpated for the assessment to be flagged `sustainable` (default `0.05`).
#'
#' @return A one-row [offtake][ot_biode] tibble summarising the study area,
#'   with the columns:
#' \describe{
#'   \item{n_settlements}{Number of settlements in `data`.}
#'   \item{area}{Size of the study area, in squared coordinate units (km^2 when
#'     coordinates are in km): the number of grid cells whose centre lies in
#'     the study area times the cell area, so it approximates the true area at
#'     the resolution of the grid.}
#'   \item{mean_density}{Mean predicted game density over the study area.}
#'   \item{min_density}{Lowest predicted density (0 where a sink is fully
#'     extirpated).}
#'   \item{area_extirpated}{Area where the density is at or below
#'     `extirpation_threshold`, in the unit of `area`.}
#'   \item{frac_extirpated}{`area_extirpated / area`: the share of the study
#'     area locally wiped out.}
#'   \item{frac_below_half_k}{Share of the study area with density below half
#'     of carrying capacity `k`, a broader depletion footprint.}
#'   \item{cpue}{Expected kills per hunt, averaged over all hunts on the grid:
#'     the effort-weighted mean density times `er * kill_rate`. This definition
#'     is specific to this package and differs in scale and weighting from the
#'     catch per unit effort reported by Levi et al. (2011), so compare
#'     settlements with each other rather than with values from that study.}
#'   \item{sustainable}{Logical verdict: `TRUE` when
#'     `frac_extirpated <= max_extirpated`.}
#' }
#'   The full predicted density surface and per-settlement CPUE are attached as
#'   attributes and are most easily retrieved with [ot_biode_surface()] and
#'   [ot_biode_cpue()].
#'
#' @references
#' Levi, T., Shepard, G. H., Ohl-Schacherer, J., Peres, C. A. & Yu, D. W. (2009)
#' Modelling the long-term sustainability of indigenous hunting in Manu National
#' Park, Peru: landscape-scale management implications for Amazonia.
#' *Journal of Applied Ecology* 46, 804-814.
#'
#' Levi, T., Shepard, G. H., Ohl-Schacherer, J., Wilmers, C. C., Peres, C. A. &
#' Yu, D. W. (2011) Spatial tools for modeling the sustainability of
#' subsistence hunting in tropical forests. *Ecological Applications* 21,
#' 1802-1818.
#'
#' Adounke, G. R. M. et al. (2026) Systematic review of sustainability
#' assessment approaches for wildlife exploitation. *Biological Conservation*
#' 313, 111606. \doi{10.1016/j.biocon.2025.111606}
#'
#' @seealso [ot_biode_surface()], [ot_biode_cpue()]
#' @examples
#' # One row per settlement. `xkm`/`ykm` = location (km), `hunters` = people.
#' settlements <- data.frame(
#'   village = c("A", "B"),
#'   xkm = c(10, 30),
#'   ykm = c(15, 20),
#'   hunters = c(80, 120)
#' )
#' # Study area: a 50 x 40 km rectangle
#' res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'                 k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'                 extent = c(-5, 45, -5, 35), resolution = 2)
#' res
#' head(ot_biode_surface(res))
#'
#' # With animal movement between cells
#' ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'          k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'          extent = c(-5, 45, -5, 35), resolution = 2, method = "dynamic")
#' @export
ot_biode <- function(data, x, y, humans, k, r, hphy, kill_rate, sigma,
                     er = 0.02, theta = 1, method = c("steady", "dynamic"),
                     extent = NULL, boundary = NULL, resolution = 1,
                     margin = 3 * sigma, years = 100, dispersal = 0.02,
                     extirpation_threshold = 1, max_extirpated = 0.05) {
  .check_data(data)
  method <- match.arg(method)
  data <- .ungroup(data)
  xs <- .pull_col(data, rlang::enquo(x), "x")
  ys <- .pull_col(data, rlang::enquo(y), "y")
  pop <- .pull_col(data, rlang::enquo(humans), "humans")
  if (length(xs) == 0L) rlang::abort("`data` has no settlements.", class = "offtake_bad_data")
  .check_range(pop, "humans", 0)
  .check_range(k, "k", 0, lower_open = TRUE)
  .check_range(r, "r", 0, lower_open = TRUE)
  .check_range(hphy, "hphy", 0)
  .check_range(kill_rate, "kill_rate", 0, 1)
  .check_range(sigma, "sigma", 0, lower_open = TRUE)
  .check_range(er, "er", 0)
  .check_range(theta, "theta", 0, lower_open = TRUE)
  .check_range(resolution, "resolution", 0, lower_open = TRUE)
  .check_range(years, "years", 1)
  .check_range(dispersal, "dispersal", 0, 0.25)
  .check_range(max_extirpated, "max_extirpated", 0, 1)

  .check_range(margin, "margin", 0)

  # Study area. The simulation grid extends `margin` beyond it, so that hunting
  # grounds and animal movement are not cut at its edge; summaries are then
  # computed inside the study area only.
  if (!is.null(boundary)) {
    boundary <- as.matrix(boundary)
    if (ncol(boundary) < 2 || nrow(boundary) < 3) {
      rlang::abort("`boundary` must have at least 3 vertices and 2 columns (x, y).",
                   class = "offtake_bad_value")
    }
    study <- c(range(boundary[, 1]), range(boundary[, 2]))
    buffer <- ceiling(margin / resolution) * resolution
  } else if (!is.null(extent)) {
    .check_range(extent, "extent")
    if (length(extent) != 4 || extent[1] >= extent[2] || extent[3] >= extent[4]) {
      rlang::abort("`extent` must be c(xmin, xmax, ymin, ymax) with xmin < xmax and ymin < ymax.",
                   class = "offtake_bad_value")
    }
    study <- extent
    buffer <- ceiling(margin / resolution) * resolution
  } else {
    study <- c(min(xs) - margin, max(xs) + margin, min(ys) - margin, max(ys) + margin)
    buffer <- 0
    rlang::warn(
      "No study area given (`extent` or `boundary`): areas and `frac_extirpated` refer to the settlements' bounding box plus `margin`, and change with it.",
      class = "offtake_no_study_area",
      .frequency = "once",
      .frequency_id = "offtake_biode_study_area"
    )
  }
  gx <- seq(study[1] - buffer + resolution / 2, study[2] + buffer, by = resolution)
  gy <- seq(study[3] - buffer + resolution / 2, study[4] + buffer, by = resolution)
  nx <- length(gx); ny <- length(gy)
  cell_x <- outer(rep(1, ny), gx)
  cell_y <- outer(gy, rep(1, nx))
  # Cells of the study area (its bounding box, then the polygon if any).
  keep_x <- gx >= study[1] & gx <= study[2]
  keep_y <- gy >= study[3] & gy <= study[4]
  inside <- outer(keep_y, keep_x, `&`)
  if (!is.null(boundary)) {
    inside <- inside & matrix(.in_polygon(as.vector(cell_x), as.vector(cell_y),
                                          boundary[, 1], boundary[, 2]), ny, nx)
  }
  if (!any(inside)) {
    rlang::abort("No grid cell falls inside the study area; check `boundary` and `resolution`.",
                 class = "offtake_bad_value")
  }

  # Hunting pressure: sum over settlements of P_c * phi(distance).
  kernel <- function(c) {
    d2 <- (cell_x - xs[c])^2 + (cell_y - ys[c])^2
    exp(-d2 / (2 * sigma^2)) / (2 * pi * sqrt(d2) + 1)
  }
  pressure <- matrix(0, ny, nx)
  for (c in seq_along(xs)) pressure <- pressure + pop[c] * kernel(c)

  dens <- if (method == "steady") {
    pmax(k^theta * (1 - er * kill_rate * hphy / r * pressure), 0)^(1 / theta)
  } else {
    .biode_dynamic(er * kill_rate * hphy * pressure, k, r, theta, dispersal, years)
  }

  # Summaries inside the study area.
  cell_area <- resolution^2
  area <- sum(inside) * cell_area
  area_ext <- sum(inside & dens <= extirpation_threshold) * cell_area
  frac_ext <- area_ext / area
  # Expected kills per hunt, weighting cells by where the hunts happen.
  cpue_global <- if (sum(pressure) > 0) {
    sum(dens * pressure) / sum(pressure) * er * kill_rate
  } else NA_real_
  local_cpue <- vapply(seq_along(xs), function(c) {
    w <- pop[c] * kernel(c)
    if (sum(w) > 0) sum(dens * w) / sum(w) * er * kill_rate else NA_real_
  }, numeric(1))

  summary_tbl <- tibble::tibble(
    n_settlements = length(xs),
    area = area,
    mean_density = mean(dens[inside]),
    min_density = min(dens[inside]),
    area_extirpated = area_ext,
    frac_extirpated = frac_ext,
    frac_below_half_k = mean(dens[inside] < k / 2),
    cpue = cpue_global,
    sustainable = frac_ext <= max_extirpated
  )

  # Outputs cover the study area only (its bounding box, NA outside a polygon).
  dens_out <- dens
  dens_out[!inside] <- NA
  dens_out <- dens_out[keep_y, keep_x, drop = FALSE]
  surface <- tibble::tibble(
    x = as.vector(cell_x[keep_y, keep_x, drop = FALSE]),
    y = as.vector(cell_y[keep_y, keep_x, drop = FALSE]),
    density = as.vector(dens_out)
  )
  cpue_tbl <- tibble::tibble(
    settlement = seq_along(xs),
    x = xs, y = ys, humans = pop,
    local_cpue = local_cpue
  )

  out <- new_offtake(
    summary_tbl,
    method = sprintf("Biode (Levi et al. spatial model, %s)", method),
    family = "model",
    reference = "Levi et al. (2009)"
  )
  attr(out, "surface") <- surface
  attr(out, "grid") <- list(x = gx[keep_x], y = gy[keep_y], density = dens_out)
  attr(out, "cpue") <- cpue_tbl
  attr(out, "params") <- list(k = k, r = r, hphy = hphy, kill_rate = kill_rate,
                              sigma = sigma, er = er, theta = theta,
                              method = method, resolution = resolution,
                              years = years, dispersal = dispersal)
  out
}

# Time-stepped model of Levi et al. (2011): movement to the four
# neighbours, logistic growth and hunting, starting from K everywhere.
# `kill` is the annual per-capita kill rate of each cell.
.biode_dynamic <- function(kill, k, r, theta, dispersal, years) {
  ny <- nrow(kill); nx <- ncol(kill)
  n <- matrix(k, ny, nx)
  for (t in seq_len(years)) {
    # Reflecting edges: a cell on the edge exchanges with itself.
    up <- n[c(1, seq_len(ny - 1)), , drop = FALSE]
    down <- n[c(seq_len(ny)[-1], ny), , drop = FALSE]
    left <- n[, c(1, seq_len(nx - 1)), drop = FALSE]
    right <- n[, c(seq_len(nx)[-1], nx), drop = FALSE]
    moved <- dispersal * (up + down + left + right) + (1 - 4 * dispersal) * n
    growth <- r * n * (1 - (n / k)^theta)
    n <- pmax(moved - kill * n + growth, 0) # matrix first, so dim is kept
  }
  n
}

# Point-in-polygon test (ray casting), vectorised over points.
.in_polygon <- function(px, py, vx, vy) {
  inside <- rep(FALSE, length(px))
  j <- length(vx)
  for (i in seq_along(vx)) {
    crosses <- ((vy[i] > py) != (vy[j] > py)) &
      (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])
    inside <- xor(inside, crosses %in% TRUE)
    j <- i
  }
  inside
}

#' Extract the predicted density surface from a Biode result
#'
#' @param x An [offtake][ot_biode] object returned by [ot_biode()].
#' @return A tibble with **one row per grid cell** and the columns:
#' \describe{
#'   \item{x}{Cell-centre x coordinate (coordinate unit of the settlements).}
#'   \item{y}{Cell-centre y coordinate.}
#'   \item{density}{Predicted game density in that cell; `NA` outside the
#'     study area when a `boundary` was given.}
#' }
#'   This long format is ready for mapping, e.g. with
#'   `ggplot2::geom_raster(ggplot2::aes(x, y, fill = density))`.
#' @seealso [ot_biode()]
#' @examples
#' settlements <- data.frame(xkm = 10, ykm = 10, hunters = 100)
#' res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'                 k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'                 extent = c(-10, 30, -10, 30), resolution = 2)
#' ot_biode_surface(res)
#' @export
ot_biode_surface <- function(x) {
  s <- attr(x, "surface")
  if (is.null(s)) rlang::abort("`x` is not an ot_biode() result.", class = "offtake_bad_object")
  s
}

#' Extract per-settlement catch per unit effort from a Biode result
#'
#' @param x An [offtake][ot_biode] object returned by [ot_biode()].
#' @return A tibble with **one row per settlement** and the columns:
#' \describe{
#'   \item{settlement}{Row index of the settlement in the input data.}
#'   \item{x, y}{Settlement coordinates.}
#'   \item{humans}{Number of hunters at the settlement.}
#'   \item{local_cpue}{Expected kills per hunt for that settlement's hunters:
#'     the density averaged over where they hunt (weighted by their own effort
#'     kernel) times `er * kill_rate`. Lower values indicate a more depleted
#'     hunting ground. See the `cpue` column of [ot_biode()] for how this
#'     relates to the catch per unit effort of Levi et al. (2011).}
#' }
#' @seealso [ot_biode()]
#' @examples
#' settlements <- data.frame(xkm = c(10, 30), ykm = c(10, 20), hunters = c(80, 120))
#' res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'                 k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'                 extent = c(-10, 50, -10, 40), resolution = 2)
#' ot_biode_cpue(res)
#' @export
ot_biode_cpue <- function(x) {
  s <- attr(x, "cpue")
  if (is.null(s)) rlang::abort("`x` is not an ot_biode() result.", class = "offtake_bad_object")
  s
}
