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
#' @section Units:
#' Distances must be in **kilometres** and densities per km^2: `x`, `y`,
#' `sigma`, `resolution`, `margin`, `extent` and `boundary` in km, `k` and
#' `extirpation_threshold` in individuals per km^2. The hunting kernel is not
#' scale-free (the `+ 1` of its denominator is one kilometre) and the encounter
#' coefficient `er` is defined for densities per km^2, as in Levi et al.
#' (2009, 2011), whose grid has 1 km cells. Projected coordinates in metres,
#' such as UTM, must therefore be divided by 1000. `mask` and `no_hunt` are the
#' exception: they may stay in a projected coordinate system in metres, and
#' the cell centres are converted to metres to be matched with them.
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
#' @section Excluded and unhunted areas (mask, no_hunt):
#' By default every cell is treated as habitat with the same carrying capacity.
#' `mask` removes the areas that are not habitat for the target species and
#' where it is therefore not hunted, such as open water, fields or, for a
#' forest species, open woodland. This option is an extension of this package;
#' it is not part of the model of Levi et al. (2009, 2011). Masked cells:
#' * hold no animals, and in `method = "dynamic"` animals do not move into
#'   them (they act like the reflecting edges of the grid);
#' * are left out of `area`, of every summary and of the verdict, and are `NA`
#'   in [ot_biode_surface()];
#' * receive no hunts: the hunting pressure on habitat cells is unchanged, and
#'   the hunts that the kernel would place in masked cells are not made (they
#'   are not moved elsewhere). Catch per unit effort is therefore averaged over
#'   habitat cells only.
#'
#' `no_hunt` describes the opposite case: habitat where the species lives but
#' is not hunted, such as a strict reserve, a sacred forest or an area hunters
#' cannot reach. These cells hold animals that grow and, in
#' `method = "dynamic"`, move to and from neighbouring cells, but receive no
#' hunting pressure, so they can resupply the hunted areas around them. They
#' are habitat, so they count in `area`, in every summary and in the verdict.
#' A cell in both `mask` and `no_hunt` is masked. Like `mask`, this option is
#' an extension of this package.
#'
#' Both apply to the whole simulation grid, including the buffer around the
#' study area. They must use a projected coordinate system, in kilometres or in
#' metres (see the Units section), whose origin is the one used for `x` and
#' `y`.
#'
#' @param data A data frame of settlements, **one row per settlement**. Required
#'   columns: the settlement coordinates (`x`, `y`) in kilometres, and the
#'   number of hunters (`humans`). The
#'   biological and hunting parameters (`k`, `r`, `hphy`, `kill_rate`, `sigma`,
#'   ...) are passed as scalar arguments, not columns, because they describe the
#'   target species and hunting technology, not individual settlements.
#' @param x,y <[`data-masked`][rlang::args_data_masking]> Settlement
#'   coordinate columns, in kilometres (see the Units section).
#' @param humans <[`data-masked`][rlang::args_data_masking]> Column giving the
#'   number of hunters (human population) at each settlement.
#' @param k Carrying-capacity density `K` (individuals per km^2).
#' @param r Intrinsic population growth rate (per year).
#' @param hphy Hunts per hunter per year (`h`).
#' @param kill_rate Per-encounter kill probability `d` (weapon efficiency, 0-1).
#' @param sigma Spatial spread of hunting effort \eqn{\sigma} (mean hunting
#'   distance from a settlement), in km.
#' @param er Coefficient converting game density to encounter rate `e`
#'   (default `0.02`, as in Levi et al. 2011).
#' @param theta Density-dependence exponent \eqn{\theta} (default `1`, logistic).
#' @param method `"steady"` (closed-form steady state, no movement) or
#'   `"dynamic"` (time-stepped, with movement between cells). See the Two
#'   methods section.
#' @param extent Optional study area as a rectangle,
#'   `c(xmin, xmax, ymin, ymax)` in km.
#' @param boundary Optional study area as a polygon: a data frame or matrix
#'   whose first two columns are the x and y coordinates of its vertices, in
#'   order, in km. Takes precedence over `extent`. (From an `sf` polygon in a
#'   projected coordinate system in metres, use
#'   `sf::st_coordinates(polygon)[, 1:2] / 1000`.)
#' @param mask Optional areas to exclude from the habitat (see the Excluded
#'   and unhunted areas section), in the projected coordinate system of `x`
#'   and `y`, in km or in metres (see the Units section). Either polygons
#'   covering the excluded areas, as an `sf`/`sfc` object or a `terra`
#'   `SpatVector`, or a `terra` `SpatRaster` whose cells are `TRUE` (or
#'   non-zero) where the area is excluded; `NA` and `0` cells, and places the
#'   raster does not cover, are kept. A cell is excluded when its centre falls
#'   in an excluded area. From a land-cover raster, build the mask with a
#'   comparison such as `landcover == 1 | landcover == 4`.
#' @param no_hunt Optional habitat areas where the species is present but not
#'   hunted (see the Excluded and unhunted areas section), in the same forms as
#'   `mask`: polygons covering these areas, or a raster that is `TRUE` (or
#'   non-zero) in them.
#' @param resolution Grid cell size in km (default `1`, the cell size of Levi
#'   et al. 2011).
#' @param margin Buffer, in km, by which the simulation grid
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
#'   locally extirpated, in individuals per km^2 (default `1`).
#' @param max_extirpated Maximum fraction of the study area allowed to be
#'   extirpated for the assessment to be flagged `sustainable` (default `0.05`).
#'
#' @return A one-row [offtake][ot_biode] tibble summarising the study area,
#'   with the columns:
#' \describe{
#'   \item{n_settlements}{Number of settlements in `data`.}
#'   \item{area}{Habitat area of the study area, in km^2: the number of grid cells whose centre
#'     lies in the study area, and outside `mask`, times the cell area, so it
#'     approximates the true area at the resolution of the grid.}
#'   \item{area_masked}{Area of the study area removed by `mask` (0 without a
#'     mask), in the unit of `area`.}
#'   \item{area_no_hunt}{Habitat area of the study area set by `no_hunt` as
#'     unhunted (0 without `no_hunt`), in the unit of `area`. It is included in
#'     `area`.}
#'   \item{mean_density}{Mean predicted game density over the habitat of the
#'     study area.}
#'   \item{min_density}{Lowest predicted density (0 where a sink is fully
#'     extirpated).}
#'   \item{area_extirpated}{Area where the density is at or below
#'     `extirpation_threshold`, in the unit of `area`.}
#'   \item{frac_extirpated}{`area_extirpated / area`: the share of the habitat
#'     locally wiped out.}
#'   \item{frac_below_half_k}{Share of the habitat with density below half
#'     of carrying capacity `k`, a broader depletion footprint.}
#'   \item{cpue}{Expected kills per hunt, averaged over all hunts made on the
#'     habitat cells of the grid: the effort-weighted mean density times
#'     `er * kill_rate`. This definition is specific to this package and
#'     differs in scale and weighting from the
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
                     extent = NULL, boundary = NULL, mask = NULL, no_hunt = NULL,
                     resolution = 1, margin = 3 * sigma, years = 100, dispersal = 0.02,
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
  # Habitat cells: everything outside `mask`.
  habitat <- if (is.null(mask)) {
    matrix(TRUE, ny, nx)
  } else {
    !matrix(.masked_cells(mask, as.vector(cell_x), as.vector(cell_y)), ny, nx)
  }
  # Unhunted habitat: animals, but no hunting pressure.
  hunted <- if (is.null(no_hunt)) {
    matrix(TRUE, ny, nx)
  } else {
    !matrix(.masked_cells(no_hunt, as.vector(cell_x), as.vector(cell_y), "no_hunt"),
            ny, nx)
  }
  if (!any(inside & habitat)) {
    rlang::abort("`mask` excludes every cell of the study area.", class = "offtake_bad_value")
  }

  # Hunting pressure: sum over settlements of P_c * phi(distance).
  kernel <- function(c) {
    d2 <- (cell_x - xs[c])^2 + (cell_y - ys[c])^2
    hunted * exp(-d2 / (2 * sigma^2)) / (2 * pi * sqrt(d2) + 1)
  }
  pressure <- matrix(0, ny, nx)
  for (c in seq_along(xs)) pressure <- pressure + pop[c] * kernel(c)

  dens <- if (method == "steady") {
    pmax(k^theta * (1 - er * kill_rate * hphy / r * pressure), 0)^(1 / theta)
  } else {
    .biode_dynamic(er * kill_rate * hphy * pressure, k, r, theta, dispersal, years,
                   habitat)
  }
  dens[!habitat] <- NA

  # Summaries inside the study area.
  cell_area <- resolution^2
  study_hab <- inside & habitat
  area <- sum(study_hab) * cell_area
  area_ext <- sum(study_hab & dens <= extirpation_threshold) * cell_area
  frac_ext <- area_ext / area
  # Expected kills per hunt, weighting habitat cells by where the hunts happen
  # (no hunt is made in masked or unhunted cells).
  cpue_of <- function(w) {
    w <- w[habitat]
    if (sum(w) > 0) sum(dens[habitat] * w) / sum(w) * er * kill_rate else NA_real_
  }
  cpue_global <- cpue_of(pressure)
  local_cpue <- vapply(seq_along(xs), function(c) cpue_of(pop[c] * kernel(c)), numeric(1))

  summary_tbl <- tibble::tibble(
    n_settlements = length(xs),
    area = area,
    area_masked = sum(inside & !habitat) * cell_area,
    area_no_hunt = sum(study_hab & !hunted) * cell_area,
    mean_density = mean(dens[study_hab]),
    min_density = min(dens[study_hab]),
    area_extirpated = area_ext,
    frac_extirpated = frac_ext,
    frac_below_half_k = mean(dens[study_hab] < k / 2),
    cpue = cpue_global,
    sustainable = frac_ext <= max_extirpated
  )

  # Outputs cover the study area only (its bounding box, NA outside a polygon
  # and in masked cells).
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
                              years = years, dispersal = dispersal,
                              masked = !is.null(mask),
                              no_hunt = !is.null(no_hunt))
  out
}

# Time-stepped model of Levi et al. (2011): movement to the four
# neighbours, logistic growth and hunting, starting from K everywhere.
# `kill` is the annual per-capita kill rate of each cell. `habitat` is a
# logical matrix; non-habitat cells hold no animals and receive none.
.biode_dynamic <- function(kill, k, r, theta, dispersal, years,
                           habitat = matrix(TRUE, nrow(kill), ncol(kill))) {
  ny <- nrow(kill); nx <- ncol(kill)
  # Reflecting edges: a cell on the edge exchanges with itself.
  rows_up <- c(1, seq_len(ny - 1))
  rows_down <- c(seq_len(ny)[-1], ny)
  cols_left <- c(1, seq_len(nx - 1))
  cols_right <- c(seq_len(nx)[-1], nx)
  h <- habitat * 1
  # Number of neighbours each cell can send animals to; movement towards a
  # non-habitat cell is cancelled, as at an edge.
  open <- h[rows_up, , drop = FALSE] + h[rows_down, , drop = FALSE] +
    h[, cols_left, drop = FALSE] + h[, cols_right, drop = FALSE]
  n <- k * h
  for (t in seq_len(years)) {
    up <- n[rows_up, , drop = FALSE]
    down <- n[rows_down, , drop = FALSE]
    left <- n[, cols_left, drop = FALSE]
    right <- n[, cols_right, drop = FALSE]
    moved <- h * (dispersal * (up + down + left + right) + (1 - open * dispersal) * n)
    growth <- r * n * (1 - (n / k)^theta)
    n <- pmax(moved - kill * n + growth, 0) # matrix first, so dim is kept
  }
  n
}

# Which points (cell centres) fall in the areas given by `mask` (or `no_hunt`,
# named by `arg`): sf/sfc or terra SpatVector polygons, or a terra SpatRaster
# that is TRUE/non-zero in these areas.
.masked_cells <- function(mask, px, py, arg = "mask") {
  bad_crs <- function() {
    rlang::abort(c(sprintf("`%s` uses longitude/latitude coordinates.", arg),
                   i = "Project it, and the settlements, to a coordinate system in the unit of `sigma` (e.g. km or m)."),
                 class = "offtake_bad_value")
  }
  # Cell centres are in km; express them in the unit of the mask's coordinate
  # system (metres for most projected systems, e.g. UTM).
  to_mask_unit <- function(metres_per_unit) {
    if (is.na(metres_per_unit) || metres_per_unit <= 0) return(1)
    1000 / metres_per_unit
  }
  hits <- if (inherits(mask, c("sf", "sfc"))) {
    rlang::check_installed("sf", sprintf("to use an sf object as `%s`.", arg))
    crs <- sf::st_crs(mask)
    if (isTRUE(sf::st_is_longlat(mask))) bad_crs()
    unit <- if (is.na(crs)) NULL else crs$units_gdal
    f <- if (is.null(unit) || is.na(unit)) 1 else switch(
      unit, metre = 1000, kilometre = 1,
      rlang::abort(sprintf("The coordinate system of `%s` is in %s; use metres or kilometres.",
                           arg, unit), class = "offtake_bad_value")
    )
    bb <- sf::st_bbox(mask)
    pts <- sf::st_as_sf(data.frame(x = px * f, y = py * f), coords = c("x", "y"), crs = crs)
    lengths(sf::st_intersects(pts, sf::st_geometry(mask))) > 0
  } else if (inherits(mask, "SpatVector")) {
    rlang::check_installed("terra", sprintf("to use a SpatVector as `%s`.", arg))
    if (terra::crs(mask) != "" && isTRUE(terra::is.lonlat(mask))) bad_crs()
    f <- if (terra::crs(mask) == "") 1 else to_mask_unit(terra::linearUnits(mask))
    bb <- as.vector(terra::ext(mask))[c(1, 3, 2, 4)]
    pts <- terra::vect(cbind(px * f, py * f), crs = terra::crs(mask))
    terra::is.related(pts, mask, "intersects")
  } else if (inherits(mask, "SpatRaster")) {
    rlang::check_installed("terra", sprintf("to use a SpatRaster as `%s`.", arg))
    if (terra::crs(mask) != "" && isTRUE(terra::is.lonlat(mask))) bad_crs()
    if (isTRUE(terra::is.factor(mask)[1])) {
      rlang::abort(c(sprintf("`%s` is a categorical raster.", arg),
                     i = "Build a TRUE/FALSE raster of the classes concerned, e.g. `landcover == 1 | landcover == 4`."),
                   class = "offtake_bad_value")
    }
    f <- if (terra::crs(mask) == "") 1 else to_mask_unit(terra::linearUnits(mask))
    bb <- as.vector(terra::ext(mask))[c(1, 3, 2, 4)]
    v <- as.numeric(terra::extract(mask[[1]], cbind(px * f, py * f))[, 1])
    !is.na(v) & v != 0
  } else {
    rlang::abort(sprintf("`%s` must be an sf or sfc object, a terra SpatVector or a terra SpatRaster.", arg),
                 class = "offtake_bad_value")
  }
  # A mask that does not overlap the grid usually means mismatched units or
  # coordinate systems.
  bb <- as.numeric(bb)
  if (!any(px * f >= bb[1] & px * f <= bb[3] & py * f >= bb[2] & py * f <= bb[4])) {
    rlang::warn(c(sprintf("`%s` does not overlap the simulation grid.", arg),
                  i = "Check that `x` and `y` are in km and use the same projection as the mask."),
                class = "offtake_mask_no_overlap")
  }
  hits
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
#'   \item{x}{Cell-centre x coordinate, in km.}
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
