#' Biodemographic spatial source-sink model (Biode)
#'
#' Predicts the steady-state spatial distribution of game density across a
#' landscape given the location and human population of hunting settlements,
#' following the spatially explicit reaction-diffusion source-sink model of
#' Levi et al. (2009). Unlike the non-spatial models ([ot_pro()], [ot_pbr()],
#' [ot_msy()]), Biode maps where hunting creates population *sinks* around
#' settlements and *sources* in the remote matrix.
#'
#' At steady state the model reduces to a closed form for the game density at
#' each point of the landscape:
#' \deqn{N(x,y) = \Big[\,K^{\theta}\Big(1 - \frac{e\, d\, h}{r}\,
#'        \sum_{c} P_c\, \varphi(\delta_c)\Big)\Big]_{+}^{1/\theta},}
#' with the Gaussian hunting-effort kernel
#' \deqn{\varphi(\delta) = \frac{\exp\!\big(-\delta^2 / 2\sigma^2\big)}
#'        {2\pi\,\delta + 1},}
#' where \eqn{\delta_c} is the distance from a point to settlement `c`, `P_c`
#' its hunter population, `K` carrying capacity, `r` the intrinsic growth rate,
#' `h` (`hphy`) hunts per hunter per year, `d` (`kill_rate`) the per-encounter
#' kill probability (weapon efficiency), `e` (`er`) the density-to-encounter
#' coefficient, \eqn{\theta} the density-dependence exponent, and
#' \eqn{[\cdot]_+} denotes truncation at zero. The form is taken directly from
#' Levi et al. (2009) and their published solver code.
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
#' @param k Carrying-capacity density `K` (individuals per cell area).
#' @param r Intrinsic population growth rate.
#' @param hphy Hunts per hunter per year (`h`).
#' @param kill_rate Per-encounter kill probability `d` (weapon efficiency, 0-1).
#' @param sigma Spatial spread of hunting effort \eqn{\sigma} (mean hunting
#'   distance from a settlement), in the coordinate unit.
#' @param er Coefficient converting game density to encounter rate `e`
#'   (default `0.02`, the value used by Levi et al. 2009).
#' @param theta Density-dependence exponent \eqn{\theta} (default `1`, logistic).
#' @param resolution Grid cell size in coordinate units (default `1`).
#' @param margin Landscape margin added around the settlement bounding box, in
#'   coordinate units (default `3 * sigma`).
#' @param extirpation_threshold Density at or below which a cell is counted as
#'   locally extirpated (default `1`).
#' @param max_extirpated Maximum fraction of the landscape allowed to be
#'   extirpated for the assessment to be flagged `sustainable` (default `0.05`).
#'
#' @return A one-row [offtake][ot_biode] tibble summarising the landscape, with
#'   the columns:
#' \describe{
#'   \item{n_settlements}{Number of settlements in `data`.}
#'   \item{mean_density}{Mean predicted game density across all grid cells.}
#'   \item{min_density}{Lowest predicted density (0 where a sink is fully
#'     extirpated).}
#'   \item{frac_extirpated}{Fraction of grid cells at or below
#'     `extirpation_threshold` -- the share of the landscape locally wiped out.}
#'   \item{frac_below_half_k}{Fraction of cells with density below half of
#'     carrying capacity `k` -- a broader depletion footprint.}
#'   \item{cpue}{Effort-weighted mean catch per unit effort predicted across the
#'     landscape.}
#'   \item{sustainable}{Logical verdict: `TRUE` when
#'     `frac_extirpated <= max_extirpated`.}
#' }
#'   The full predicted density surface and per-settlement local CPUE are
#'   attached as attributes and are most easily retrieved with
#'   [ot_biode_surface()] and [ot_biode_cpue()].
#'
#' @references
#' Levi, T., Shepard, G. H., Ohl-Schacherer, J., Peres, C. A. & Yu, D. W. (2009)
#' Modelling the long-term sustainability of indigenous hunting in Manu National
#' Park, Peru: landscape-scale management implications for Amazonia.
#' *Journal of Applied Ecology* 46, 804-814.
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
#' res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'                 k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'                 resolution = 2)
#' res
#' head(ot_biode_surface(res))
#' @export
ot_biode <- function(data, x, y, humans, k, r, hphy, kill_rate, sigma,
                     er = 0.02, theta = 1, resolution = 1, margin = 3 * sigma,
                     extirpation_threshold = 1, max_extirpated = 0.05) {
  .check_data(data)
  xs <- .pull_col(data, rlang::enquo(x), "x")
  ys <- .pull_col(data, rlang::enquo(y), "y")
  pop <- .pull_col(data, rlang::enquo(humans), "humans")

  if (length(xs) == 0L) rlang::abort("`data` has no settlements.", class = "offtake_bad_data")

  # Landscape grid (cell centres).
  gx <- seq(min(xs) - margin, max(xs) + margin, by = resolution)
  gy <- seq(min(ys) - margin, max(ys) + margin, by = resolution)
  nx <- length(gx); ny <- length(gy)

  coef <- er * kill_rate * hphy / r

  # Accumulate hunting pressure = sum_c P_c * phi(distance) over settlements.
  pressure <- matrix(0, nrow = ny, ncol = nx)
  for (c in seq_along(xs)) {
    dx <- outer(rep(1, ny), gx) - xs[c]           # ny x nx difference in x
    dy <- outer(gy, rep(1, nx)) - ys[c]           # ny x nx difference in y
    d2 <- dx^2 + dy^2
    phi <- exp(-d2 / (2 * sigma^2)) / (2 * pi * sqrt(d2) + 1)
    pressure <- pressure + pop[c] * phi
  }

  dens <- pmax(0, k^theta * (1 - coef * pressure))^(1 / theta)

  frac_ext <- mean(dens <= extirpation_threshold)
  cpue_global <- if (sum(pressure) > 0) {
    sum(dens * pressure) / sum(pressure) * er * kill_rate
  } else NA_real_

  # Per-settlement local CPUE: effort-weighted mean density using that
  # settlement's own kernel (cf. Levi et al. 2009 solver).
  local_cpue <- vapply(seq_along(xs), function(c) {
    dx <- outer(rep(1, ny), gx) - xs[c]
    dy <- outer(gy, rep(1, nx)) - ys[c]
    d2 <- dx^2 + dy^2
    w <- pop[c] * exp(-d2 / (2 * sigma^2)) / (2 * pi * sqrt(d2) + 1)
    if (sum(w) > 0) sum(dens * w) / sum(w) * er * kill_rate else NA_real_
  }, numeric(1))

  summary_tbl <- tibble::tibble(
    n_settlements = length(xs),
    mean_density = mean(dens),
    min_density = min(dens),
    frac_extirpated = frac_ext,
    frac_below_half_k = mean(dens < k / 2),
    cpue = cpue_global,
    sustainable = frac_ext <= max_extirpated
  )

  surface <- tibble::tibble(
    x = rep(gx, each = ny),
    y = rep(gy, times = nx),
    density = as.vector(dens)
  )
  cpue_tbl <- tibble::tibble(
    settlement = seq_along(xs),
    x = xs, y = ys, humans = pop,
    local_cpue = local_cpue
  )

  out <- new_offtake(
    summary_tbl,
    method = "Biode (Levi et al. spatial source-sink model)",
    family = "model",
    reference = "Levi et al. (2009)"
  )
  attr(out, "surface") <- surface
  attr(out, "grid") <- list(x = gx, y = gy, density = dens)
  attr(out, "cpue") <- cpue_tbl
  attr(out, "params") <- list(k = k, r = r, hphy = hphy, kill_rate = kill_rate,
                              sigma = sigma, er = er, theta = theta,
                              resolution = resolution)
  out
}

#' Extract the predicted density surface from a Biode result
#'
#' @param x An [offtake][ot_biode] object returned by [ot_biode()].
#' @return A tibble with **one row per grid cell** and the columns:
#' \describe{
#'   \item{x}{Cell-centre x coordinate (coordinate unit of the settlements).}
#'   \item{y}{Cell-centre y coordinate.}
#'   \item{density}{Predicted steady-state game density in that cell.}
#' }
#'   This long format is ready for mapping, e.g. with
#'   `ggplot2::geom_raster(ggplot2::aes(x, y, fill = density))`.
#' @seealso [ot_biode()]
#' @examples
#' settlements <- data.frame(xkm = 10, ykm = 10, hunters = 100)
#' res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'                 k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'                 resolution = 2)
#' ot_biode_surface(res)
#' @export
ot_biode_surface <- function(x) {
  s <- attr(x, "surface")
  if (is.null(s)) rlang::abort("`x` is not an ot_biode() result.", class = "offtake_bad_object")
  s
}

#' Extract per-settlement local CPUE from a Biode result
#'
#' @param x An [offtake][ot_biode] object returned by [ot_biode()].
#' @return A tibble with **one row per settlement** and the columns:
#' \describe{
#'   \item{settlement}{Row index of the settlement in the input data.}
#'   \item{x, y}{Settlement coordinates.}
#'   \item{humans}{Number of hunters at the settlement.}
#'   \item{local_cpue}{Predicted catch per unit effort in that settlement's
#'     own hunting neighbourhood (effort-weighted mean density times the
#'     encounter and kill coefficients). Lower values indicate a more depleted
#'     local catchment.}
#' }
#' @seealso [ot_biode()]
#' @examples
#' settlements <- data.frame(xkm = c(10, 30), ykm = c(10, 20), hunters = c(80, 120))
#' res <- ot_biode(settlements, x = xkm, y = ykm, humans = hunters,
#'                 k = 25, r = 0.07, hphy = 40, kill_rate = 0.1, sigma = 6,
#'                 resolution = 2)
#' ot_biode_cpue(res)
#' @export
ot_biode_cpue <- function(x) {
  s <- attr(x, "cpue")
  if (is.null(s)) rlang::abort("`x` is not an ot_biode() result.", class = "offtake_bad_object")
  s
}
