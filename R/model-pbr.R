#' Potential biological removal (PBR)
#'
#' Estimates the maximum number of animals that may be removed from a population
#' per year while allowing it to reach or stay near its optimal size, and
#' compares it with the observed human-caused removal. Originally developed for
#' marine mammal management (Wade 1998) and increasingly applied to terrestrial
#' harvest (Adounke et al. 2026).
#'
#' \deqn{PBR = N_{min}\; \tfrac{1}{2} R_{max}\; F_R}
#' where `N_min` is a conservative (minimum) population estimate,
#' \eqn{R_{max}} is the maximum annual net productivity rate (the discrete rate
#' \eqn{\lambda_{max} - 1}) and \eqn{F_R} is a recovery factor between 0.1 and
#' 1. The term \eqn{\tfrac{1}{2} R_{max}} is the per-capita growth rate of a
#' logistic population at half its carrying capacity, where it is most
#' productive. **A removal exceeding PBR indicates unsustainability.**
#'
#' `N_min` may be supplied directly through `nmin`, or computed from an
#' abundance estimate `n` and its coefficient of variation `cv` as the lower
#' tail of a log-normal distribution,
#' \eqn{N_{min} = N / \exp\!\big(z\sqrt{\ln(1+CV^2)}\big)}, with `z` the
#' standard-normal quantile for the chosen percentile (`z = 0.842` for the 20th
#' percentile, as recommended by Wade 1998).
#'
#' Grouped data frames (from `dplyr::group_by()`) are accepted; the group
#' columns are kept in the output.
#'
#' @inheritSection ot_pro Uncertainty
#' @inheritParams ot_pro
#' @param data A data frame with **one row per population/stock**. Required
#'   columns: the maximum growth rate (`rmax`) and the observed annual removal
#'   (`removal`). You must also supply the population size, either as a minimum
#'   estimate (`nmin`) or as a point estimate (`n`) with an optional coefficient
#'   of variation (`cv`) from which `nmin` is derived. May be grouped with
#'   `dplyr::group_by()`.
#' @param rmax <[`data-masked`][rlang::args_data_masking]> Maximum annual net
#'   productivity rate \eqn{R_{max}}, i.e. \eqn{\lambda_{max} - 1} (e.g. `0.04`
#'   for many large mammals, `0.12` for fast-breeding species). An instantaneous
#'   rate `r` can be converted with `rmax = exp(r) - 1`.
#' @param removal <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   human-caused removal (offtake), on the same basis as the population size.
#' @param nmin <[`data-masked`][rlang::args_data_masking]> Minimum population
#'   estimate. Optional if `n` (and optionally `cv`) are supplied.
#' @param n <[`data-masked`][rlang::args_data_masking]> Point abundance estimate,
#'   used with `cv` to derive `nmin`. Use the current abundance, not the
#'   carrying capacity.
#' @param cv <[`data-masked`][rlang::args_data_masking]> Coefficient of variation
#'   of `n` (e.g. `0.3` for a 30% CV). Treated as 0 (so `nmin = n`) when omitted.
#' @param fr Recovery factor \eqn{F_R} in `[0.1, 1]` (default `0.5`). May be a
#'   single value or an embraced column. Use 0.1 for endangered or poorly known
#'   populations.
#' @param percentile Lower percentile used to convert `n`/`cv` to `nmin`
#'   (default `0.20`, giving `z = 0.842`).
#' @param uncertainty Optional named list of coefficients of variation (numbers
#'   or columns of `data`). Allowed names: `n` or `nmin` (whichever is used),
#'   `rmax` and `removal`. If `cv` already describes the uncertainty of `n`, do
#'   not add `n` here as well.
#'
#' @return An [offtake][ot_pbr] tibble with one row per input row and the
#'   columns:
#' \describe{
#'   \item{nmin}{Minimum population estimate used (supplied, or derived from
#'     `n` and `cv`).}
#'   \item{rmax}{Maximum annual net productivity rate used.}
#'   \item{fr}{Recovery factor applied.}
#'   \item{limit}{Potential biological removal, the maximum sustainable annual
#'     removal.}
#'   \item{observed}{The observed annual removal (`removal`).}
#'   \item{ratio}{Exploitation ratio, `observed / limit`.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `observed <= limit`.}
#'   \item{limit_lo, limit_hi, p_unsustainable}{Only with `uncertainty`; see
#'     the Uncertainty section.}
#' }
#'   Group columns come first when `data` is grouped.
#'
#' @references
#' Wade, P. R. (1998) Calculating limits to the allowable human-caused mortality
#' of cetaceans and pinnipeds. *Marine Mammal Science* 14, 1-37.
#' \doi{10.1111/j.1748-7692.1998.tb00688.x}
#'
#' @seealso [ot_pro()], [ot_samse()], [ot_msy()]
#' @examples
#' # One row per stock. Columns:
#' #   abund = point abundance estimate (animals)
#' #   cv = coefficient of variation of that estimate
#' #   rmax = maximum annual net productivity rate
#' #   take = observed annual removal (animals)
#' d <- data.frame(
#'   stock = c("A", "B"),
#'   abund = c(1200, 800),
#'   cv = c(0.3, 0.2),
#'   rmax = c(0.04, 0.12),
#'   take = c(15, 40)
#' )
#' ot_pbr(d, rmax = rmax, removal = take, n = abund, cv = cv, fr = 0.5)
#' @export
ot_pbr <- function(data, rmax, removal, nmin = NULL, n = NULL, cv = NULL,
                   fr = 0.5, percentile = 0.20,
                   uncertainty = NULL, n_sim = 1000, level = 0.95, seed = NULL) {
  .check_data(data)
  groups <- .group_info(data)
  data <- .ungroup(data)
  Rmax <- .pull_col(data, rlang::enquo(rmax), "rmax")
  .check_range(Rmax, "rmax", 0, lower_open = TRUE)
  R <- .pull_col(data, rlang::enquo(removal), "removal")
  .check_range(R, "removal", 0)
  .check_range(percentile, "percentile", 0, 0.5, lower_open = TRUE)
  z <- stats::qnorm(percentile, lower.tail = FALSE)

  Nmin <- .pull_col(data, rlang::enquo(nmin), "nmin", required = FALSE)
  use_n <- is.null(Nmin)
  if (use_n) {
    N <- .pull_col(data, rlang::enquo(n), "n", required = FALSE)
    if (is.null(N)) {
      rlang::abort("Provide either `nmin`, or `n` (with optional `cv`) to compute it.",
                   class = "offtake_missing_arg")
    }
    .check_range(N, "n", 0, lower_open = TRUE)
    CV <- .pull_col(data, rlang::enquo(cv), "cv", required = FALSE)
    if (is.null(CV)) CV <- rep(0, length(N))
    .check_range(CV, "cv", 0)
    nmin_of <- function(N) N / exp(z * sqrt(log(1 + CV^2)))
    Nmin <- nmin_of(N)
  } else {
    .check_range(Nmin, "nmin", 0, lower_open = TRUE)
  }

  # `fr` may be a constant or an embraced column.
  FR <- rlang::eval_tidy(rlang::enquo(fr), data)
  FR <- if (length(FR) == 1L) rep(FR, nrow(data)) else FR
  .check_range(FR, "fr", 0, 1, lower_open = TRUE)
  if (any(FR < 0.1, na.rm = TRUE)) {
    rlang::warn("Some `fr` values are below 0.1, the lowest recovery factor considered by Wade (1998).",
                class = "offtake_unusual_value")
  }

  PBR <- Nmin * 0.5 * Rmax * FR
  res <- tibble::tibble(nmin = Nmin, rmax = Rmax, fr = FR)

  allowed <- c(if (use_n) "n" else "nmin", "rmax", "removal")
  unc <- .eval_uncertainty(rlang::enquo(uncertainty), data, allowed)
  if (is.null(unc)) {
    res <- .finish_model(res, PBR, R)
  } else {
    .check_range(n_sim, "n_sim", 1)
    .check_range(level, "level", 0, 1, lower_open = TRUE)
    if (!is.null(seed)) set.seed(seed)
    nmin_draws <- if (use_n) nmin_of(.draw_lnorm(N, unc$n, n_sim)) else .draw_lnorm(Nmin, unc$nmin, n_sim)
    rmax_draws <- .draw_lnorm(Rmax, unc$rmax, n_sim)
    r_draws <- .draw_lnorm(R, unc$removal, n_sim)
    pbr_draws <- nmin_draws * 0.5 * rmax_draws * FR
    res <- .finish_model(res, PBR, R, pbr_draws, r_draws, level)
  }
  if (!is.null(groups)) res <- cbind(.row_keys(groups, nrow(data)), res)

  new_offtake(
    res,
    method = "PBR (potential biological removal)",
    family = "model",
    reference = "Wade (1998)"
  )
}
