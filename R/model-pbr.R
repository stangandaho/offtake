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
#' \eqn{R_{max}} is the maximum per-capita annual rate of increase and
#' \eqn{F_R} is a recovery factor between 0.1 and 1. **A removal exceeding PBR
#' indicates unsustainability.**
#'
#' `N_min` may be supplied directly through `nmin`, or computed from an
#' abundance estimate `n` and its coefficient of variation `cv` as the lower
#' tail of a log-normal distribution,
#' \eqn{N_{min} = N / \exp\!\big(z\sqrt{\ln(1+CV^2)}\big)}, with `z` the
#' standard-normal quantile for the chosen percentile (`z = 0.842` for the 20th
#' percentile, as recommended by Wade 1998).
#'
#' @param data A data frame with **one row per population/stock**. Required
#'   columns: the maximum growth rate (`rmax`) and the observed annual removal
#'   (`removal`). You must also supply the population size, either as a minimum
#'   estimate (`nmin`) or as a point estimate (`n`) with an optional coefficient
#'   of variation (`cv`) from which `nmin` is derived.
#' @param rmax <[`data-masked`][rlang::args_data_masking]> Maximum per-capita
#'   annual rate of increase \eqn{R_{max}} (e.g. `0.04` for many large mammals,
#'   `0.12` for fast-breeding species).
#' @param removal <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   human-caused removal (offtake), in numbers of animals.
#' @param nmin <[`data-masked`][rlang::args_data_masking]> Minimum population
#'   estimate. Optional if `n` (and optionally `cv`) are supplied.
#' @param n <[`data-masked`][rlang::args_data_masking]> Point abundance estimate,
#'   used with `cv` to derive `nmin`.
#' @param cv <[`data-masked`][rlang::args_data_masking]> Coefficient of variation
#'   of `n` (e.g. `0.3` for a 30% CV). Treated as 0 (so `nmin = n`) when omitted.
#' @param fr Recovery factor \eqn{F_R} in `[0.1, 1]` (default `0.5`). May be a
#'   single value or an embraced column.
#' @param percentile Lower percentile used to convert `n`/`cv` to `nmin`
#'   (default `0.20`, giving `z = 0.842`).
#'
#' @return An [offtake][ot_pbr] tibble with **one row per input row** and the
#'   columns:
#' \describe{
#'   \item{nmin}{Minimum population estimate used (supplied, or derived from
#'     `n` and `cv`).}
#'   \item{rmax}{Maximum per-capita rate of increase used.}
#'   \item{fr}{Recovery factor applied.}
#'   \item{pbr}{Potential biological removal -- the maximum sustainable annual
#'     removal, in numbers of animals.}
#'   \item{removal}{The observed annual removal, echoed back for comparison.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `removal <= pbr`, `FALSE`
#'     when the observed removal exceeds PBR.}
#' }
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
#' #   rmax = maximum per-capita rate of increase (per year)
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
                   fr = 0.5, percentile = 0.20) {
  .check_data(data)
  Rmax <- .pull_col(data, rlang::enquo(rmax), "rmax")
  R <- .pull_col(data, rlang::enquo(removal), "removal")

  Nmin <- .pull_col(data, rlang::enquo(nmin), "nmin", required = FALSE)
  if (is.null(Nmin)) {
    N <- .pull_col(data, rlang::enquo(n), "n", required = FALSE)
    if (is.null(N)) {
      rlang::abort("Provide either `nmin`, or `n` (with optional `cv`) to compute it.",
                   class = "offtake_missing_arg")
    }
    CV <- .pull_col(data, rlang::enquo(cv), "cv", required = FALSE)
    if (is.null(CV)) CV <- rep(0, length(N))
    z <- stats::qnorm(percentile, lower.tail = FALSE)
    Nmin <- N / exp(z * sqrt(log(1 + CV^2)))
  }

  # `fr` may be a constant or an embraced column.
  frq <- rlang::enquo(fr)
  FR <- if (rlang::quo_is_null(frq)) rep(0.5, nrow(data)) else {
    val <- rlang::eval_tidy(frq, data)
    if (length(val) == 1L) rep(val, nrow(data)) else val
  }

  PBR <- Nmin * 0.5 * Rmax * FR
  res <- tibble::tibble(
    nmin = Nmin,
    rmax = Rmax,
    fr = FR,
    pbr = PBR,
    removal = R,
    sustainable = R <= PBR
  )
  new_offtake(
    res,
    method = "PBR (potential biological removal)",
    family = "model",
    reference = "Wade (1998)"
  )
}
