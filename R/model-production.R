#' Robinson & Redford production model (Pro)
#'
#' Estimates the maximum sustainable harvest (production) of a population from
#' its carrying capacity, maximum rate of increase and a longevity-based
#' mortality factor, and compares it with the observed annual offtake. This is
#' the most widely used model in bushmeat sustainability assessments
#' (Adounke et al. 2026).
#'
#' The production (maximum sustainable number that can be taken per year) is
#' \deqn{P = 0.6\,K\,(\lambda_{max} - 1)\,F}
#' where `K` is carrying-capacity density (or population size), \eqn{\lambda_{max}}
#' is the maximum finite rate of increase and `F` is the mortality factor set by
#' longevity: `F = 0.6` for short-lived species (lifespan ~ 5 years), `0.4` for
#' medium-lived (5-10 years) and `0.2` for long-lived species (> 10 years)
#' (Robinson & Redford 1991). The factor `0.6 K` is the density at which
#' production is assumed maximal. **If the observed harvest exceeds `P` the
#' harvest is considered unsustainable.**
#'
#' Supply \eqn{\lambda_{max}} through `lambda` directly, or leave it `NULL` and
#' provide the life-history columns `b`, `a`, `w` so it is computed with
#' [ot_lambda_max()] (Cole's equation). Provide the mortality factor either
#' through `f` directly or through a `longevity` column.
#'
#' Grouped data frames (from `dplyr::group_by()`) are accepted; the group
#' columns are kept in the output.
#'
#' @section Uncertainty:
#' Model inputs are rarely known precisely: \eqn{\lambda_{max}} from Cole's
#' equation and from the Caughley & Krebs (1983) equation can differ by up to
#' four times for the same species, and densities depend on the survey method
#' (Adounke et al. 2026). With `uncertainty`, each named input is drawn
#' `n_sim` times from a lognormal distribution with the given value as its mean
#' and the given coefficient of variation (CV), the safe limit is recomputed for
#' every draw, and the output gains three columns: `limit_lo` and `limit_hi`
#' (the central `level` interval of the limit) and `p_unsustainable` (the share
#' of draws in which the observed offtake exceeds the limit). This turns the
#' TRUE/FALSE verdict into a probability. For `lambda`, the CV applies to the
#' annual surplus \eqn{\lambda_{max} - 1}, so that drawn growth rates stay above
#' 1.
#'
#' @param data A data frame with **one row per species or population unit**.
#'   Required columns: carrying capacity (`k`) and observed annual offtake
#'   (`harvest`), on the *same basis* (both per km^2, or both absolute counts).
#'   You must also supply the growth rate -- either a `lambda` column or the
#'   life-history columns `b`, `a`, `w` -- and the mortality factor -- either an
#'   `f` column or a `longevity` column. May be grouped with
#'   `dplyr::group_by()`.
#' @param k <[`data-masked`][rlang::args_data_masking]> Carrying capacity `K`
#'   (density per km^2, or absolute population size), strictly positive.
#' @param harvest <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   offtake, on the same basis as `k`.
#' @param lambda <[`data-masked`][rlang::args_data_masking]> Maximum finite rate
#'   of increase \eqn{\lambda_{max}} (dimensionless, > 1). Optional if `b`, `a`,
#'   `w` are given.
#' @param f <[`data-masked`][rlang::args_data_masking]> Mortality factor `F`
#'   (0.2, 0.4 or 0.6). Optional if `longevity` is given.
#' @param longevity <[`data-masked`][rlang::args_data_masking]> Age at last
#'   reproduction / lifespan (years), used to derive `F` when `f` is not given.
#' @param b,a,w <[`data-masked`][rlang::args_data_masking]> Life-history columns
#'   passed to [ot_lambda_max()] when `lambda` is not given: annual female
#'   offspring per female (`b`), age at first reproduction (`a`) and age at last
#'   reproduction (`w`), all in years.
#' @param uncertainty Optional named list of coefficients of variation, one per
#'   uncertain input; each value can be a number or a column of `data`. Here the
#'   allowed names are `k`, `lambda` and `harvest`, e.g.
#'   `list(k = 0.3, lambda = 0.2)`. See the Uncertainty section.
#' @param n_sim Number of draws used when `uncertainty` is given (default
#'   `1000`).
#' @param level Width of the interval `limit_lo` to `limit_hi` (default
#'   `0.95`).
#' @param seed Optional integer seed, for reproducible draws.
#'
#' @return An [offtake][ot_pro] tibble with **one row per input row** and the
#'   columns:
#' \describe{
#'   \item{lambda_max}{Maximum finite rate of increase used, whether supplied or
#'     computed from `b`, `a`, `w`.}
#'   \item{f}{Mortality factor `F` used (0.2, 0.4 or 0.6).}
#'   \item{limit}{Estimated production `P`, the maximum sustainable harvest, on
#'     the same basis as `k` and `harvest` (e.g. individuals per km^2 per
#'     year).}
#'   \item{observed}{The observed annual offtake (`harvest`).}
#'   \item{ratio}{Exploitation ratio, `observed / limit`. Above 1 means the
#'     offtake exceeds the production; 3 means three times too much.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `observed <= limit`.}
#'   \item{limit_lo, limit_hi, p_unsustainable}{Only with `uncertainty`: the
#'     interval of the limit and the probability that the offtake exceeds it.}
#' }
#'   Group columns come first when `data` is grouped.
#'
#' @references
#' Robinson, J. G. & Redford, K. H. (1991) Sustainable harvest of Neotropical
#' forest mammals. In *Neotropical Wildlife Use and Conservation* (eds Robinson
#' & Redford), 415-429. University of Chicago Press.
#'
#' Adounke, G. R. M. et al. (2026) Systematic review of sustainability
#' assessment approaches for wildlife exploitation. *Biological Conservation*
#' 313, 111606. \doi{10.1016/j.biocon.2025.111606}
#'
#' @seealso [ot_pbr()], [ot_msy()], [ot_lambda_max()]
#' @examples
#' # One row per species. Columns:
#' #   dens = carrying-capacity density (ind / km^2)
#' #   offtake = observed annual harvest (ind / km^2 / year)
#' #   lifespan = age at last reproduction (years) -> sets F
#' #   lam = maximum finite rate of increase (lambda_max)
#' d <- data.frame(
#'   species = c("red duiker", "blue duiker"),
#'   dens = c(10, 25),
#'   offtake = c(3.0, 6.0),
#'   lifespan = c(9, 6),
#'   lam = c(1.35, 1.55)
#' )
#' ot_pro(d, k = dens, harvest = offtake, lambda = lam, longevity = lifespan)
#'
#' # With 30% uncertainty on K and 20% on the growth surplus
#' ot_pro(d, k = dens, harvest = offtake, lambda = lam, longevity = lifespan,
#'        uncertainty = list(k = 0.3, lambda = 0.2), seed = 1)
#' @export
ot_pro <- function(data, k, harvest, lambda = NULL, f = NULL, longevity = NULL,
                   b = NULL, a = NULL, w = NULL,
                   uncertainty = NULL, n_sim = 1000, level = 0.95, seed = NULL) {
  .check_data(data)
  groups <- .group_info(data)
  data <- .ungroup(data)
  K <- .pull_col(data, rlang::enquo(k), "k")
  .check_range(K, "k", 0, lower_open = TRUE)
  H <- .pull_col(data, rlang::enquo(harvest), "harvest")
  .check_range(H, "harvest", 0)

  lam <- .pull_col(data, rlang::enquo(lambda), "lambda", required = FALSE)
  if (is.null(lam)) {
    bb <- .pull_col(data, rlang::enquo(b), "b", required = FALSE)
    aa <- .pull_col(data, rlang::enquo(a), "a", required = FALSE)
    ww <- .pull_col(data, rlang::enquo(w), "w", required = FALSE)
    if (is.null(bb) || is.null(aa) || is.null(ww)) {
      rlang::abort("Provide either `lambda`, or all of `b`, `a` and `w` to compute it.",
                   class = "offtake_missing_arg")
    }
    lam <- ot_lambda_max(bb, aa, ww)
  }
  .check_range(lam, "lambda", 0, lower_open = TRUE)
  if (any(lam <= 1, na.rm = TRUE)) {
    rlang::warn("Some `lambda` values are <= 1: the population has no surplus, so the production is zero or negative.",
                class = "offtake_no_surplus")
  }

  ff <- .pull_col(data, rlang::enquo(f), "f", required = FALSE)
  if (is.null(ff)) {
    lon <- .pull_col(data, rlang::enquo(longevity), "longevity", required = FALSE)
    if (is.null(lon)) {
      rlang::abort("Provide either `f` (mortality factor) or `longevity`.",
                   class = "offtake_missing_arg")
    }
    .check_range(lon, "longevity", 0, lower_open = TRUE)
    ff <- .longevity_factor(lon)
  }
  .check_range(ff, "f", 0, 1, lower_open = TRUE)

  P <- 0.6 * K * (lam - 1) * ff
  res <- tibble::tibble(lambda_max = lam, f = ff)

  unc <- .eval_uncertainty(rlang::enquo(uncertainty), data, c("k", "lambda", "harvest"))
  if (is.null(unc)) {
    res <- .finish_model(res, P, H)
  } else {
    .check_range(n_sim, "n_sim", 1)
    .check_range(level, "level", 0, 1, lower_open = TRUE)
    if (!is.null(seed)) set.seed(seed)
    k_draws <- .draw_lnorm(K, unc$k, n_sim)
    surplus_draws <- .draw_lnorm(lam - 1, unc$lambda, n_sim)
    h_draws <- .draw_lnorm(H, unc$harvest, n_sim)
    p_draws <- 0.6 * k_draws * surplus_draws * ff
    res <- .finish_model(res, P, H, p_draws, h_draws, level)
  }
  if (!is.null(groups)) res <- cbind(.row_keys(groups, nrow(data)), res)

  new_offtake(
    res,
    method = "Pro (Robinson & Redford production model)",
    family = "model",
    reference = "Robinson & Redford (1991)"
  )
}
