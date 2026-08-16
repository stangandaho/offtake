#' Sustainable anthropogenic mortality in stochastic environments (SAMSE)
#'
#' Estimates, by Monte-Carlo simulation, the largest constant annual removal
#' that does **not** drive a negative long-run stochastic growth rate, given
#' environmental variability in the population growth rate, and compares it with
#' the observed removal. SAMSE is a stochastic successor to [ot_pbr()] proposed
#' by Manlik et al. (2022) and highlighted by Adounke et al. (2026).
#'
#' @section Implementation note:
#' The original SAMSE limit of Manlik et al. (2022) is obtained with an
#' individual-based population viability analysis run in the *Vortex* software,
#' iterating removal levels until the forecast stochastic growth rate is no
#' longer negative. This function reproduces that **principle** with a
#' transparent, self-contained stochastic projection; it is **not** a
#' re-implementation of Vortex and does not include age structure, inbreeding or
#' catastrophes. Use it as a precautionary screening tool and cite Manlik et al.
#' (2022) for the concept.
#'
#' The projection is a stochastic exponential model with a constant annual take
#' `H`:
#' \deqn{N_{t+1} = \max\!\big(0,\; N_t\, e^{r_t} - H\big), \quad
#'       r_t \sim \mathrm{Normal}(r_{max},\, \sigma_e).}
#' The stochastic growth rate is \eqn{\rho(H) = \mathrm{mean}\,
#' \log(N_{t}/N_{0})/T}. The SAMSE limit is the largest `H` with
#' \eqn{\rho(H) \ge 0}, found by bisection. **An observed removal above the
#' SAMSE limit indicates unsustainability.**
#'
#' @param data A data frame with **one row per population/stock**. Required
#'   columns: mean maximum growth rate (`rmax`), its environmental standard
#'   deviation (`sd_env`), the observed annual removal (`removal`) and the
#'   starting population size (`n0`).
#' @param rmax <[`data-masked`][rlang::args_data_masking]> Mean maximum annual
#'   growth rate \eqn{r_{max}} on the log scale (so \eqn{\lambda = e^{r_{max}}};
#'   e.g. `0.10` for ~10% growth).
#' @param sd_env <[`data-masked`][rlang::args_data_masking]> Environmental
#'   standard deviation of the annual growth rate \eqn{\sigma_e} (larger =
#'   more year-to-year variability, lower sustainable removal).
#' @param removal <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   human-caused removal (animals per year).
#' @param n0 <[`data-masked`][rlang::args_data_masking]> Starting (ideally
#'   minimum, `N_min`-style) population size.
#' @param years Projection horizon in years (default `50`).
#' @param nsims Number of Monte-Carlo trajectories per candidate removal
#'   (default `500`).
#' @param tol Convergence tolerance of the bisection, as a fraction of `n0`
#'   (default `0.001`).
#' @param seed Optional integer seed for reproducibility (uses common random
#'   numbers across candidate removals).
#'
#' @return An [offtake][ot_samse] tibble with **one row per input row** and the
#'   columns:
#' \describe{
#'   \item{n0}{Starting population size used.}
#'   \item{rmax}{Mean maximum growth rate used.}
#'   \item{sd_env}{Environmental standard deviation used.}
#'   \item{samse_limit}{Estimated SAMSE limit -- the largest constant annual
#'     removal (animals per year) keeping the stochastic growth rate
#'     non-negative.}
#'   \item{removal}{The observed annual removal, echoed back for comparison.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `removal <= samse_limit`.}
#' }
#'
#' @references
#' Manlik, O., Lacy, R. C., Sherwin, W. B., Finn, H., Loneragan, N. R. & Allen,
#' S. J. (2022) A stochastic model for estimating sustainable limits to wildlife
#' mortality in a changing world. *Conservation Biology* 36, e13897.
#' \doi{10.1111/cobi.13897}
#'
#' @seealso [ot_pbr()]
#' @examples
#' # One row per stock. Columns:
#' #   rmax = mean max growth rate (log scale, per year)
#' #   sd = environmental SD of the annual growth rate
#' #   n0 = starting population size
#' #   take = observed annual removal (animals)
#' d <- data.frame(stock = "A", rmax = 0.10, sd = 0.25, n0 = 500, take = 20)
#' ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0,
#'          nsims = 200, years = 40, seed = 1)
#' @export
ot_samse <- function(data, rmax, sd_env, removal, n0,
                     years = 50, nsims = 500, tol = 0.001, seed = NULL) {
  .check_data(data)
  Rmax <- .pull_col(data, rlang::enquo(rmax), "rmax")
  Sd <- .pull_col(data, rlang::enquo(sd_env), "sd_env")
  R <- .pull_col(data, rlang::enquo(removal), "removal")
  N0 <- .pull_col(data, rlang::enquo(n0), "n0")

  limit <- vapply(seq_len(nrow(data)), function(i) {
    .samse_limit(Rmax[i], Sd[i], N0[i], years, nsims, tol, seed)
  }, numeric(1))

  res <- tibble::tibble(
    n0 = N0,
    rmax = Rmax,
    sd_env = Sd,
    samse_limit = limit,
    removal = R,
    sustainable = R <= limit
  )
  new_offtake(
    res,
    method = "SAMSE (sustainable anthropogenic mortality, stochastic)",
    family = "model",
    reference = "Manlik et al. (2022)",
    notes = "Simplified stochastic re-implementation of the SAMSE principle; not the original Vortex-based procedure."
  )
}

# Mean stochastic growth rate for a constant take H (common random numbers).
.samse_rho <- function(H, rmax, sd_env, n0, years, nsims, rand) {
  N <- rep(n0, nsims)
  alive <- rep(TRUE, nsims)
  for (t in seq_len(years)) {
    lam <- exp(rmax + sd_env * rand[, t])
    N <- pmax(0, N * lam - H)
    alive <- alive & N > 0
  }
  # Extinct trajectories contribute a strongly negative growth rate.
  ratio <- ifelse(N > 0, N / n0, .Machine$double.eps)
  mean(log(ratio)) / years
}

# Largest H with non-negative mean stochastic growth rate, by bisection.
.samse_limit <- function(rmax, sd_env, n0, years, nsims, tol, seed) {
  if (anyNA(c(rmax, sd_env, n0)) || n0 <= 0) return(NA_real_)
  if (!is.null(seed)) set.seed(seed)
  rand <- matrix(stats::rnorm(nsims * years), nrow = nsims, ncol = years)

  rho <- function(H) .samse_rho(H, rmax, sd_env, n0, years, nsims, rand)

  lo <- 0
  hi <- n0 * (exp(rmax) - 1)               # deterministic surplus at n0
  if (!is.finite(hi) || hi <= 0) return(0) # no surplus -> no sustainable take
  # Expand hi until rho(hi) < 0 (or give up).
  it <- 0
  while (rho(hi) >= 0 && it < 40) { hi <- hi * 1.5; it <- it + 1 }
  if (rho(hi) >= 0) return(hi)             # essentially unlimited within horizon
  # Bisection between lo (rho >= 0) and hi (rho < 0).
  eps <- tol * n0
  while (hi - lo > eps) {
    mid <- (lo + hi) / 2
    if (rho(mid) >= 0) lo <- mid else hi <- mid
  }
  lo
}
