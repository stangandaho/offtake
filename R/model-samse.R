#' Sustainable anthropogenic mortality in stochastic environments (SAMSE)
#'
#' Estimates, by Monte-Carlo simulation, the largest constant annual removal
#' that does **not** produce a negative stochastic growth rate from the current
#' population size, given environmental variability in the growth rate and
#' density dependence, and compares it with the observed removal. SAMSE is a
#' stochastic successor to [ot_pbr()] proposed by Manlik et al. (2022) and
#' highlighted by Adounke et al. (2026).
#'
#' @section Model:
#' The population is projected with a theta-logistic (Ricker type) model with
#' environmental noise and a constant annual take `H`:
#' \deqn{N_{t+1} = \max\!\Big(0,\; N_t \exp\!\big[r_{max}\big(1 - (N_t/K)^{\theta}\big)
#'       + \varepsilon_t\big] - H\Big), \quad \varepsilon_t \sim \mathrm{Normal}(0, \sigma_e).}
#' Growth slows as the population approaches its carrying capacity `K`, so a
#' population near `K` has little surplus to give. For each candidate `H`,
#' `nsims` trajectories of `years` years are run from `n0`. The stochastic
#' growth rate is the mean of the annual \eqn{\log(N_{t+1}/N_t)} over all
#' trajectories and years in which the population is still present, and the
#' extinction probability is the share of trajectories that reach zero. The
#' SAMSE limit is the largest `H` for which the stochastic growth rate is not
#' negative **and** the extinction probability does not exceed
#' `max_extinction`, found by bisection with common random numbers. **An
#' observed removal above the SAMSE limit indicates unsustainability.**
#'
#' Because the criterion is "no decline from the current size", the limit
#' depends on `n0`: it is close to the surplus the population produces at that
#' size, reduced by environmental variability. Use the current abundance for
#' `n0`; a population at its carrying capacity has (almost) no surplus.
#'
#' @section Implementation note:
#' Manlik et al. (2022) obtain the SAMSE limit with an individual-based
#' population viability analysis in the *Vortex* software. This function
#' reproduces the principle (largest removal without a negative stochastic
#' growth rate) with a transparent count-based projection. It does not include
#' age structure, demographic stochasticity, inbreeding or catastrophes, and the
#' extinction constraint (`max_extinction`) is an addition of this package.
#'
#' Grouped data frames (from `dplyr::group_by()`) are accepted; the group
#' columns are kept in the output.
#'
#' @param data A data frame with **one row per population/stock**. Required
#'   columns: mean maximum growth rate (`rmax`), its environmental standard
#'   deviation (`sd_env`), the observed annual removal (`removal`), the current
#'   population size (`n0`) and the carrying capacity (`k`). May be grouped
#'   with `dplyr::group_by()`.
#' @param rmax <[`data-masked`][rlang::args_data_masking]> Maximum annual growth
#'   rate \eqn{r_{max}} on the log scale (growth of a small population is
#'   about \eqn{e^{r_{max}}} per year; e.g. `0.10` for ~10% growth).
#' @param sd_env <[`data-masked`][rlang::args_data_masking]> Environmental
#'   standard deviation of the annual growth rate \eqn{\sigma_e} (larger =
#'   more year-to-year variability, lower sustainable removal).
#' @param removal <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   human-caused removal, on the same basis as `n0` and `k`.
#' @param n0 <[`data-masked`][rlang::args_data_masking]> Current population size
#'   (or density), the starting point of the projections.
#' @param k <[`data-masked`][rlang::args_data_masking]> Carrying capacity, on the
#'   same basis as `n0`.
#' @param theta Shape of density dependence \eqn{\theta} (default `1`,
#'   logistic; values above 1 keep growth high until the population is close
#'   to `K`).
#' @param years Projection horizon in years (default `50`).
#' @param nsims Number of Monte-Carlo trajectories per candidate removal
#'   (default `500`).
#' @param max_extinction Largest acceptable probability that a trajectory goes
#'   extinct within `years` (default `0.05`).
#' @param tol Convergence tolerance of the bisection, as a fraction of `n0`
#'   (default `0.001`).
#' @param seed Optional integer seed for reproducibility (common random numbers
#'   are used across candidate removals).
#'
#' @return An [offtake][ot_samse] tibble with one row per input row and the
#'   columns:
#' \describe{
#'   \item{n0}{Current population size used.}
#'   \item{k}{Carrying capacity used.}
#'   \item{rmax}{Maximum growth rate used.}
#'   \item{sd_env}{Environmental standard deviation used.}
#'   \item{limit}{Estimated SAMSE limit: the largest constant annual removal
#'     keeping the stochastic growth rate non-negative and the extinction
#'     probability at or below `max_extinction`.}
#'   \item{observed}{The observed annual removal (`removal`).}
#'   \item{ratio}{Exploitation ratio, `observed / limit`.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `observed <= limit`.}
#'   \item{p_extinct}{Probability of extinction within `years` if the observed
#'     removal continues.}
#' }
#'   Group columns come first when `data` is grouped.
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
#' #   rmax = maximum growth rate (log scale, per year)
#' #   sd = environmental SD of the annual growth rate
#' #   n0 = current population size
#' #   K = carrying capacity
#' #   take = observed annual removal (animals)
#' d <- data.frame(stock = "A", rmax = 0.10, sd = 0.25, n0 = 500, K = 1000, take = 20)
#' ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
#'          nsims = 200, years = 40, seed = 1)
#' @export
ot_samse <- function(data, rmax, sd_env, removal, n0, k, theta = 1,
                     years = 50, nsims = 500, max_extinction = 0.05,
                     tol = 0.001, seed = NULL) {
  .check_data(data)
  groups <- .group_info(data)
  data <- .ungroup(data)
  Rmax <- .pull_col(data, rlang::enquo(rmax), "rmax")
  .check_range(Rmax, "rmax", 0, lower_open = TRUE)
  Sd <- .pull_col(data, rlang::enquo(sd_env), "sd_env")
  .check_range(Sd, "sd_env", 0)
  R <- .pull_col(data, rlang::enquo(removal), "removal")
  .check_range(R, "removal", 0)
  N0 <- .pull_col(data, rlang::enquo(n0), "n0")
  .check_range(N0, "n0", 0, lower_open = TRUE)
  K <- .pull_col(data, rlang::enquo(k), "k")
  .check_range(K, "k", 0, lower_open = TRUE)
  .check_range(theta, "theta", 0, lower_open = TRUE)
  .check_range(years, "years", 1)
  .check_range(nsims, "nsims", 1)
  .check_range(max_extinction, "max_extinction", 0, 1)
  .check_range(tol, "tol", 0, lower_open = TRUE)
  if (any(N0 > K, na.rm = TRUE)) {
    rlang::warn("Some `n0` values are above `k`: those populations decline even without removals.",
                class = "offtake_unusual_value")
  }

  out <- lapply(seq_len(nrow(data)), function(i) {
    .samse_one(Rmax[i], Sd[i], N0[i], K[i], theta, R[i], years, nsims,
               max_extinction, tol, seed)
  })
  limit <- vapply(out, `[[`, numeric(1), "limit")
  p_ext <- vapply(out, `[[`, numeric(1), "p_extinct")

  res <- tibble::tibble(n0 = N0, k = K, rmax = Rmax, sd_env = Sd)
  res <- .finish_model(res, limit, R)
  res$p_extinct <- p_ext
  if (!is.null(groups)) res <- cbind(.row_keys(groups, nrow(data)), res)

  new_offtake(
    res,
    method = "SAMSE (sustainable anthropogenic mortality, stochastic)",
    family = "model",
    reference = "Manlik et al. (2022)",
    notes = "Count-based re-implementation of the SAMSE principle with density dependence; not the original Vortex-based procedure."
  )
}

# Project `nsims` trajectories with a constant take H. Returns the stochastic
# growth rate (mean annual log growth while the population is present) and the
# share of trajectories that went extinct.
.samse_sim <- function(H, rmax, sd_env, n0, k, theta, rand) {
  N <- rep(n0, nrow(rand))
  sum_log <- 0
  n_steps <- 0
  for (t in seq_len(ncol(rand))) {
    alive <- N > 0
    if (!any(alive)) break
    Na <- N[alive]
    growth <- rmax * (1 - (Na / k)^theta) + sd_env * rand[alive, t]
    nxt <- pmax(0, Na * exp(growth) - H)
    present <- nxt > 0
    sum_log <- sum_log + sum(log(nxt[present] / Na[present]))
    n_steps <- n_steps + sum(present)
    N[alive] <- nxt
  }
  list(r = if (n_steps > 0) sum_log / n_steps else -Inf,
       p_extinct = mean(N <= 0))
}

.samse_one <- function(rmax, sd_env, n0, k, theta, removal, years, nsims,
                       max_extinction, tol, seed) {
  if (anyNA(c(rmax, sd_env, n0, k))) return(list(limit = NA_real_, p_extinct = NA_real_))
  if (!is.null(seed)) set.seed(seed)
  rand <- matrix(stats::rnorm(nsims * years), nrow = nsims, ncol = years)
  ok <- function(H) {
    s <- .samse_sim(H, rmax, sd_env, n0, k, theta, rand)
    s$r >= 0 && s$p_extinct <= max_extinction
  }
  p_ext <- if (is.na(removal)) NA_real_ else
    .samse_sim(removal, rmax, sd_env, n0, k, theta, rand)$p_extinct

  if (!ok(0)) return(list(limit = 0, p_extinct = p_ext))
  lo <- 0
  hi <- max(rmax * k, n0)
  it <- 0
  while (ok(hi) && it < 40) { lo <- hi; hi <- hi * 2; it <- it + 1 }
  eps <- tol * n0
  while (hi - lo > eps) {
    mid <- (lo + hi) / 2
    if (ok(mid)) lo <- mid else hi <- mid
  }
  list(limit = lo, p_extinct = p_ext)
}
