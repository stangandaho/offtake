#' Maximum sustainable yield of the logistic model (MSY)
#'
#' Estimates the maximum sustainable yield from the logistic (surplus
#' production) model and compares it with the observed annual harvest. This is
#' the stock-recruitment / surplus-production benchmark used across fisheries
#' and, increasingly, wildlife harvest (Adounke et al. 2026).
#'
#' The logistic model of population growth is
#' \deqn{\frac{dN}{dt} = r N \left(1 - \frac{N}{K}\right)}
#' whose surplus production is maximised at \eqn{N = K/2}, giving
#' \deqn{MSY = \frac{rK}{4}.}
#' **A harvest above MSY is considered unsustainable.** When a current
#' abundance `n` is supplied, the instantaneous surplus production at that
#' abundance, \eqn{rN(1 - N/K)}, is also returned as `surplus_at_n`, which is
#' the sustainable yield at the *current* (not optimal) population size. A
#' population already below \eqn{K/2} cannot sustain MSY.
#'
#' Grouped data frames (from `dplyr::group_by()`) are accepted; the group
#' columns are kept in the output.
#'
#' @inheritSection ot_pro Uncertainty
#' @inheritParams ot_pro
#' @param data A data frame with **one row per population/stock**. Required
#'   columns: intrinsic growth rate (`r`), carrying capacity (`k`) and observed
#'   annual harvest (`harvest`); optionally current abundance (`n`). May be
#'   grouped with `dplyr::group_by()`.
#' @param r <[`data-masked`][rlang::args_data_masking]> Intrinsic rate of
#'   natural increase (per year), strictly positive.
#' @param k <[`data-masked`][rlang::args_data_masking]> Carrying capacity `K`
#'   (population size or density), strictly positive.
#' @param harvest <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   harvest, on the same basis as `k`.
#' @param n <[`data-masked`][rlang::args_data_masking]> Optional current
#'   abundance, used to compute `surplus_at_n`.
#' @param uncertainty Optional named list of coefficients of variation (numbers
#'   or columns of `data`). Allowed names: `r`, `k` and `harvest`.
#'
#' @return An [offtake][ot_msy] tibble with one row per input row and the
#'   columns:
#' \describe{
#'   \item{r}{Intrinsic growth rate used.}
#'   \item{k}{Carrying capacity used.}
#'   \item{limit}{Maximum sustainable yield, `(r * k) / 4`.}
#'   \item{observed}{The observed annual harvest (`harvest`).}
#'   \item{ratio}{Exploitation ratio, `observed / limit`.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `observed <= limit`.}
#'   \item{limit_lo, limit_hi, p_unsustainable}{Only with `uncertainty`; see
#'     the Uncertainty section.}
#'   \item{surplus_at_n}{Only present when `n` is supplied: the sustainable
#'     yield at the *current* abundance, `(r * n) * (1 - n / k)`. Equals `limit`
#'     when `n = k / 2`.}
#' }
#'   Group columns come first when `data` is grouped.
#'
#' @references
#' Schaefer, M. B. (1954) Some aspects of the dynamics of populations important
#' to the management of the commercial marine fisheries. *Bulletin of the
#' Inter-American Tropical Tuna Commission* 1, 27-56.
#'
#' Adounke, G. R. M. et al. (2026) Systematic review of sustainability
#' assessment approaches for wildlife exploitation. *Biological Conservation*
#' 313, 111606. \doi{10.1016/j.biocon.2025.111606}
#'
#' @seealso [ot_pro()], [ot_pbr()]
#' @examples
#' # One row per stock. Columns:
#' #   r = intrinsic growth rate (per year)
#' #   K = carrying capacity
#' #   take = observed annual harvest
#' #   now = current abundance (optional, for surplus_at_n)
#' d <- data.frame(stock = "A", r = 0.4, K = 1000, take = 80, now = 600)
#' ot_msy(d, r = r, k = K, harvest = take, n = now)
#' @export
ot_msy <- function(data, r, k, harvest, n = NULL,
                   uncertainty = NULL, n_sim = 1000, level = 0.95, seed = NULL) {
  .check_data(data)
  groups <- .group_info(data)
  data <- .ungroup(data)
  rr <- .pull_col(data, rlang::enquo(r), "r")
  .check_range(rr, "r", 0, lower_open = TRUE)
  K <- .pull_col(data, rlang::enquo(k), "k")
  .check_range(K, "k", 0, lower_open = TRUE)
  H <- .pull_col(data, rlang::enquo(harvest), "harvest")
  .check_range(H, "harvest", 0)
  N <- .pull_col(data, rlang::enquo(n), "n", required = FALSE)
  if (!is.null(N)) .check_range(N, "n", 0)

  MSY <- rr * K / 4
  res <- tibble::tibble(r = rr, k = K)

  unc <- .eval_uncertainty(rlang::enquo(uncertainty), data, c("r", "k", "harvest"))
  if (is.null(unc)) {
    res <- .finish_model(res, MSY, H)
  } else {
    .check_range(n_sim, "n_sim", 1)
    .check_range(level, "level", 0, 1, lower_open = TRUE)
    if (!is.null(seed)) set.seed(seed)
    msy_draws <- .draw_lnorm(rr, unc$r, n_sim) * .draw_lnorm(K, unc$k, n_sim) / 4
    h_draws <- .draw_lnorm(H, unc$harvest, n_sim)
    res <- .finish_model(res, MSY, H, msy_draws, h_draws, level)
  }
  if (!is.null(N)) res$surplus_at_n <- rr * N * (1 - N / K)
  if (!is.null(groups)) res <- cbind(.row_keys(groups, nrow(data)), res)

  new_offtake(
    res,
    method = "MSY (logistic maximum sustainable yield)",
    family = "model",
    reference = "Schaefer (1954)"
  )
}
