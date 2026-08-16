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
#' the sustainable yield at the *current* (not optimal) population size.
#'
#' @param data A data frame with **one row per population/stock**. Required
#'   columns: intrinsic growth rate (`r`), carrying capacity (`k`) and observed
#'   annual harvest (`harvest`); optionally current abundance (`n`).
#' @param r <[`data-masked`][rlang::args_data_masking]> Intrinsic rate of
#'   natural increase (per year).
#' @param k <[`data-masked`][rlang::args_data_masking]> Carrying capacity `K`
#'   (population size or density).
#' @param harvest <[`data-masked`][rlang::args_data_masking]> Observed annual
#'   harvest, on the same basis as `k`.
#' @param n <[`data-masked`][rlang::args_data_masking]> Optional current
#'   abundance, used to compute `surplus_at_n`.
#'
#' @return An [offtake][ot_msy] tibble with **one row per input row** and the
#'   columns:
#' \describe{
#'   \item{r}{Intrinsic growth rate used.}
#'   \item{k}{Carrying capacity used.}
#'   \item{msy}{Maximum sustainable yield, `(r * k) / 4`.}
#'   \item{harvest}{The observed annual harvest, echoed back for comparison.}
#'   \item{sustainable}{Logical verdict: `TRUE` when `harvest <= msy`.}
#'   \item{surplus_at_n}{Only present when `n` is supplied: the sustainable
#'     yield at the *current* abundance, `(r * n) * (1 - n / k)`. Equals `msy`
#'     when `n = k / 2`.}
#' }
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
ot_msy <- function(data, r, k, harvest, n = NULL) {
  .check_data(data)
  rr <- .pull_col(data, rlang::enquo(r), "r")
  K <- .pull_col(data, rlang::enquo(k), "k")
  H <- .pull_col(data, rlang::enquo(harvest), "harvest")
  MSY <- rr * K / 4

  res <- tibble::tibble(
    r = rr,
    k = K,
    msy = MSY,
    harvest = H,
    sustainable = H <= MSY
  )
  N <- .pull_col(data, rlang::enquo(n), "n", required = FALSE)
  if (!is.null(N)) {
    res$surplus_at_n <- rr * N * (1 - N / K)
  }
  new_offtake(
    res,
    method = "MSY (logistic maximum sustainable yield)",
    family = "model",
    reference = "Schaefer (1954)"
  )
}
