# Internal helpers -------------------------------------------------------------

# Evaluate an embraced column argument against `data` and return a vector.
# `quo` is the result of rlang::enquo(); `arg` is used for error messages.
.pull_col <- function(data, quo, arg, required = TRUE) {
  if (rlang::quo_is_null(quo)) {
    if (required) {
      rlang::abort(
        sprintf("`%s` must name a column of `data`.", arg),
        class = "offtake_missing_arg"
      )
    }
    return(NULL)
  }
  out <- rlang::eval_tidy(quo, data)
  if (length(out) != nrow(data)) {
    rlang::abort(
      sprintf("`%s` must resolve to a column of length %d.", arg, nrow(data)),
      class = "offtake_bad_col"
    )
  }
  out
}

.check_data <- function(data) {
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame or tibble.", class = "offtake_bad_data")
  }
  invisible(data)
}

# Robinson & Redford (1991) mortality / longevity factor F.
# Longevity is the age (years) at which the species typically stops reproducing
# (age of last reproduction), used as a proxy for lifespan class.
.longevity_factor <- function(longevity) {
  vapply(longevity, function(x) {
    if (is.na(x)) return(NA_real_)
    if (x <= 5) return(0.6)      # short-lived, lifespan ~ 5 years
    if (x <= 10) return(0.4)     # medium-lived, 5-10 years
    0.2                          # long-lived, > 10 years
  }, numeric(1))
}

# Constructor for the classed result returned by every offtake function.
new_offtake <- function(x, method, family, reference, notes = NULL) {
  x <- tibble::as_tibble(x)
  structure(
    x,
    class = c("offtake", class(x)),
    method = method,
    family = family,
    reference = reference,
    notes = notes
  )
}

#' @export
print.offtake <- function(x, ...) {
  method <- attr(x, "method")
  family <- attr(x, "family")
  ref <- attr(x, "reference")
  notes <- attr(x, "notes")
  cat(sprintf("<offtake: %s>  (%s-based)\n", method, family))
  if (!is.null(ref)) cat("Reference:", ref, "\n")
  if (!is.null(notes)) cat("Note:", notes, "\n")
  cat("\n")
  # Print the underlying tibble without the classed attributes.
  NextMethod()
  invisible(x)
}

# Robinson & Redford / Cole solver --------------------------------------------

#' Maximum finite rate of population increase (Cole's equation)
#'
#' Solves Cole's (1954) characteristic equation for the maximum intrinsic rate
#' of natural increase, returning the finite (per-year) multiplication rate
#' \eqn{\lambda_{max} = e^{r_{max}}}. This is the rate required by the
#' production model [ot_pro()] when it is not supplied directly.
#'
#' Cole's equation is
#' \deqn{1 = e^{-r} + b\,e^{-r a} - b\,e^{-r (w + 1)}}
#' where `b` is the mean annual number of *female* offspring per female, `a`
#' is the age at first reproduction and `w` is the age at last reproduction
#' (all in years). The equation is solved numerically for \eqn{r} and
#' \eqn{\lambda_{max}=e^{r}} is returned.
#'
#' @param b Mean annual number of female offspring per female (i.e.
#'   litter size \eqn{\times} litters per year \eqn{\times} proportion female).
#' @param a Age at first reproduction (years).
#' @param w Age at last reproduction (years); often approximated by lifespan.
#' @param interval Search interval for \eqn{r_{max}} passed to
#'   [stats::uniroot()].
#'
#' @return A numeric vector, recycled to the length of the longest of `b`, `a`
#'   and `w`, giving \eqn{\lambda_{max}} -- the maximum finite (per-year) rate
#'   of increase. It is dimensionless and greater than 1 for a growing
#'   population (e.g. `1.30` means the population can grow by at most 30% per
#'   year). `NA` is returned for elements with missing inputs or with no
#'   positive root inside `interval`.
#'
#' @references
#' Cole, L. C. (1954) The population consequences of life history phenomena.
#' *The Quarterly Review of Biology* 29, 103-137.
#'
#' Robinson, J. G. & Redford, K. H. (1986) Intrinsic rate of natural increase
#' in Neotropical forest mammals. *Oecologia* 68, 516-520.
#'
#' @examples
#' # Inputs are life-history traits (not columns of a data frame here):
#' #   b = female young per female per year, a = age at first reproduction,
#' #   w = age at last reproduction. A duiker breeding from age 1 to 8:
#' ot_lambda_max(b = 0.6, a = 1, w = 8)
#'
#' # Vectorised over several species at once:
#' ot_lambda_max(b = c(0.6, 0.5), a = c(1, 2), w = c(8, 12))
#' @export
ot_lambda_max <- function(b, a, w, interval = c(1e-6, 5)) {
  n <- max(length(b), length(a), length(w))
  b <- rep_len(b, n); a <- rep_len(a, n); w <- rep_len(w, n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    if (anyNA(c(b[i], a[i], w[i]))) { out[i] <- NA_real_; next }
    f <- function(r) exp(-r) + b[i] * exp(-r * a[i]) -
      b[i] * exp(-r * (w[i] + 1)) - 1
    # f is decreasing in r; require a sign change over `interval`.
    lo <- f(interval[1]); hi <- f(interval[2])
    if (is.na(lo) || is.na(hi) || lo * hi > 0) { out[i] <- NA_real_; next }
    r <- stats::uniroot(f, interval = interval)$root
    out[i] <- exp(r)
  }
  out
}
