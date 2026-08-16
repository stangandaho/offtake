#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang enquo eval_tidy quo_is_null abort warn %||%
#' @importFrom stats prop.test t.test uniroot chisq.test rnorm quantile sd qnorm
#' @importFrom tibble tibble as_tibble
## usethis namespace: end
NULL

# Quiet R CMD check notes about tidy-evaluation pronouns / NSE columns.
utils::globalVariables(c(".data"))
