# Shared comparison engine for the index-based approaches ----------------------
#
# All three indices catalogued by Adounke et al. (2026) that are implemented
# here (ASC, HYCo, PDC) share the same logic: a metric measured at one or more
# *hunted* sites is compared with the same metric at a *reference*
# (unhunted / lightly hunted) site, and a lower value at the hunted site is
# interpreted as unsustainable (Weinbaum et al. 2013, Table 1).

# Compare a continuous metric between hunted and reference groups.
# `lower_is_unsustainable = TRUE` means: hunted significantly below reference
# -> flagged unsustainable.
.compare_groups <- function(value, group, reference, alpha,
                            lower_is_unsustainable = TRUE) {
  group <- as.character(group)
  levs <- unique(group)
  if (!reference %in% levs) {
    rlang::abort(
      sprintf("`reference` = \"%s\" is not a level of the grouping column (levels: %s).",
              reference, paste(levs, collapse = ", ")),
      class = "offtake_bad_reference"
    )
  }
  ref_val <- value[group == reference]
  hunted_levs <- setdiff(levs, reference)
  if (length(hunted_levs) == 0L) {
    rlang::abort("No hunted group found: the grouping column only contains the reference level.",
                 class = "offtake_bad_reference")
  }

  rows <- lapply(hunted_levs, function(g) {
    hv <- value[group == g]
    mh <- mean(hv, na.rm = TRUE)
    mr <- mean(ref_val, na.rm = TRUE)
    pct <- (mh - mr) / mr * 100

    # A one-sided Welch test is only possible with replicates in both groups.
    testable <- sum(!is.na(hv)) >= 2 && sum(!is.na(ref_val)) >= 2 &&
      stats::sd(hv, na.rm = TRUE) + stats::sd(ref_val, na.rm = TRUE) > 0
    if (testable) {
      alt <- if (lower_is_unsustainable) "less" else "greater"
      p <- tryCatch(
        stats::t.test(hv, ref_val, alternative = alt)$p.value,
        error = function(e) NA_real_
      )
    } else {
      p <- NA_real_
    }

    # Direction of concern.
    concerning <- if (lower_is_unsustainable) mh < mr else mh > mr
    if (is.na(p)) {
      # No statistical test possible; fall back to the raw direction and warn.
      sustainable <- !concerning
    } else {
      sustainable <- !(p < alpha && concerning)
    }
    tibble::tibble(
      hunted_site = g,
      reference_site = reference,
      hunted_value = mh,
      reference_value = mr,
      pct_change = pct,
      p_value = p,
      sustainable = sustainable
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "any_untestable") <- any(is.na(out$p_value))
  out
}

#' Population density comparison (PDC)
#'
#' Compares wildlife abundance/density between hunted and reference
#' (unhunted or lightly hunted) sites. Lower density at the hunted site(s) is
#' interpreted as unsustainable (Adounke et al. 2026; Weinbaum et al. 2013).
#'
#' Provide **one row per density estimate** (e.g. per line transect or
#' camera-trap station). With two or more replicates per site a one-sided Welch
#' *t*-test is used (hunted < reference); with a single value per site the
#' verdict falls back to the sign of the difference and a warning is issued.
#'
#' @param data A data frame in *long* format with one row per density estimate.
#'   It must contain at least: a numeric density column (`density`) and a
#'   site/treatment column (`group`) whose values include the reference level
#'   and one or more hunted levels. Extra columns are ignored.
#' @param density <[`data-masked`][rlang::args_data_masking]> Column of
#'   abundance or density estimates (e.g. individuals per km^2). Higher = more
#'   animals.
#' @param group <[`data-masked`][rlang::args_data_masking]> Column identifying
#'   the site or hunting treatment. Must contain the `reference` level and at
#'   least one hunted level.
#' @param reference Character scalar naming the level of `group` treated as the
#'   unhunted / lightly hunted reference (e.g. `"reference"`).
#' @param alpha Significance level for the one-sided test (default `0.05`).
#'
#' @return An [offtake][ot_pdc] tibble with **one row per hunted site** and the
#'   columns:
#' \describe{
#'   \item{hunted_site}{The hunted level of `group` being assessed.}
#'   \item{reference_site}{The reference level it is compared against.}
#'   \item{hunted_value}{Mean density at the hunted site.}
#'   \item{reference_value}{Mean density at the reference site.}
#'   \item{pct_change}{Percentage change of the hunted site relative to the
#'     reference, `(hunted - reference) / reference * 100`. Negative means the
#'     hunted site is depleted.}
#'   \item{p_value}{P-value of the one-sided Welch *t*-test that the hunted
#'     site has a *lower* mean than the reference. `NA` when a site has fewer
#'     than two replicates.}
#'   \item{sustainable}{Logical verdict. `FALSE` when the hunted site is
#'     significantly lower than the reference (`p_value < alpha`); otherwise
#'     `TRUE`.}
#' }
#'
#' @references
#' Adounke, G. R. M. et al. (2026) Systematic review of sustainability
#' assessment approaches for wildlife exploitation. *Biological Conservation*
#' 313, 111606. \doi{10.1016/j.biocon.2025.111606}
#'
#' Weinbaum, K. Z., Brashares, J. S., Golden, C. D. & Getz, W. M. (2013)
#' Searching for sustainability: are assessments of wildlife harvests behind
#' the times? *Ecology Letters* 16, 99-111. \doi{10.1111/ele.12008}
#'
#' @seealso [ot_hyco()], [ot_asc()]
#' @examples
#' # One row per transect. `site` = treatment, `dens` = animals per km^2.
#' d <- data.frame(
#'   site = rep(c("control", "hunted"), each = 4), # reference vs hunted
#'   dens = c(12, 14, 11, 13, 6, 7, 5, 8) # density per transect
#' )
#' ot_pdc(d, density = dens, group = site, reference = "control")
#' @export
ot_pdc <- function(data, density, group, reference, alpha = 0.05) {
  .check_data(data)
  value <- .pull_col(data, rlang::enquo(density), "density")
  grp <- .pull_col(data, rlang::enquo(group), "group")
  res <- .compare_groups(value, grp, reference, alpha,
                         lower_is_unsustainable = TRUE)
  if (isTRUE(attr(res, "any_untestable"))) {
    rlang::warn(
      "Some sites had < 2 replicates; verdict for those is based on the sign of the difference only, with no significance test.",
      class = "offtake_untestable"
    )
  }
  new_offtake(
    res,
    method = "PDC (population density comparison)",
    family = "index",
    reference = "Adounke et al. (2026); Weinbaum et al. (2013)"
  )
}

#' Hunting yield comparison (HYCo)
#'
#' Compares harvested biomass (or catch-per-unit-effort, CPUE) between more- and
#' less-hunted sites. Lower yields at the more-hunted site(s) are interpreted as
#' unsustainable (Adounke et al. 2026; Weinbaum et al. 2013). If an `effort`
#' column is supplied the metric compared is CPUE (`yield / effort`); otherwise
#' raw yield is compared.
#'
#' @param data A data frame in *long* format with one row per harvest record.
#'   It must contain at least a yield column (`yield`) and a site column
#'   (`group`); optionally a hunting-effort column (`effort`).
#' @param yield <[`data-masked`][rlang::args_data_masking]> Column of harvested
#'   biomass (e.g. kg) or number of animals taken.
#' @param group <[`data-masked`][rlang::args_data_masking]> Column identifying
#'   the site / hunting intensity.
#' @param reference Character scalar naming the level of `group` treated as the
#'   less-hunted reference (e.g. `"low"`).
#' @param effort <[`data-masked`][rlang::args_data_masking]> Optional column of
#'   hunting effort (e.g. hunter-days). When supplied, `yield / effort` (CPUE)
#'   is compared instead of raw yield.
#' @param alpha Significance level for the one-sided test (default `0.05`).
#'
#' @return An [offtake][ot_hyco] tibble with **one row per hunted site**. The
#'   columns are the same as for [ot_pdc()], except that `hunted_value` and
#'   `reference_value` hold the mean yield (or mean CPUE when `effort` is
#'   supplied) at each site rather than density; `sustainable` is `FALSE` when
#'   the more-hunted site yields significantly less.
#'
#' @inherit ot_pdc references
#' @seealso [ot_pdc()], [ot_asc()]
#' @examples
#' # One row per harvest record. `site` = intensity, `kg` = biomass taken,
#' # `days` = hunter-days of effort.
#' d <- data.frame(
#'   site = rep(c("low", "high"), each = 3),  # less- vs more-hunted
#'   kg = c(40, 45, 38, 20, 25, 18), # harvested biomass (kg)
#'   days = c(10, 11, 9, 10, 12, 9) # hunting effort (hunter-days)
#' )
#' # Compare CPUE (kg per hunter-day):
#' ot_hyco(d, yield = kg, group = site, reference = "low", effort = days)
#' @export
ot_hyco <- function(data, yield, group, reference, effort = NULL, alpha = 0.05) {
  .check_data(data)
  y <- .pull_col(data, rlang::enquo(yield), "yield")
  grp <- .pull_col(data, rlang::enquo(group), "group")
  eff <- .pull_col(data, rlang::enquo(effort), "effort", required = FALSE)
  value <- if (is.null(eff)) y else y / eff
  res <- .compare_groups(value, grp, reference, alpha,
                         lower_is_unsustainable = TRUE)
  if (isTRUE(attr(res, "any_untestable"))) {
    rlang::warn(
      "Some sites had < 2 replicates; verdict for those is based on the sign of the difference only, with no significance test.",
      class = "offtake_untestable"
    )
  }
  new_offtake(
    res,
    method = if (is.null(eff)) "HYCo (hunting yield comparison)"
             else "HYCo (hunting yield comparison, CPUE)",
    family = "index",
    reference = "Adounke et al. (2026); Weinbaum et al. (2013)"
  )
}

#' Age structure comparison (ASC)
#'
#' Compares the age/sex structure of a population between hunted and reference
#' (unhunted or lightly hunted) sites. Following Adounke et al. (2026), a
#' **lower proportion of the juvenile class at the hunted site is interpreted as
#' unsustainable** (a signal of recruitment failure). Optionally the full
#' class-by-site frequency distribution can be compared with a chi-squared test.
#'
#' @param data A data frame with **one row per sampled/harvested individual**
#'   (or one row per class if `count` is supplied). It must contain an age/sex
#'   class column (`stage`) and a site column (`group`). If your data are
#'   already tallied, add a `count` column and pass it to `count`.
#' @param stage <[`data-masked`][rlang::args_data_masking]> Column giving the
#'   age or sex class of each individual (e.g. `"juvenile"`/`"adult"`).
#' @param group <[`data-masked`][rlang::args_data_masking]> Column identifying
#'   the site / hunting treatment.
#' @param reference Character scalar naming the level of `group` treated as the
#'   unhunted / lightly hunted reference.
#' @param juvenile Character vector of the level(s) of `stage` that represent
#'   the juvenile (pre-reproductive) class. Required for
#'   `method = "proportion"`.
#' @param count <[`data-masked`][rlang::args_data_masking]> Optional column of
#'   counts, used when `data` is already aggregated to class \eqn{\times} site.
#' @param method `"proportion"` (default) compares the juvenile proportion
#'   between the two sites with a two-sample test of proportions; `"distribution"`
#'   compares the whole class \eqn{\times} site table with a chi-squared test.
#' @param alpha Significance level (default `0.05`).
#'
#' @return An [offtake][ot_asc] tibble with **one row per hunted site**. Columns
#'   depend on `method`.
#'
#'   For `method = "proportion"`:
#' \describe{
#'   \item{hunted_site, reference_site}{The two levels being compared.}
#'   \item{juv_prop_hunted}{Proportion of individuals in the juvenile class at
#'     the hunted site (0-1).}
#'   \item{juv_prop_reference}{Juvenile proportion at the reference site.}
#'   \item{pct_change}{Percentage change in juvenile proportion relative to the
#'     reference; negative means fewer juveniles at the hunted site.}
#'   \item{p_value}{One-sided [stats::prop.test()] p-value for a *lower*
#'     juvenile proportion at the hunted site.}
#'   \item{sustainable}{`FALSE` when the juvenile share is significantly lower
#'     at the hunted site.}
#' }
#'
#'   For `method = "distribution"`:
#' \describe{
#'   \item{hunted_site, reference_site}{The two levels being compared.}
#'   \item{statistic}{Pearson chi-squared statistic for the class
#'     \eqn{\times} site table.}
#'   \item{df}{Degrees of freedom of the test.}
#'   \item{p_value}{Chi-squared p-value.}
#'   \item{sustainable}{`FALSE` when the two age structures differ
#'     significantly (`p_value < alpha`); interpret alongside the direction of
#'     the shift (Weinbaum et al. 2013).}
#' }
#'
#' @inherit ot_pdc references
#' @seealso [ot_pdc()], [ot_hyco()]
#' @examples
#' # Individual-level data: one row per animal, `class` = age class of that
#' # animal, `site` = where it was sampled.
#' set.seed(1)
#' d <- data.frame(
#'   site = rep(c("control", "hunted"), c(60, 60)),
#'   class = c(sample(c("juv", "adult"), 60, TRUE, c(0.40, 0.60)),
#'             sample(c("juv", "adult"), 60, TRUE, c(0.15, 0.85)))
#' )
#' ot_asc(d, stage = class, group = site, reference = "control", juvenile = "juv")
#' @export
ot_asc <- function(data, stage, group, reference, juvenile = NULL,
                   count = NULL, method = c("proportion", "distribution"),
                   alpha = 0.05) {
  .check_data(data)
  method <- match.arg(method)
  st <- as.character(.pull_col(data, rlang::enquo(stage), "stage"))
  grp <- as.character(.pull_col(data, rlang::enquo(group), "group"))
  cnt <- .pull_col(data, rlang::enquo(count), "count", required = FALSE)
  if (is.null(cnt)) cnt <- rep(1, length(st))

  levs <- unique(grp)
  if (!reference %in% levs) {
    rlang::abort(
      sprintf("`reference` = \"%s\" is not a level of the grouping column.", reference),
      class = "offtake_bad_reference"
    )
  }
  hunted_levs <- setdiff(levs, reference)

  if (method == "proportion") {
    if (is.null(juvenile)) {
      rlang::abort("`juvenile` must be supplied when method = \"proportion\".",
                   class = "offtake_missing_arg")
    }
    is_juv <- st %in% juvenile
    tab <- function(g) {
      idx <- grp == g
      c(juv = sum(cnt[idx & is_juv]), tot = sum(cnt[idx]))
    }
    rt <- tab(reference)
    rows <- lapply(hunted_levs, function(g) {
      ht <- tab(g)
      # One-sided test: is the juvenile proportion LOWER at the hunted site?
      p <- tryCatch(
        stats::prop.test(c(ht["juv"], rt["juv"]), c(ht["tot"], rt["tot"]),
                         alternative = "less")$p.value,
        error = function(e) NA_real_
      )
      ph <- unname(ht["juv"] / ht["tot"]); pr <- unname(rt["juv"] / rt["tot"])
      concerning <- ph < pr
      sustainable <- if (is.na(p)) !concerning else !(p < alpha && concerning)
      tibble::tibble(
        hunted_site = g,
        reference_site = reference,
        juv_prop_hunted = ph,
        juv_prop_reference = pr,
        pct_change = (ph - pr) / pr * 100,
        p_value = p,
        sustainable = sustainable
      )
    })
    res <- do.call(rbind, rows)
    return(new_offtake(
      res,
      method = "ASC (age structure comparison, juvenile proportion)",
      family = "index",
      reference = "Adounke et al. (2026); Weinbaum et al. (2013)"
    ))
  }

  # method == "distribution": chi-squared of class x site (reference vs each hunted)
  rows <- lapply(hunted_levs, function(g) {
    idx <- grp %in% c(reference, g)
    m <- tapply(cnt[idx], list(st[idx], grp[idx]), sum)
    m[is.na(m)] <- 0
    ct <- suppressWarnings(stats::chisq.test(m))
    sustainable <- !(ct$p.value < alpha)  # significant difference => flagged
    tibble::tibble(
      hunted_site = g,
      reference_site = reference,
      statistic = unname(ct$statistic),
      df = unname(ct$parameter),
      p_value = ct$p.value,
      sustainable = sustainable
    )
  })
  res <- do.call(rbind, rows)
  new_offtake(
    res,
    method = "ASC (age structure comparison, distribution)",
    family = "index",
    reference = "Adounke et al. (2026); Weinbaum et al. (2013)",
    notes = "A significant class x site difference is flagged as unsustainable; interpret alongside the direction of the shift (Weinbaum et al. 2013)."
  )
}
