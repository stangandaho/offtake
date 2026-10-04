# Shared comparison engine for the index-based approaches
#
# All three indices catalogued by Adounke et al. (2026) that are implemented
# here (ASC, HYCo, PDC) share the same logic: a metric measured at one or more
# *hunted* sites is compared with the same metric at a *reference*
# (unhunted / lightly hunted) site, and a lower value at the hunted site is
# interpreted as unsustainable (Weinbaum et al. 2013, Table 1).

# Three-level outcome shared by the indices.
# - "unsustainable": significantly lower at the hunted site (adjusted p < alpha);
# - "no substantial decline": not significant AND the confidence interval of the
#   log response ratio excludes a decline of `decline` or more;
# - "inconclusive": everything else (no test possible, or an interval too wide
#   to rule out a substantial decline).
.index_outcome <- function(p_adj, lnrr_lo, concerning, alpha, decline) {
  unsus <- !is.na(p_adj) & p_adj < alpha & concerning
  no_dec <- !unsus & !is.na(lnrr_lo) & lnrr_lo >= log(1 - decline)
  outcome <- rep("inconclusive", length(p_adj))
  outcome[unsus] <- "unsustainable"
  outcome[no_dec] <- "no substantial decline"
  sustainable <- rep(NA, length(p_adj))
  sustainable[unsus] <- FALSE
  sustainable[no_dec] <- TRUE
  list(outcome = outcome, sustainable = sustainable)
}

.check_index_args <- function(alpha, decline, p_adjust) {
  .check_range(alpha, "alpha", 0, 1, lower_open = TRUE)
  .check_range(decline, "decline", 0, 1, lower_open = TRUE)
  if (!p_adjust %in% stats::p.adjust.methods) {
    rlang::abort(sprintf("`p_adjust` must be one of: %s.",
                         paste(stats::p.adjust.methods, collapse = ", ")),
                 class = "offtake_bad_value")
  }
}

.check_reference <- function(group, reference) {
  levs <- unique(group)
  if (!reference %in% levs) {
    rlang::abort(
      sprintf("`reference` = \"%s\" is not a level of the grouping column (levels: %s).",
              reference, paste(levs, collapse = ", ")),
      class = "offtake_bad_reference"
    )
  }
  hunted <- setdiff(levs, reference)
  if (length(hunted) == 0L) {
    rlang::abort("No hunted group found: the grouping column only contains the reference level.",
                 class = "offtake_bad_reference")
  }
  hunted
}

# Compare a continuous metric between each hunted group and the reference.
.compare_groups <- function(value, group, reference, alpha, p_adjust, decline) {
  group <- as.character(group)
  hunted_levs <- .check_reference(group, reference)
  rv <- value[group == reference & !is.na(value)]
  z <- stats::qnorm(1 - alpha / 2)

  rows <- lapply(hunted_levs, function(g) {
    hv <- value[group == g & !is.na(value)]
    nh <- length(hv); nr <- length(rv)
    mh <- mean(hv); mr <- mean(rv)
    sh <- if (nh >= 2) stats::sd(hv) else NA_real_
    sr <- if (nr >= 2) stats::sd(rv) else NA_real_

    # A one-sided Welch test is only possible with replicates in both groups.
    testable <- nh >= 2 && nr >= 2 && (sh + sr) > 0
    p <- if (testable) {
      tryCatch(stats::t.test(hv, rv, alternative = "less")$p.value,
               error = function(e) NA_real_)
    } else NA_real_

    # Log response ratio and its delta-method standard error (Hedges et al. 1999).
    lnrr <- if (mh > 0 && mr > 0) log(mh / mr) else NA_real_
    se <- if (testable && !is.na(lnrr)) sqrt(sh^2 / (nh * mh^2) + sr^2 / (nr * mr^2)) else NA_real_

    tibble::tibble(
      hunted_site = g,
      reference_site = reference,
      n_hunted = nh,
      n_reference = nr,
      hunted_value = mh,
      reference_value = mr,
      pct_change = (mh - mr) / mr * 100,
      lnrr = lnrr,
      lnrr_lo = lnrr - z * se,
      lnrr_hi = lnrr + z * se,
      p_value = p
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  verdict <- .index_outcome(out$p_adj, out$lnrr_lo, out$hunted_value < out$reference_value,
                            alpha, decline)
  out$outcome <- verdict$outcome
  out$sustainable <- verdict$sustainable
  out
}

.warn_untestable <- function(res) {
  if (any(is.na(res$p_value))) {
    rlang::warn(
      "Some sites had fewer than 2 replicates (or no variation); no test was possible and their outcome is \"inconclusive\".",
      class = "offtake_untestable"
    )
  }
}

#' Population density comparison (PDC)
#'
#' Compares wildlife abundance/density between hunted and reference
#' (unhunted or lightly hunted) sites. Lower density at the hunted site(s) is
#' interpreted as unsustainable (Adounke et al. 2026; Weinbaum et al. 2013).
#'
#' Provide one row per density estimate (e.g. per line transect or
#' camera-trap station). For each hunted site the function reports the size of
#' the difference (percent change and log response ratio with its confidence
#' interval) and a one-sided Welch *t*-test (hunted < reference). When several
#' hunted sites are compared with the same reference, p-values are adjusted for
#' multiple comparisons (`p_adjust`).
#'
#' @section Outcome:
#' The verdict has three levels, so that "no difference found" is not mistaken
#' for evidence of sustainability:
#' * `"unsustainable"`: the hunted site is significantly lower (adjusted
#'   p-value below `alpha`); `sustainable = FALSE`.
#' * `"no substantial decline"`: not significant, and the confidence interval
#'   of the response ratio rules out a decline of `decline` (default 20%) or
#'   more; `sustainable = TRUE`.
#' * `"inconclusive"`: no test was possible (fewer than two replicates per
#'   site) or the interval is too wide to rule out a substantial decline;
#'   `sustainable = NA`.
#'
#' @section Pseudoreplication:
#' Replicates taken inside one hunted site and one reference site (transects,
#' camera stations) are not independent replicates of hunting. The test then
#' shows that the two *sites* differ, not that *hunting* causes the
#' difference (Hurlbert 1984). Prefer several hunted and several reference
#' sites, keep the sites otherwise comparable, and read the result as a
#' warning sign rather than a proof.
#'
#' Grouped data frames (from `dplyr::group_by()`) are analysed group by group,
#' for example one comparison per species, and the group columns are kept in
#' the output.
#'
#' @param data A data frame in *long* format with one row per density estimate.
#'   It must contain at least: a numeric density column (`density`) and a
#'   site/treatment column (`group`) whose values include the reference level
#'   and one or more hunted levels. Extra columns are ignored. May be grouped
#'   with `dplyr::group_by()`.
#' @param density <[`data-masked`][rlang::args_data_masking]> Column of
#'   abundance or density estimates (e.g. individuals per km^2). Higher = more
#'   animals.
#' @param group <[`data-masked`][rlang::args_data_masking]> Column identifying
#'   the site or hunting treatment. Must contain the `reference` level and at
#'   least one hunted level.
#' @param reference Character scalar naming the level of `group` treated as the
#'   unhunted / lightly hunted reference (e.g. `"reference"`).
#' @param alpha Significance level of the one-sided test; the confidence
#'   intervals have level `1 - alpha` (default `0.05`).
#' @param p_adjust Method used to adjust p-values when several hunted sites are
#'   compared with the same reference, passed to [stats::p.adjust()] (default
#'   `"holm"`; `"none"` for no adjustment).
#' @param decline Smallest relative decline considered substantial, between 0
#'   and 1 (default `0.2`, i.e. 20%). Used to tell "no substantial decline"
#'   from "inconclusive".
#'
#' @return An [offtake][ot_pdc] tibble with one row per hunted site (and per
#'   group, if `data` is grouped) and the columns:
#' \describe{
#'   \item{hunted_site}{The hunted level of `group` being assessed.}
#'   \item{reference_site}{The reference level it is compared against.}
#'   \item{n_hunted, n_reference}{Number of non-missing values at each site.}
#'   \item{hunted_value}{Mean density at the hunted site.}
#'   \item{reference_value}{Mean density at the reference site.}
#'   \item{pct_change}{Percentage change of the hunted site relative to the
#'     reference, `(hunted - reference) / reference * 100`. Negative means the
#'     hunted site is depleted.}
#'   \item{lnrr}{Log response ratio, `log(hunted / reference)`. A value of
#'     `-0.69` means the hunted site has half the reference value.}
#'   \item{lnrr_lo, lnrr_hi}{Confidence interval of `lnrr` (level
#'     `1 - alpha`), from the delta-method variance of Hedges et al. (1999).
#'     `NA` when a site has fewer than two replicates.}
#'   \item{p_value}{P-value of the one-sided Welch *t*-test that the hunted
#'     site has a *lower* mean than the reference. `NA` when no test is
#'     possible.}
#'   \item{p_adj}{`p_value` adjusted across the hunted sites with `p_adjust`.}
#'   \item{outcome}{`"unsustainable"`, `"no substantial decline"` or
#'     `"inconclusive"` (see the Outcome section).}
#'   \item{sustainable}{`FALSE`, `TRUE` or `NA`, matching `outcome`.}
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
#' Hedges, L. V., Gurevitch, J. & Curtis, P. S. (1999) The meta-analysis of
#' response ratios in experimental ecology. *Ecology* 80, 1150-1156.
#'
#' Hurlbert, S. H. (1984) Pseudoreplication and the design of ecological field
#' experiments. *Ecological Monographs* 54, 187-211.
#'
#' @seealso [ot_hyco()], [ot_asc()]
#' @examples
#' # One row per transect. `site` = treatment, `dens` = animals per km^2.
#' d <- data.frame(
#'   site = rep(c("control", "hunted"), each = 4), # reference vs hunted
#'   dens = c(12, 14, 11, 13, 6, 7, 5, 8) # density per transect
#' )
#' ot_pdc(d, density = dens, group = site, reference = "control")
#'
#' # The bundled example data
#' ot_pdc(bushmeat_sites, density = density, group = site_type,
#'        reference = "reference")
#' @export
ot_pdc <- function(data, density, group, reference, alpha = 0.05,
                   p_adjust = "holm", decline = 0.2) {
  .check_data(data)
  .check_index_args(alpha, decline, p_adjust)
  groups <- .group_info(data)
  data <- .ungroup(data)
  value <- .pull_col(data, rlang::enquo(density), "density")
  .check_range(value, "density", 0)
  grp <- as.character(.pull_col(data, rlang::enquo(group), "group"))

  run <- function(rows) .compare_groups(value[rows], grp[rows], reference,
                                        alpha, p_adjust, decline)
  res <- if (is.null(groups)) run(seq_len(nrow(data))) else .by_group(groups, run)
  .warn_untestable(res)
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
#' column is supplied the metric compared is CPUE (`yield / effort`, one value
#' per record); otherwise raw yield is compared. The statistics, the outcome
#' and the caveats are those of [ot_pdc()].
#'
#' @inheritParams ot_pdc
#' @inheritSection ot_pdc Outcome
#' @inheritSection ot_pdc Pseudoreplication
#' @param data A data frame in *long* format with one row per harvest record.
#'   It must contain at least a yield column (`yield`) and a site column
#'   (`group`); optionally a hunting-effort column (`effort`). May be grouped
#'   with `dplyr::group_by()`.
#' @param yield <[`data-masked`][rlang::args_data_masking]> Column of harvested
#'   biomass (e.g. kg) or number of animals taken.
#' @param group <[`data-masked`][rlang::args_data_masking]> Column identifying
#'   the site / hunting intensity.
#' @param reference Character scalar naming the level of `group` treated as the
#'   less-hunted reference (e.g. `"low"`).
#' @param effort <[`data-masked`][rlang::args_data_masking]> Optional column of
#'   hunting effort (e.g. hunter-days), strictly positive. When supplied,
#'   `yield / effort` (CPUE) is compared instead of raw yield.
#'
#' @return An [offtake][ot_hyco] tibble with one row per hunted site. The
#'   columns are the same as for [ot_pdc()], except that `hunted_value` and
#'   `reference_value` hold the mean yield (or mean CPUE when `effort` is
#'   supplied) at each site rather than density.
#'
#' @inherit ot_pdc references
#' @seealso [ot_pdc()], [ot_asc()]
#' @examples
#' # One row per harvest record. `site` = intensity, `kg` = biomass taken,
#' # `days` = hunter-days of effort.
#' d <- data.frame(
#'   site = rep(c("low", "high"), each = 3), # less- vs more-hunted
#'   kg = c(40, 45, 38, 20, 25, 18), # harvested biomass (kg)
#'   days = c(10, 11, 9, 10, 12, 9) # hunting effort (hunter-days)
#' )
#' # Compare CPUE (kg per hunter-day):
#' ot_hyco(d, yield = kg, group = site, reference = "low", effort = days)
#' @export
ot_hyco <- function(data, yield, group, reference, effort = NULL, alpha = 0.05,
                    p_adjust = "holm", decline = 0.2) {
  .check_data(data)
  .check_index_args(alpha, decline, p_adjust)
  groups <- .group_info(data)
  data <- .ungroup(data)
  y <- .pull_col(data, rlang::enquo(yield), "yield")
  .check_range(y, "yield", 0)
  grp <- as.character(.pull_col(data, rlang::enquo(group), "group"))
  eff <- .pull_col(data, rlang::enquo(effort), "effort", required = FALSE)
  if (!is.null(eff)) .check_range(eff, "effort", 0, lower_open = TRUE)
  value <- if (is.null(eff)) y else y / eff

  run <- function(rows) .compare_groups(value[rows], grp[rows], reference,
                                        alpha, p_adjust, decline)
  res <- if (is.null(groups)) run(seq_len(nrow(data))) else .by_group(groups, run)
  .warn_untestable(res)
  new_offtake(
    res,
    method = if (is.null(eff)) "HYCo (hunting yield comparison)"
             else "HYCo (hunting yield comparison, CPUE)",
    family = "index",
    reference = "Adounke et al. (2026); Weinbaum et al. (2013)"
  )
}

# Juvenile share at each hunted site against the reference.
.asc_proportion <- function(st, grp, cnt, reference, juvenile, alpha, p_adjust, decline) {
  hunted_levs <- .check_reference(grp, reference)
  is_juv <- st %in% juvenile
  tab <- function(g) {
    idx <- grp == g
    c(juv = sum(cnt[idx & is_juv]), tot = sum(cnt[idx]))
  }
  rt <- tab(reference)
  z <- stats::qnorm(1 - alpha / 2)
  rows <- lapply(hunted_levs, function(g) {
    ht <- tab(g)
    # One-sided test: is the juvenile proportion LOWER at the hunted site?
    p <- tryCatch(
      stats::prop.test(c(ht[["juv"]], rt[["juv"]]), c(ht[["tot"]], rt[["tot"]]),
                       alternative = "less")$p.value,
      error = function(e) NA_real_
    )
    ph <- ht[["juv"]] / ht[["tot"]]
    pr <- rt[["juv"]] / rt[["tot"]]
    # Log ratio of the two proportions, with a 0.5 correction when a juvenile
    # count is zero.
    a <- ht[["juv"]]; n1 <- ht[["tot"]]; cc <- rt[["juv"]]; n2 <- rt[["tot"]]
    if (n1 > 0 && n2 > 0 && (a == 0 || cc == 0)) {
      a <- a + 0.5; cc <- cc + 0.5; n1 <- n1 + 1; n2 <- n2 + 1
    }
    ok <- n1 > 0 && n2 > 0
    lnrr <- if (ok) log((a / n1) / (cc / n2)) else NA_real_
    se <- if (ok) sqrt(1 / a - 1 / n1 + 1 / cc - 1 / n2) else NA_real_
    tibble::tibble(
      hunted_site = g,
      reference_site = reference,
      n_hunted = ht[["tot"]],
      n_reference = rt[["tot"]],
      juv_prop_hunted = ph,
      juv_prop_reference = pr,
      pct_change = (ph - pr) / pr * 100,
      lnrr = lnrr,
      lnrr_lo = lnrr - z * se,
      lnrr_hi = lnrr + z * se,
      p_value = p
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  verdict <- .index_outcome(out$p_adj, out$lnrr_lo,
                            out$juv_prop_hunted < out$juv_prop_reference, alpha, decline)
  out$outcome <- verdict$outcome
  out$sustainable <- verdict$sustainable
  out
}

# Whole class x site distribution, chi-squared test per hunted site.
.asc_distribution <- function(st, grp, cnt, reference, alpha, p_adjust) {
  hunted_levs <- .check_reference(grp, reference)
  rows <- lapply(hunted_levs, function(g) {
    idx <- grp %in% c(reference, g)
    m <- tapply(cnt[idx], list(st[idx], grp[idx]), sum)
    m[is.na(m)] <- 0
    ct <- suppressWarnings(stats::chisq.test(m))
    tibble::tibble(
      hunted_site = g,
      reference_site = reference,
      statistic = unname(ct$statistic),
      df = unname(ct$parameter),
      p_value = ct$p.value
    )
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  differs <- !is.na(out$p_adj) & out$p_adj < alpha
  out$outcome <- ifelse(differs, "unsustainable", "inconclusive")
  out$sustainable <- ifelse(differs, FALSE, NA)
  out
}

#' Age structure comparison (ASC)
#'
#' Compares the age/sex structure of a population between hunted and reference
#' (unhunted or lightly hunted) sites. Following Adounke et al. (2026), a
#' **lower proportion of the juvenile class at the hunted site is interpreted as
#' unsustainable** (a signal of recruitment failure). Optionally the full
#' class-by-site frequency distribution can be compared with a chi-squared test.
#'
#' With `method = "proportion"` the size of the difference is reported as a log
#' ratio of the two juvenile proportions with its confidence interval, the test
#' is a one-sided two-sample test of proportions, and the three-level outcome
#' is that of [ot_pdc()]. With `method = "distribution"` a chi-squared test
#' cannot tell the direction of the shift, nor show that two structures are
#' similar: a significant difference gives `"unsustainable"` (as in Weinbaum
#' et al. 2013) and anything else `"inconclusive"`.
#'
#' @inheritParams ot_pdc
#' @inheritSection ot_pdc Pseudoreplication
#' @param data A data frame with one row per sampled/harvested individual
#'   (or one row per class if `count` is supplied). It must contain an age/sex
#'   class column (`stage`) and a site column (`group`). If your data are
#'   already tallied, add a `count` column and pass it to `count`. May be
#'   grouped with `dplyr::group_by()`.
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
#'   between the two sites; `"distribution"` compares the whole class
#'   \eqn{\times} site table with a chi-squared test.
#'
#' @return An [offtake][ot_asc] tibble with one row per hunted site (and per
#'   group, if `data` is grouped). Columns depend on `method`.
#'
#'   For `method = "proportion"`:
#' \describe{
#'   \item{hunted_site, reference_site}{The two levels being compared.}
#'   \item{n_hunted, n_reference}{Number of individuals at each site.}
#'   \item{juv_prop_hunted}{Proportion of individuals in the juvenile class at
#'     the hunted site (0-1).}
#'   \item{juv_prop_reference}{Juvenile proportion at the reference site.}
#'   \item{pct_change}{Percentage change in juvenile proportion relative to the
#'     reference; negative means fewer juveniles at the hunted site.}
#'   \item{lnrr, lnrr_lo, lnrr_hi}{Log ratio of the two juvenile proportions and
#'     its confidence interval (level `1 - alpha`).}
#'   \item{p_value}{One-sided [stats::prop.test()] p-value for a *lower*
#'     juvenile proportion at the hunted site.}
#'   \item{p_adj}{`p_value` adjusted across the hunted sites with `p_adjust`.}
#'   \item{outcome, sustainable}{Three-level verdict, as in [ot_pdc()].}
#' }
#'
#'   For `method = "distribution"`:
#' \describe{
#'   \item{hunted_site, reference_site}{The two levels being compared.}
#'   \item{statistic}{Pearson chi-squared statistic for the class
#'     \eqn{\times} site table.}
#'   \item{df}{Degrees of freedom of the test.}
#'   \item{p_value, p_adj}{Chi-squared p-value, raw and adjusted.}
#'   \item{outcome}{`"unsustainable"` when the two structures differ
#'     significantly, otherwise `"inconclusive"`.}
#'   \item{sustainable}{`FALSE` or `NA`, matching `outcome`.}
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
                   alpha = 0.05, p_adjust = "holm", decline = 0.2) {
  .check_data(data)
  method <- match.arg(method)
  .check_index_args(alpha, decline, p_adjust)
  groups <- .group_info(data)
  data <- .ungroup(data)
  st <- as.character(.pull_col(data, rlang::enquo(stage), "stage"))
  grp <- as.character(.pull_col(data, rlang::enquo(group), "group"))
  cnt <- .pull_col(data, rlang::enquo(count), "count", required = FALSE)
  if (is.null(cnt)) cnt <- rep(1, length(st))
  .check_range(cnt, "count", 0)
  if (method == "proportion" && is.null(juvenile)) {
    rlang::abort("`juvenile` must be supplied when method = \"proportion\".",
                 class = "offtake_missing_arg")
  }

  run <- function(rows) {
    if (method == "proportion") {
      .asc_proportion(st[rows], grp[rows], cnt[rows], reference, juvenile,
                      alpha, p_adjust, decline)
    } else {
      .asc_distribution(st[rows], grp[rows], cnt[rows], reference, alpha, p_adjust)
    }
  }
  res <- if (is.null(groups)) run(seq_len(nrow(data))) else .by_group(groups, run)

  if (method == "proportion") {
    return(new_offtake(
      res,
      method = "ASC (age structure comparison, juvenile proportion)",
      family = "index",
      reference = "Adounke et al. (2026); Weinbaum et al. (2013)"
    ))
  }
  new_offtake(
    res,
    method = "ASC (age structure comparison, distribution)",
    family = "index",
    reference = "Adounke et al. (2026); Weinbaum et al. (2013)",
    notes = "A chi-squared test cannot give the direction of the shift nor show similarity: a significant difference is flagged, anything else is inconclusive."
  )
}
