test_that("ot_pdc() flags lower density at hunted sites as unsustainable", {
  d <- data.frame(site = rep(c("control", "hunted"), each = 4),
                  dens = c(12, 14, 11, 13, 6, 7, 5, 8))
  res <- ot_pdc(d, density = dens, group = site, reference = "control")
  expect_s3_class(res, "offtake")
  expect_equal(res$hunted_site, "hunted")
  expect_lt(res$pct_change, 0)
  expect_lt(res$p_value, 0.05)
  expect_equal(res$outcome, "unsustainable")
  expect_false(res$sustainable)
})

test_that("ot_pdc() reports the log response ratio and its interval", {
  d <- data.frame(site = rep(c("control", "hunted"), each = 4),
                  dens = c(12, 14, 11, 13, 6, 7, 5, 8))
  res <- ot_pdc(d, density = dens, group = site, reference = "control")
  expect_equal(res$lnrr, log(6.5 / 12.5))
  # Delta-method standard error (Hedges et al. 1999).
  se <- sqrt(var(c(6, 7, 5, 8)) / (4 * 6.5^2) + var(c(12, 14, 11, 13)) / (4 * 12.5^2))
  expect_equal(res$lnrr_lo, log(6.5 / 12.5) - qnorm(0.975) * se)
  expect_equal(res$lnrr_hi, log(6.5 / 12.5) + qnorm(0.975) * se)
})

test_that("a precise absence of difference gives 'no substantial decline'", {
  d <- data.frame(site = rep(c("control", "hunted"), each = 6),
                  dens = c(20, 21, 19, 20, 21, 19, 20, 21, 19, 20, 21, 19))
  res <- ot_pdc(d, density = dens, group = site, reference = "control")
  expect_equal(res$outcome, "no substantial decline")
  expect_true(res$sustainable)
})

test_that("few or noisy replicates give 'inconclusive', not 'sustainable'", {
  # Two noisy replicates: not significant, but a large decline cannot be ruled out.
  d <- data.frame(site = c("r", "r", "h", "h"), v = c(10, 14, 11, 9))
  res <- ot_pdc(d, density = v, group = site, reference = "r")
  expect_equal(res$outcome, "inconclusive")
  expect_true(is.na(res$sustainable))
  # One value per site: no test at all.
  d1 <- data.frame(site = c("r", "h"), v = c(10, 4))
  expect_warning(res1 <- ot_pdc(d1, density = v, group = site, reference = "r"),
                 class = "offtake_untestable")
  expect_equal(res1$outcome, "inconclusive")
})

test_that("p-values are adjusted across several hunted sites", {
  set.seed(3)
  d <- data.frame(
    site = rep(c("ref", "h1", "h2", "h3"), each = 5),
    dens = c(rnorm(5, 20, 2), rnorm(5, 17, 2), rnorm(5, 18, 2), rnorm(5, 15, 2))
  )
  res <- ot_pdc(d, density = dens, group = site, reference = "ref")
  expect_equal(res$p_adj, p.adjust(res$p_value, "holm"))
  none <- ot_pdc(d, density = dens, group = site, reference = "ref", p_adjust = "none")
  expect_equal(none$p_adj, none$p_value)
})

test_that("ot_pdc() analyses grouped data group by group", {
  skip_if_not_installed("dplyr")
  d <- rbind(
    data.frame(sp = "A", site = rep(c("ref", "hunt"), each = 4),
               dens = c(12, 14, 11, 13, 6, 7, 5, 8)),
    data.frame(sp = "B", site = rep(c("ref", "hunt"), each = 4),
               dens = c(20, 21, 19, 20, 20, 21, 19, 20))
  )
  res <- ot_pdc(dplyr::group_by(d, sp), density = dens, group = site, reference = "ref")
  expect_equal(res$sp, c("A", "B"))
  expect_equal(res$outcome[1], "unsustainable")
})

test_that("ot_hyco() compares CPUE when effort is supplied", {
  d <- data.frame(site = rep(c("low", "high"), each = 3),
                  kg = c(40, 45, 38, 20, 25, 18),
                  days = c(10, 11, 9, 10, 12, 9))
  res_raw <- ot_hyco(d, yield = kg, group = site, reference = "low")
  res_cpue <- ot_hyco(d, yield = kg, group = site, reference = "low", effort = days)
  expect_false(res_raw$sustainable)
  expect_false(res_cpue$sustainable)
  expect_equal(res_cpue$hunted_value, mean(c(20 / 10, 25 / 12, 18 / 9)))
})

test_that("indices check their inputs", {
  d <- data.frame(site = rep(c("a", "b"), each = 3), v = 1:6)
  expect_error(ot_pdc(d, density = v, group = site, reference = "zzz"),
               class = "offtake_bad_reference")
  expect_error(ot_pdc(d, density = -v, group = site, reference = "a"),
               class = "offtake_bad_value")
  expect_error(ot_pdc(d, density = v, group = site, reference = "a", decline = 2),
               class = "offtake_bad_value")
  expect_error(ot_hyco(d, yield = v, group = site, reference = "a", effort = v * 0),
               class = "offtake_bad_value")
})

test_that("ot_asc() flags a lower juvenile proportion at the hunted site", {
  res <- ot_asc(bushmeat_ages, stage = age_class, group = site_type,
                reference = "reference", juvenile = "juvenile")
  expect_lt(res$juv_prop_hunted, res$juv_prop_reference)
  expect_equal(res$outcome, "unsustainable")
  # Log ratio of the two proportions.
  expect_equal(res$lnrr, log(res$juv_prop_hunted / res$juv_prop_reference))
})

test_that("ot_asc() distribution method never claims sustainability", {
  set.seed(2)
  d <- data.frame(
    site = rep(c("control", "hunted"), c(80, 80)),
    class = c(sample(c("juv", "sub", "adult"), 80, TRUE),
              sample(c("juv", "sub", "adult"), 80, TRUE))
  )
  res <- ot_asc(d, stage = class, group = site, reference = "control",
                method = "distribution")
  expect_true(all(c("statistic", "df", "p_value", "p_adj") %in% names(res)))
  expect_true(res$outcome %in% c("unsustainable", "inconclusive"))
  expect_false(isTRUE(res$sustainable))
})

test_that("ot_asc() accepts tallied counts", {
  tallied <- data.frame(site = c("ref", "ref", "hunt", "hunt"),
                        class = c("juv", "adult", "juv", "adult"),
                        n = c(45, 75, 19, 101))
  res <- ot_asc(tallied, stage = class, group = site, reference = "ref",
                juvenile = "juv", count = n)
  expect_equal(res$juv_prop_hunted, 19 / 120)
  expect_equal(res$n_reference, 120)
})
