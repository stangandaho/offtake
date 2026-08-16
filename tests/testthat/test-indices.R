test_that("ot_pdc() flags lower density at hunted sites as unsustainable", {
  d <- data.frame(site = rep(c("control", "hunted"), each = 4),
                  dens = c(12, 14, 11, 13, 6, 7, 5, 8))
  res <- ot_pdc(d, density = dens, group = site, reference = "control")
  expect_s3_class(res, "offtake")
  expect_equal(res$hunted_site, "hunted")
  expect_lt(res$pct_change, 0)          # hunted lower
  expect_lt(res$p_value, 0.05)
  expect_false(res$sustainable)
})

test_that("ot_pdc() does not flag when hunted site is not lower", {
  d <- data.frame(site = rep(c("control", "hunted"), each = 4),
                  dens = c(12, 14, 11, 13, 12, 14, 11, 13))
  res <- ot_pdc(d, density = dens, group = site, reference = "control")
  expect_true(res$sustainable)
})

test_that("ot_hyco() compares CPUE when effort is supplied", {
  d <- data.frame(site = rep(c("low", "high"), each = 3),
                  kg = c(40, 45, 38, 20, 25, 18),
                  days = c(10, 11, 9, 10, 12, 9))
  res_raw <- ot_hyco(d, yield = kg, group = site, reference = "low")
  res_cpue <- ot_hyco(d, yield = kg, group = site, reference = "low", effort = days)
  expect_false(res_raw$sustainable)
  expect_false(res_cpue$sustainable)
  # CPUE value differs from raw-yield value.
  expect_false(isTRUE(all.equal(res_raw$hunted_value, res_cpue$hunted_value)))
})

test_that("ot_hyco() warns and still returns when a site lacks replicates", {
  d <- data.frame(site = c("low", "high"), kg = c(40, 20))
  expect_warning(
    res <- ot_hyco(d, yield = kg, group = site, reference = "low"),
    class = "offtake_untestable"
  )
  expect_true(is.na(res$p_value))
  expect_false(res$sustainable)         # falls back to sign of difference
})

test_that("ot_asc() flags a lower juvenile proportion at the hunted site", {
  set.seed(1)
  d <- data.frame(
    site = rep(c("control", "hunted"), c(80, 80)),
    class = c(sample(c("juv", "adult"), 80, TRUE, c(0.45, 0.55)),
              sample(c("juv", "adult"), 80, TRUE, c(0.12, 0.88)))
  )
  res <- ot_asc(d, stage = class, group = site, reference = "control",
                juvenile = "juv")
  expect_lt(res$juv_prop_hunted, res$juv_prop_reference)
  expect_false(res$sustainable)
})

test_that("ot_asc() distribution method returns a chi-squared test", {
  set.seed(2)
  d <- data.frame(
    site = rep(c("control", "hunted"), c(80, 80)),
    class = c(sample(c("juv", "sub", "adult"), 80, TRUE),
              sample(c("juv", "sub", "adult"), 80, TRUE, c(0.1, 0.2, 0.7)))
  )
  res <- ot_asc(d, stage = class, group = site, reference = "control",
                method = "distribution")
  expect_true(all(c("statistic", "df", "p_value") %in% names(res)))
})

test_that("bad reference level is an error", {
  d <- data.frame(site = rep(c("a", "b"), each = 3), v = 1:6)
  expect_error(ot_pdc(d, density = v, group = site, reference = "zzz"),
               class = "offtake_bad_reference")
})
