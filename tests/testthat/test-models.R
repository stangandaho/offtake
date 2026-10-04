test_that("ot_lambda_max solves Cole's equation and is monotonic in fecundity", {
  # Higher fecundity -> higher lambda_max.
  expect_gt(ot_lambda_max(0.8, 1, 8), ot_lambda_max(0.4, 1, 8))
  # A growing population has lambda > 1.
  expect_gt(ot_lambda_max(0.6, 1, 8), 1)
  # Cole's equation should be (near) satisfied at the returned root.
  lam <- ot_lambda_max(0.6, 1, 8); r <- log(lam)
  f <- exp(-r) + 0.6 * exp(-r * 1) - 0.6 * exp(-r * (8 + 1)) - 1
  expect_lt(abs(f), 1e-6)
})

test_that("ot_lambda_max matches the analytic limit of Cole's equation", {
  # With a = 1 and a very late last reproduction, Cole's equation reduces to
  # 1 = e^-r (1 + b), so lambda_max = 1 + b.
  expect_equal(ot_lambda_max(0.6, 1, 500), 1.6, tolerance = 1e-6)
})

test_that("ot_pro() matches the Robinson & Redford formula and flags overharvest", {
  d <- data.frame(K = 100, lam = 1.5, F = 0.4, H = 50)
  res <- ot_pro(d, k = K, harvest = H, lambda = lam, f = F)
  # P = 0.6 * 100 * (1.5 - 1) * 0.4 = 12
  expect_equal(res$limit, 12)
  expect_equal(res$observed, 50)
  expect_equal(res$ratio, 50 / 12)
  expect_false(res$sustainable)
  res2 <- ot_pro(data.frame(K = 100, lam = 1.5, F = 0.4, H = 10),
                 k = K, harvest = H, lambda = lam, f = F)
  expect_true(res2$sustainable)
})

test_that("ot_pro() derives F from longevity", {
  d <- data.frame(K = 100, lam = 1.5, H = 1, life = c(4, 8, 15))
  res <- ot_pro(d, k = K, harvest = H, lambda = lam, longevity = life)
  expect_equal(res$f, c(0.6, 0.4, 0.2))
})

test_that("ot_pro() checks its inputs", {
  d <- data.frame(K = 100, lam = 1.5, F = 0.4, H = 10)
  expect_error(ot_pro(d, k = -K, harvest = H, lambda = lam, f = F), class = "offtake_bad_value")
  expect_error(ot_pro(d, k = K, harvest = H, lambda = lam, f = 2), class = "offtake_bad_value")
  expect_warning(ot_pro(d, k = K, harvest = H, lambda = 0.9, f = F), class = "offtake_no_surplus")
  expect_error(ot_pro(d, k = K, harvest = H, lambda = lam, f = F, uncertainty = list(zz = 0.1)),
               class = "offtake_bad_value")
})

test_that("uncertainty propagation gives a probability and an interval", {
  d <- data.frame(K = 100, lam = 1.5, F = 0.4, H = 12)
  res <- ot_pro(d, k = K, harvest = H, lambda = lam, f = F,
                uncertainty = list(k = 0.3, harvest = 0.3), n_sim = 4000, seed = 1)
  expect_true(all(c("limit_lo", "limit_hi", "p_unsustainable") %in% names(res)))
  expect_lt(res$limit_lo, res$limit)
  expect_gt(res$limit_hi, res$limit)
  # Observed = limit and symmetric uncertainty: about half the draws exceed.
  expect_gt(res$p_unsustainable, 0.35)
  expect_lt(res$p_unsustainable, 0.65)
  # Reproducible with a seed.
  res2 <- ot_pro(d, k = K, harvest = H, lambda = lam, f = F,
                 uncertainty = list(k = 0.3, harvest = 0.3), n_sim = 4000, seed = 1)
  expect_equal(res$p_unsustainable, res2$p_unsustainable)
  # Coefficients of variation can come from a column.
  d$k_cv <- 0.2
  res3 <- ot_pro(d, k = K, harvest = H, lambda = lam, f = F,
                 uncertainty = list(k = k_cv), n_sim = 500, seed = 1)
  expect_true(is.numeric(res3$p_unsustainable))
})

test_that("model functions keep group columns of grouped data", {
  skip_if_not_installed("dplyr")
  res <- ot_pro(dplyr::group_by(duiker_demography, species),
                k = density_k, harvest = annual_take, b = b, a = a, w = w,
                longevity = lifespan)
  expect_equal(res$species, duiker_demography$species)
})

test_that("ot_pbr() matches Wade (1998) and derives Nmin from n/cv", {
  res <- ot_pbr(data.frame(N = 1000, rmax = 0.04, take = 5), rmax = rmax,
                removal = take, n = N, cv = 0, fr = 0.5)
  # With cv = 0, Nmin = N; PBR = 1000 * 0.5 * 0.04 * 0.5 = 10
  expect_equal(res$limit, 10)
  expect_true(res$sustainable)
  # Wade's 20th percentile: Nmin = N / exp(0.842 * sqrt(log(1 + CV^2))).
  res_cv <- ot_pbr(data.frame(N = 1000, rmax = 0.04, take = 5), rmax = rmax,
                   removal = take, n = N, cv = 0.3, fr = 0.5)
  expect_equal(res_cv$nmin, 1000 / exp(0.8416212 * sqrt(log(1 + 0.3^2))), tolerance = 1e-6)
})

test_that("ot_pbr() equals ot_msy() when N = K, CV = 0 and FR = 0.5", {
  d <- data.frame(K = 40, r = 0.3, take = 2)
  pbr <- ot_pbr(d, rmax = r, removal = take, n = K, fr = 0.5)
  msy <- ot_msy(d, r = r, k = K, harvest = take)
  expect_equal(pbr$limit, msy$limit)
})

test_that("ot_msy() equals rK/4 and reports surplus at n", {
  res <- ot_msy(data.frame(r = 0.4, K = 1000, H = 80, N = 500),
                r = r, k = K, harvest = H, n = N)
  expect_equal(res$limit, 100)
  expect_equal(res$surplus_at_n, 100)
  expect_true(res$sustainable)
})

test_that("ot_samse() without noise equals the deterministic surplus at n0", {
  # With sd_env = 0 the population does not decline as long as H is below the
  # Ricker surplus at n0: n0 * (exp(r * (1 - n0 / K)) - 1).
  d <- data.frame(rmax = 0.3, sd = 0, n0 = 20, K = 40, take = 1)
  res <- ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
                  nsims = 20, years = 40, tol = 1e-4, seed = 1)
  expected <- 20 * (exp(0.3 * (1 - 20 / 40)) - 1)
  expect_equal(res$limit, expected, tolerance = 1e-3)
})

test_that("ot_samse() is precautionary: no surplus at K, less with more noise", {
  at_k <- ot_samse(data.frame(rmax = 0.3, sd = 0, n0 = 40, K = 40, take = 0),
                   rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
                   nsims = 20, years = 40, seed = 1)
  expect_equal(at_k$limit, 0, tolerance = 0.05)
  calm <- ot_samse(data.frame(rmax = 0.3, sd = 0, n0 = 20, K = 40, take = 1),
                   rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
                   nsims = 200, years = 40, seed = 1)
  noisy <- ot_samse(data.frame(rmax = 0.3, sd = 0.3, n0 = 20, K = 40, take = 1),
                    rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
                    nsims = 200, years = 40, seed = 1)
  expect_lt(noisy$limit, calm$limit)
  # Reproducible with a seed.
  again <- ot_samse(data.frame(rmax = 0.3, sd = 0.3, n0 = 20, K = 40, take = 1),
                    rmax = rmax, sd_env = sd, removal = take, n0 = n0, k = K,
                    nsims = 200, years = 40, seed = 1)
  expect_equal(noisy$limit, again$limit)
  expect_true(noisy$p_extinct >= 0 && noisy$p_extinct <= 1)
})

# Benchmarks against Levi et al. (2011) ----------------------------------------
# The reference numbers were computed with an independent, line-by-line
# transcription of the steady-state and time-stepped implementations published
# with Levi et al. (2011, Ecological Applications 21: 1802-1818; Ecological
# Archives A021-081-S1), using the exact value of pi.

test_that("ot_biode(method = 'steady') reproduces Levi et al. (2011)", {
  st <- data.frame(x = c(10, 28, 33), y = c(12, 8, 22), pop = c(50, 90, 30))
  # Reference grid: cell (i, j) has its centre at x = j, y = i.
  run <- function(theta) {
    res <- ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07, hphy = 40,
                    kill_rate = 0.1, sigma = 6, er = 0.02, theta = theta,
                    extent = c(0.5, 40.5, 0.5, 30.5), resolution = 1)
    attr(res, "grid")$density
  }
  n1 <- run(1)
  expect_equal(sum(n1), 13182.7352771376, tolerance = 1e-10)
  expect_equal(sum(n1 == 0), 369)
  expect_equal(n1[15, 20], 11.9412480545, tolerance = 1e-9)
  expect_equal(n1[1, 1], 24.0430104391, tolerance = 1e-9)
  expect_equal(n1[26, 36], 6.2206158534, tolerance = 1e-9)
  n2 <- run(2)
  expect_equal(sum(n2), 15949.3364017728, tolerance = 1e-10)
  expect_equal(n2[15, 20], 17.2780554856, tolerance = 1e-9)
  expect_equal(n2[26, 36], 12.4705812349, tolerance = 1e-9)
})

test_that("ot_biode(method = 'dynamic') reproduces Levi et al. (2011) away from the edges", {
  # Reference 0-based cell (row, col) has its centre at y = row + 1, x = col + 1.
  st <- data.frame(x = c(23, 29), y = c(23, 28), pop = c(60, 40))
  res <- ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07, hphy = 40,
                  kill_rate = 0.1, sigma = 6, er = 0.02,
                  extent = c(0.5, 50.5, 0.5, 50.5), method = "dynamic",
                  years = 60, dispersal = 0.02)
  n <- attr(res, "grid")$density
  # The published implementation treats the grid edges differently; compare
  # the central block.
  expect_equal(sum(n[16:35, 16:35]), 2543.7089769185, tolerance = 1e-8)
  expect_equal(n[31, 21], 7.5332948498, tolerance = 1e-8)
  expect_equal(n[26, 26], 0.0004621074, tolerance = 1e-6)
})

test_that("ot_biode(method = 'dynamic') without movement converges to the steady state", {
  st <- data.frame(x = 10, y = 10, pop = 80)
  args <- list(st, x = quote(x), y = quote(y), humans = quote(pop), k = 25, r = 0.3,
               hphy = 40, kill_rate = 0.1, sigma = 4, extent = c(0, 20, 0, 20))
  steady <- attr(do.call(ot_biode, args), "grid")$density
  dynamic <- attr(do.call(ot_biode, c(args, method = "dynamic", dispersal = 0, years = 400)),
                  "grid")$density
  # Cells right at the extinction threshold converge slowly, so compare the
  # average difference.
  expect_lt(mean(abs(dynamic - steady)), 0.01)
})

test_that("ot_biode() summaries use the study area, not an arbitrary margin", {
  st <- data.frame(x = c(10, 30), y = c(15, 20), pop = c(80, 120))
  base <- function(...) ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07,
                                 hphy = 40, kill_rate = 0.1, sigma = 6, resolution = 2, ...)
  a <- base(extent = c(0, 40, 0, 36))
  b <- base(extent = c(0, 40, 0, 36), margin = 100) # margin is ignored with extent
  expect_equal(a$frac_extirpated, b$frac_extirpated)
  expect_equal(a$area, 40 * 36)
  expect_equal(a$area_extirpated, a$frac_extirpated * a$area)
  # Without a study area a warning explains the dependence on `margin`.
  rlang::reset_warning_verbosity("offtake_biode_study_area")
  expect_warning(base(), class = "offtake_no_study_area")
})

test_that("ot_biode() restricts summaries to a boundary polygon", {
  st <- data.frame(x = 10, y = 10, pop = 50)
  square <- data.frame(x = c(0, 20, 20, 0), y = c(0, 0, 20, 20))
  # No cell centre lies on the hypotenuse x + y = 20.2.
  triangle <- data.frame(x = c(0, 20.2, 0), y = c(0, 0, 20.2))
  run <- function(b) ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07,
                              hphy = 40, kill_rate = 0.1, sigma = 4, boundary = b)
  sq <- run(square)
  tr <- run(triangle)
  expect_equal(sq$area, 400)
  # Centres (i + 0.5, j + 0.5) with i + j <= 19: 20 * 21 / 2 = 210 cells.
  expect_equal(tr$area, 210)
  expect_true(any(is.na(ot_biode_surface(tr)$density)))
})

test_that("ot_biode() does not cut the hunting ground of a settlement at the study edge", {
  # A village 3 km inside the study area hunts well beyond it; thanks to the
  # buffer its catch per effort matches a run on a much larger area.
  st <- data.frame(x = 3, y = 20, pop = 60)
  run <- function(ext) ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07,
                                hphy = 40, kill_rate = 0.1, sigma = 6, extent = ext)
  tight <- ot_biode_cpue(run(c(0, 40, 0, 40)))$local_cpue
  wide <- ot_biode_cpue(run(c(-40, 40, 0, 40)))$local_cpue
  expect_equal(tight, wide, tolerance = 0.01)
})

test_that("ot_biode() returns per-settlement CPUE", {
  st <- data.frame(x = c(10, 30), y = c(15, 20), pop = c(80, 120))
  res <- ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07, hphy = 40,
                  kill_rate = 0.1, sigma = 6, extent = c(0, 40, 0, 36))
  cp <- ot_biode_cpue(res)
  expect_equal(nrow(cp), 2)
  # Expected kills per hunt cannot exceed er * kill_rate * K.
  expect_true(all(cp$local_cpue <= 0.02 * 0.1 * 25))
})
