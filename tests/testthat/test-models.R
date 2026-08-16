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

test_that("ot_pro() matches the Robinson & Redford formula and flags overharvest", {
  d <- data.frame(K = 100, lam = 1.5, F = 0.4, H = 50)
  res <- ot_pro(d, k = K, harvest = H, lambda = lam, f = F)
  # P = 0.6 * 100 * (1.5 - 1) * 0.4 = 12
  expect_equal(res$production, 12)
  expect_false(res$sustainable)          # 50 > 12
  res2 <- ot_pro(data.frame(K = 100, lam = 1.5, F = 0.4, H = 10),
                 k = K, harvest = H, lambda = lam, f = F)
  expect_true(res2$sustainable)          # 10 < 12
})

test_that("ot_pro() derives F from longevity", {
  d <- data.frame(K = 100, lam = 1.5, H = 1, life = c(4, 8, 15))
  res <- ot_pro(d, k = K, harvest = H, lambda = lam, longevity = life)
  expect_equal(res$f, c(0.6, 0.4, 0.2))
})

test_that("ot_pbr() matches Wade (1998) and derives Nmin from n/cv", {
  res <- ot_pbr(data.frame(N = 1000, rmax = 0.04, take = 5), rmax = rmax,
                removal = take, nmin = NULL, n = N, cv = 0, fr = 0.5)
  # With cv = 0, Nmin = N; PBR = 1000 * 0.5 * 0.04 * 0.5 = 10
  expect_equal(res$pbr, 10)
  expect_true(res$sustainable)           # 5 < 10
  # A positive cv must shrink Nmin.
  res_cv <- ot_pbr(data.frame(N = 1000, rmax = 0.04, take = 5), rmax = rmax,
                   removal = take, n = N, cv = 0.3, fr = 0.5)
  expect_lt(res_cv$nmin, 1000)
})

test_that("ot_msy() equals rK/4 and reports surplus at n", {
  res <- ot_msy(data.frame(r = 0.4, K = 1000, H = 80, N = 500),
                r = r, k = K, harvest = H, n = N)
  expect_equal(res$msy, 100)             # 0.4 * 1000 / 4
  expect_equal(res$surplus_at_n, 100)    # at N = K/2 surplus == MSY
  expect_true(res$sustainable)
})

test_that("ot_samse() returns a limit below the deterministic surplus and is reproducible", {
  d <- data.frame(rmax = 0.1, sd = 0.2, n0 = 500, take = 1)
  r1 <- ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0,
                 nsims = 150, years = 30, seed = 42)
  r2 <- ot_samse(d, rmax = rmax, sd_env = sd, removal = take, n0 = n0,
                 nsims = 150, years = 30, seed = 42)
  expect_equal(r1$samse_limit, r2$samse_limit)          # reproducible
  det_surplus <- 500 * (exp(0.1) - 1)
  expect_lt(r1$samse_limit, det_surplus)                # stochastic <= deterministic
  expect_true(r1$sustainable)                           # take = 1 is tiny
})

test_that("ot_biode() reproduces the Levi et al. steady-state closed form", {
  st <- data.frame(x = 0, y = 0, pop = 100)
  res <- ot_biode(st, x = x, y = y, humans = pop, k = 25, r = 0.07, hphy = 40,
                  kill_rate = 0.1, sigma = 6, er = 0.02, theta = 1,
                  resolution = 1, margin = 18)
  surf <- ot_biode_surface(res)
  # Recompute the density at a specific cell by hand and compare.
  cell <- surf[surf$x == 6 & surf$y == 0, ]
  d2 <- 6^2
  phi <- exp(-d2 / (2 * 6^2)) / (2 * pi * sqrt(d2) + 1)
  coef <- 0.02 * 0.1 * 40 / 0.07
  expected <- max(0, 25 * (1 - coef * 100 * phi))
  expect_equal(cell$density, expected, tolerance = 1e-8)
  # Density approaches K far from the settlement (a source) but never exceeds
  # it, and is depleted near the settlement (a sink).
  expect_lte(max(surf$density), 25)
  expect_gt(max(surf$density), 24.9)
  centre <- surf[surf$x == 0 & surf$y == 0, ]
  expect_lt(centre$density, max(surf$density))
})
