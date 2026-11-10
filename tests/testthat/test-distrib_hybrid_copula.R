r1 <- c(A = 0.20, B = 0.15)
r2 <- c(A = 0.312, B = 0.23925)
d_ok <- matrix(c(0.20, 0.05, 0.05, 0.15), 2)

test_that("one and two unit blocks reproduce r1 and r2 exactly (both families)", {
  for (ty in c("bbd", "nbd")) {
    f1 <- distrib_hybrid_copula(c(V = 0.2), c(V = 0.312), 1, ty)
    f2 <- distrib_hybrid_copula(c(V = 0.2), c(V = 0.312), 2, ty)
    expect_equal(f1$metrics$reach_1plus, 0.2, tolerance = 1e-10)
    expect_equal(f2$metrics$reach_1plus, 0.312, tolerance = 1e-10)
  }
})

test_that("NBD marginal equals the exact Negative Binomial", {
  f <- distrib_hybrid_copula(c(D = 0.02), c(D = 0.039), 50, "nbd")
  g <- f$marginals
  expect_equal(f$distribution$pmf[1:15],
               dnbinom(0:14, size = g$alpha, mu = g$mean_contacts),
               tolerance = 1e-12)
  expect_true(f$metrics$open_last_bin)
  expect_equal(sum(f$distribution$pmf), 1)
})

test_that("BBD marginal equals the exact Beta-Binomial", {
  f <- distrib_hybrid_copula(c(V = 0.2), c(V = 0.312), 5, "bbd")
  a <- f$marginals$alpha
  b <- f$marginals$beta
  k <- 0:5
  ref <- choose(5, k) * beta(k + a, 5 - k + b) / beta(a, b)
  expect_equal(f$distribution$pmf, ref, tolerance = 1e-12)
  expect_false(f$metrics$open_last_bin)
})

test_that("independent vehicles give the convolution of the marginals", {
  f <- distrib_gcbbd(r1, r2, c(5, 5), n_qmc = 4096)
  m1 <- distrib_gcbbd(r1[1], r2[1], 5)$distribution$pmf
  m2 <- distrib_gcbbd(r1[2], r2[2], 5)$distribution$pmf
  ref <- as.numeric(tapply(outer(m1, m2), outer(0:5, 0:5, "+"), sum))
  expect_equal(f$distribution$pmf, ref, tolerance = 2e-3)
  expect_equal(f$Sigma_latent[1, 2], 0, tolerance = 1e-3)
})

test_that("GRPs, reach of one block and duplication are conserved", {
  f <- distrib_hybrid_copula(
    c(TV = 0.2, Dig = 0.02), c(0.312, 0.039), c(5, 50), c("bbd", "nbd"),
    dup_matrix = matrix(c(0.2, 0.006, 0.006, 0.02), 2), n_qmc = 4096)
  expect_equal(f$metrics$grp,
               100 * (5 * 0.2 + f$marginals$mean_contacts[2]), tolerance = 1e-12)
  expect_lt(f$diagnostics$reach1_error, 1e-3)
  expect_lt(f$diagnostics$dup_error, 1e-4)
  expect_equal(unname(f$dup_fitted[1, 2]), 0.006, tolerance = 2e-2)
  expect_equal(sum(f$distribution$pmf), 1)
  expect_true(all(f$distribution$pmf >= 0))
  expect_equal(f$distribution$cdf[nrow(f$distribution)], 1)
})

test_that("percent and proportion scales are equivalent", {
  d_pct <- matrix(c(20, 5, 5, 15), 2)
  a <- distrib_gcbbd(r1, r2, c(5, 5), d_ok, n_qmc = 1024)
  b <- distrib_gcbbd(100 * r1, 100 * r2, c(5, 5), d_pct, scale = "percent",
                     n_qmc = 1024)
  expect_equal(b$distribution$pmf, 100 * a$distribution$pmf, tolerance = 1e-10)
  expect_equal(b$metrics$reach_1plus, 100 * a$metrics$reach_1plus,
               tolerance = 1e-10)
})

test_that("dup_type = 'reach' is equivalent to duplication", {
  reach <- outer(r1, r1, "+") - d_ok
  a <- distrib_gcbbd(r1, r2, c(5, 5), d_ok, n_qmc = 1024)
  b <- distrib_gcbbd(r1, r2, c(5, 5), reach, dup_type = "reach", n_qmc = 1024)
  expect_equal(a$Sigma_latent, b$Sigma_latent, tolerance = 1e-10)
})

test_that("an unattainable duplication warns and uses the nearest bound", {
  d_high <- matrix(c(0.20, 0.0908, 0.0908, 0.15), 2)
  expect_warning(f <- distrib_gcbbd(r1, r2, c(5, 5), d_high, n_qmc = 512),
                 "not attainable")
  expect_identical(f$diagnostics$pairs$status, "above_attainable")
  expect_lt(f$diagnostics$pairs$attained, 0.0908)
})

test_that("a duplication impossible for the reaches is an error", {
  d_bad <- matrix(c(0.20, 0.30, 0.30, 0.15), 2)
  expect_error(distrib_gcbbd(r1, r2, c(5, 5), d_bad), "impossible")
})

test_that("input errors are informative", {
  expect_error(distrib_gcbbd(c(20, 15), c(31.2, 23.9), c(5, 5)), "percent")
  expect_error(distrib_gcbbd(r1, c(0.2, 0.2), c(5, 5)), "strictly between")
  expect_error(distrib_gcbbd(r1, r2, c(5, 5.5)), "integer valued")
  expect_error(distrib_gcbbd(r1, r2, c(5, 5), matrix(1, 3, 3)), "m x m")
  expect_error(distrib_gcbbd(r1, r2, c(5, 5), matrix(c(.2, .05, .06, .15), 2)),
               "symmetric")
  expect_error(distrib_hybrid_copula(r1, r2, c(5, 5), "poisson"), "model_type")
})

test_that("the global random state is untouched", {
  set.seed(42)
  a <- runif(1)
  set.seed(42)
  invisible(distrib_gcbbd(r1, r2, c(5, 5), d_ok, n_qmc = 512))
  expect_identical(runif(1), a)
})

test_that("results are deterministic", {
  a <- distrib_gcbbd(r1, r2, c(5, 5), d_ok, n_qmc = 512)
  b <- distrib_gcbbd(r1, r2, c(5, 5), d_ok, n_qmc = 512)
  expect_identical(a$distribution, b$distribution)
})

test_that("k_max truncates with an open last bin that keeps the mass", {
  f <- distrib_hybrid_copula(
    c(TV = 0.2, Dig = 0.02), c(0.312, 0.039), c(5, 50), c("bbd", "nbd"),
    dup_matrix = matrix(c(0.2, 0.006, 0.006, 0.02), 2), k_max = 6, n_qmc = 1024)
  expect_equal(nrow(f$distribution), 7L)
  expect_true(f$metrics$open_last_bin)
  expect_equal(sum(f$distribution$pmf), 1)
  g <- distrib_gcbbd(r1, r2, c(5, 5), d_ok, k_max = 4, n_qmc = 1024)
  expect_equal(nrow(g$distribution), 5L)
  expect_true(g$metrics$open_last_bin)
})

test_that("three vehicles run and print", {
  d3 <- matrix(c(0.20, 0.045, 0.03, 0.045, 0.15, 0.02, 0.03, 0.02, 0.10), 3)
  expect_warning(
    f <- distrib_gcbbd(c(0.2, 0.15, 0.1), c(0.312, 0.23925, 0.1855), c(4, 4, 4),
                       d3, n_qmc = 1024),
    NA)
  expect_equal(dim(f$Sigma_latent), c(3L, 3L))
  expect_equal(nrow(f$distribution), 13L)
  expect_lt(f$diagnostics$dup_error, 1e-3)
  expect_false(f$diagnostics$psd_adjusted)
  expect_output(print(f), "GRPs")
})

test_that("an inconsistent set of duplications triggers the PSD warning", {
  d_inc <- matrix(c(.2, .065, .025, .065, .2, .065, .025, .065, .2), 3)
  msgs <- character()
  f <- withCallingHandlers(
    distrib_gcbbd(rep(0.2, 3), rep(0.312, 3), c(3, 3, 3), d_inc, n_qmc = 512),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("positive definite", msgs)))
  expect_true(f$diagnostics$psd_adjusted)
})

test_that("a strongly polarized marginal triggers the quadrature warning", {
  r2_pol <- c(A = 0.20 + 1e-3 * 0.20 * 0.80, B = 0.23925)
  expect_warning(distrib_gcbbd(r1, r2_pol, c(5, 5), n_qmc = 256),
                 "strongly polarized")
  f <- distrib_gcbbd(r1, r2, c(5, 5), n_qmc = 256)
  expect_lt(max(f$diagnostics$quadrature_error), 1e-3)
})
