# Independent numerical validation (not part of the package build).
suppressMessages(library(mediaPlanCopula))
ns <- asNamespace("mediaPlanCopula")
cat("== 1. Calibración exacta de r1 y r2 con un solo bloque/dos bloques ==\n")
for (ty in c("bbd", "nbd")) {
  f <- distrib_hybrid_copula(c(V = .2), c(V = .312), 2, ty)
  cat(ty, "reach(n=2) =", f$metrics$reach_1plus, " r2 =", .312,
      " reach(n=1) =", distrib_hybrid_copula(c(V=.2), c(V=.312), 1, ty)$metrics$reach_1plus, "\n")
}
cat("\n== 2. NBD: pmf total = dnbinom exacta ==\n")
f <- distrib_hybrid_copula(c(D = .02), c(D = .039), 50, "nbd")
mg <- f$marginals
cat("size", mg$alpha, "mean", mg$mean_contacts, " max|pmf - dnbinom| =",
    max(abs(f$distribution$pmf[1:20] - dnbinom(0:19, size = mg$alpha, mu = mg$mean_contacts))), "\n")

cat("\n== 3. Precisión de E[g] (rejilla normal + spline) ==\n")
gh <- ns$.normal_grid()
for (ty in c("bbd", "nbd")) for (q in c(0.21, 0.30, 0.355)) {
  mg <- if (ty == "bbd") ns$.marginal_bbd(.2, q, 5, "x") else ns$.marginal_nbd(.2, q, 5, "x")
  g <- ns$.propensity_curve(mg)
  cat(ty, "q=", q, " E[g] =", sum(gh$w * g(gh$x)), " p = 0.2\n")
}

cat("\n== 4. Pmf QMC vs. integración determinista 2D (Gauss-Hermite anidada) ==\n")
gh_nodes <- function(n) { i <- seq_len(n - 1L); J <- matrix(0, n, n); J[cbind(i, i + 1L)] <- sqrt(i); J[cbind(i + 1L, i)] <- sqrt(i); e <- eigen(J, symmetric = TRUE); o <- order(e$values); list(x = e$values[o], w = e$vectors[1L, o]^2) }
exact2 <- function(margs, rho, kc, n = 80L) {
  gh <- gh_nodes(n); x <- gh$x; w <- gh$w
  pm <- function(mg, z, kc) {
    lat <- ns$.latent(mg, z)
    if (mg$type == "bbd") t(vapply(lat, function(l) dbinom(0:(kc-1), mg$n, l), numeric(kc)))
    else t(vapply(lat, function(l) dpois(0:(kc-1), mg$n * l), numeric(kc)))
  }
  A <- pm(margs[[1]], x, kc)
  out <- numeric(kc)
  for (a in seq_len(n)) {
    B <- pm(margs[[2]], rho * x[a] + sqrt(1 - rho^2) * x, kc)
    wb <- colSums(B * w)               # E_b[pmf_2 | a]
    out <- out + w[a] * stats::convolve(A[a, ], rev(wb), type = "open")[1:kc]
  }
  out
}
cases <- list(
  bbd_bbd = list(r1 = c(A = .2, B = .15), r2 = c(.312, .23925), n = c(5, 5), ty = c("bbd", "bbd"), d = 0.05),
  bbd_nbd = list(r1 = c(A = .2, B = .02), r2 = c(.312, .039), n = c(5, 50), ty = c("bbd", "nbd"), d = 0.006)
)
for (nmc in names(cases)) {
  cs <- cases[[nmc]]
  D <- matrix(c(cs$r1[1], cs$d, cs$d, cs$r1[2]), 2)
  kk <- if (all(cs$ty == "bbd")) sum(cs$n) + 1L else 26L      # columnas exactas
  open <- any(cs$ty == "nbd")
  fit <- suppressWarnings(distrib_hybrid_copula(cs$r1, cs$r2, cs$n, cs$ty, D, k_max = if (open) 25 else NULL, n_qmc = 4096))
  rho <- fit$Sigma_latent[1, 2]
  margs <- lapply(1:2, function(i) if (cs$ty[i] == "bbd") ns$.marginal_bbd(cs$r1[i], cs$r2[i], cs$n[i], "x") else ns$.marginal_nbd(cs$r1[i], cs$r2[i], cs$n[i], "x"))
  ex <- exact2(margs, rho, if (open) 25 else kk)
  if (open) ex <- c(ex, 1 - sum(ex))
  cat(nmc, ": rho =", round(rho, 4), " max|pmf_QMC - pmf_exacta| (pp) =",
      round(100 * max(abs(fit$distribution$pmf - ex)), 4),
      " | err dup =", signif(fit$diagnostics$dup_error, 3), " err r1 =", signif(fit$diagnostics$reach1_error, 3), "\n")
  fit_hi <- suppressWarnings(distrib_hybrid_copula(cs$r1, cs$r2, cs$n, cs$ty, D, k_max = if (open) 25 else NULL, n_qmc = 32768))
  cat("   con n_qmc=32768: max|pmf - exacta| (pp) =", round(100 * max(abs(fit_hi$distribution$pmf - ex)), 4), "\n")
}

cat("\n== 5. Caso del usuario (n2 = 200): ¿rho calibrado? ==\n")
mu <- 1.2; rd <- .8; n2 <- 200; p2 <- mu / n2; r11 <- p2^2 * (1 + 1 / rd); r2d <- 2 * p2 - r11
D <- matrix(c(.2, 1.35 * .2 * p2, 1.35 * .2 * p2, p2), 2)
fit <- distrib_hybrid_copula(c(TV = .2, Dig = p2), c(.312, r2d), c(5, n2), c("bbd", "nbd"), D, k_max = 15)
print(round(fit$Sigma_latent, 4)); print(fit$diagnostics$pairs)

cat("\n== 6. Estado del RNG intacto ==\n")
set.seed(1); a <- runif(1); set.seed(1); invisible(distrib_gcbbd(c(.2,.15), c(.312,.23925), c(5,5), n_qmc = 512)); b <- runif(1); cat("identico:", identical(a, b), "\n")
cat("\n== 7. Tiempos ==\n")
print(system.time(distrib_hybrid_copula(c(TV = .2, Dig = p2), c(.312, r2d), c(5, n2), c("bbd", "nbd"), D)))
