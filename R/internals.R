# Internal helpers (not exported).
# Notation. For vehicle i, p = R1 is the reach of one insertion (or one unit
# block), q = R2 the reach of two, n the number of insertions (or blocks).

# --- Marginal calibration ---------------------------------------------------

# Beta-Binomial: X | P ~ Bin(n, P), P ~ Beta(a, b), E[P] = p.
# R2 = 2p - E[P^2] identifies the concentration a + b exactly.
.marginal_bbd <- function(p, q, n, name) {
  if (!(q > p && q < 2 * p - p^2)) {
    stop(sprintf(paste0("Vehicle '%s': r2 must lie strictly between r1 and ",
                        "2*r1 - r1^2 (in proportions: (%.6g, %.6g)); got %.6g."),
                 name, p, 2 * p - p^2, q), call. = FALSE)
  }
  r11 <- 2 * p - q                         # P(exposed to both insertions)
  theta <- (r11 - p^2) / (p * (1 - p))     # intra-vehicle correlation
  conc <- (1 - theta) / theta              # a + b
  list(type = "bbd", name = name, p = p, q = q, n = n,
       alpha = p * conc, beta = (1 - p) * conc, mean_contacts = n * p)
}

# Negative Binomial as a Poisson-Gamma mixture: X | L ~ Poisson(L),
# L = n * G, G ~ Gamma(shape = r, rate = r / s), so that a unit block has
# reach 1 - (1 + s / r)^(-r) and two blocks 1 - (1 + 2 s / r)^(-r).
# Both reaches are matched exactly by solving for (r, s).
.marginal_nbd <- function(p, q, n, name) {
  if (!(q > p && q < 2 * p - p^2)) {
    stop(sprintf(paste0("Vehicle '%s': r2 must lie strictly between r1 and ",
                        "2*r1 - r1^2 (in proportions: (%.6g, %.6g)); got %.6g."),
                 name, p, 2 * p - p^2, q), call. = FALSE)
  }
  target <- log1p(-p) - log1p(-q)          # log((1 - p) / (1 - q)) > 0
  f <- function(log_r) {
    r <- exp(log_r)
    big_l <- -log1p(-p) / r                # log of (1 - p)^(-1 / r)
    r * log1p(-expm1(-big_l)) - target     # increasing in r
  }
  lo <- log(5e-3)
  hi <- log(1e6)
  if (!(f(lo) < 0 && f(hi) > 0)) {
    stop(sprintf(paste0("Vehicle '%s': no Negative Binomial dispersion in ",
                        "[0.005, 1e6] reproduces r1 and r2."), name),
         call. = FALSE)
  }
  r <- exp(stats::uniroot(f, c(lo, hi), tol = 1e-14)$root)
  s <- r * expm1(-log1p(-p) / r)           # mean contacts per unit block
  list(type = "nbd", name = name, p = p, q = q, n = n,
       alpha = r, beta = r / s, mean_contacts = n * s)
}

# --- Latent transforms ------------------------------------------------------

# Latent per-person quantity as a function of a standard normal z:
# BBD -> propensity P; NBD -> intensity per unit block (Gamma).
.latent <- function(mg, z) {
  u <- pnorm(z, lower.tail = FALSE)        # upper-tail form keeps precision
  if (mg$type == "bbd") {
    qbeta(u, mg$alpha, mg$beta, lower.tail = FALSE)
  } else {
    qgamma(u, shape = mg$alpha, rate = mg$beta, lower.tail = FALSE)
  }
}

# Probability of being exposed to one unit block, given the latent quantity.
.propensity <- function(mg, v) {
  if (mg$type == "bbd") v else -expm1(-v)
}

# --- Quadrature -------------------------------------------------------------

# Equispaced grid on [-8, 8] with standard normal weights (trapezoid rule,
# spectrally accurate for these smooth, rapidly decaying integrands; weights
# are normalised to sum to one).
.normal_grid <- function(by = 0.05, lim = 8) {
  x <- seq(-lim, lim, by = by)
  w <- dnorm(x)
  list(x = x, w = w / sum(w))
}

# Monotone interpolant of the propensity as a function of z.
.propensity_curve <- function(mg) {
  zg <- seq(-10, 10, by = 0.02)
  f <- splinefun(zg, .propensity(mg, .latent(mg, zg)), method = "monoH.FC")
  function(z) f(pmax(pmin(z, 10), -10))
}

# Duplication D(rho) = E[g1(Z1) g2(rho Z1 + sqrt(1 - rho^2) Z2)] for a pair.
.dup_function <- function(g1, g2, grid) {
  x <- grid$x
  w <- grid$w
  a <- g1(x)
  function(rho) {
    z2 <- outer(rho * x, sqrt(1 - rho^2) * x, "+")
    sum(w * a * (matrix(g2(z2), nrow(z2)) %*% w))
  }
}

# Latent correlation reproducing the target duplication. Returns the root, or
# the nearest attainable bound with a status flag; never changes silently.
.calibrate_rho <- function(g1, g2, target, grid, bound = 0.9999) {
  d_fun <- .dup_function(g1, g2, grid)
  f <- function(rho) d_fun(rho) - target
  f_lo <- f(-bound)
  f_hi <- f(bound)
  if (f_lo >= 0) {
    return(list(rho = -bound, status = "below_attainable", attained = d_fun(-bound)))
  }
  if (f_hi <= 0) {
    return(list(rho = bound, status = "above_attainable", attained = d_fun(bound)))
  }
  rho <- uniroot(f, c(-bound, bound), tol = 1e-12)$root
  list(rho = rho, status = "ok", attained = d_fun(rho))
}

# --- Quasi-Monte Carlo ------------------------------------------------------

.radical_inverse <- function(i, base) {
  out <- numeric(length(i))
  f <- 1 / base
  k <- i
  while (any(k > 0)) {
    out <- out + f * (k %% base)
    k <- k %/% base
    f <- f / base
  }
  out
}

# Deterministic shifted Halton points in (0, 1)^dim (no use of the RNG).
.halton <- function(n, dim) {
  primes <- c(2L, 3L, 5L, 7L, 11L, 13L, 17L, 19L, 23L, 29L, 31L, 37L)
  if (dim > length(primes)) {
    stop("At most 12 vehicles are supported with more than one vehicle.",
         call. = FALSE)
  }
  idx <- seq_len(n)
  shift <- (seq_len(dim) * 0.6180339887498949) %% 1
  u <- vapply(seq_len(dim), function(j) {
    (.radical_inverse(idx, primes[j]) + shift[j]) %% 1
  }, numeric(n))
  u <- matrix(u, nrow = n, ncol = dim)
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}

# Convolution, across the rows of two matrices, truncated to `kc` columns.
.row_convolve <- function(a, b, kc) {
  out <- matrix(0, nrow(a), kc)
  for (j in seq_len(min(ncol(a), kc))) {
    len <- min(ncol(b), kc - j + 1L)
    cols <- j:(j + len - 1L)
    out[, cols] <- out[, cols] + a[, j] * b[, seq_len(len), drop = FALSE]
  }
  out
}
