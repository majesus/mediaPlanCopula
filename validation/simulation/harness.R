# Simulation harness: copula model vs. existing exposure-distribution models.
suppressMessages({library(mediaPlanR); library(mediaPlanCopula)})
NPOP <- 20000L

# ---- exact distribution of total contacts: mean of Poisson-binomial pmfs ----
pb_pmf <- function(Pins) {
  n <- nrow(Pins); N <- ncol(Pins)
  pm <- matrix(0, n, N + 1L); pm[, 1] <- 1
  for (j in seq_len(N)) {
    p <- Pins[, j]
    nw <- pm * (1 - p)
    nw[, 2:(N + 1L)] <- nw[, 2:(N + 1L)] + pm[, 1:N] * p
    pm <- nw
  }
  colMeans(pm)
}

equi <- function(m, rho) (1 - rho) * diag(m) + rho * matrix(1, m, m)
beta_par <- function(r1, th) { conc <- (1 - th) / th; list(a = r1 * conc, b = (1 - r1) * conc) }

# logit-normal location giving mean r1 for a given scale
logitnorm_mu <- function(r1, sig) {
  zz <- qnorm((1:400 - .5) / 400)
  uniroot(function(mu) mean(plogis(mu + sig * zz)) - r1, c(-12, 8))$root
}

# ---- Part A: BBD-only plans ----
gen_A <- function(dgp) {
  m <- sample(2:4, 1); ni <- sample(2:6, m, TRUE)
  r1 <- runif(m, .06, .35); th <- runif(m, .15, .6); rho <- runif(1, .1, .7)
  decay <- 1
  Z <- matrix(rnorm(NPOP * m), NPOP, m)
  if (dgp == "G6_negative") {
    a <- runif(1, .2, .6); v <- sample(c(-1, 1), m, TRUE); if (all(v == v[1])) v[1] <- -v[1]
    Sig <- (1 - a) * diag(m) + a * outer(v, v)
  } else Sig <- equi(m, rho)
  Zc <- Z %*% chol(Sig)
  P <- matrix(0, NPOP, m)
  if (dgp %in% c("G1_gauss_beta", "G5_decay", "G6_negative")) {
    for (i in 1:m) { bp <- beta_par(r1[i], th[i]); P[, i] <- qbeta(pnorm(Zc[, i]), bp$a, bp$b) }
    if (dgp == "G5_decay") decay <- runif(1, .75, .95)
  } else if (dgp == "G2_gauss_logitnormal") {
    for (i in 1:m) { sg <- runif(1, .8, 2); P[, i] <- plogis(logitnorm_mu(r1[i], sg) + sg * Zc[, i]) }
  } else if (dgp == "G3_t_copula") {
    nu <- 4; Zt <- Zc / sqrt(rchisq(NPOP, nu) / nu)
    for (i in 1:m) { bp <- beta_par(r1[i], th[i]); P[, i] <- qbeta(pt(Zt[, i], nu), bp$a, bp$b) }
  } else if (dgp == "G4_latent_class") {
    w <- c(.25, .35, .40); base <- c(2.5, 1, .25); cls <- sample(1:3, NPOP, TRUE, prob = w)
    for (i in 1:m) {
      mult <- base * exp(.3 * rnorm(3)); mult <- mult / sum(w * mult)
      P[, i] <- pmin(r1[i] * mult[cls], .97)
    }
  } else stop("dgp?")
  P <- pmin(pmax(P, 1e-9), 1 - 1e-9)
  veh <- rep(seq_len(m), ni); pos <- unlist(lapply(ni, seq_len))
  Pins <- P[, veh, drop = FALSE] * matrix(rep(decay^(pos - 1), each = NPOP), NPOP, length(pos))
  list(m = m, ni = ni, P = P, decay = decay, Pins = Pins, dgp = dgp)
}

truth_A <- function(sc) pb_pmf(sc$Pins)

inputs_A <- function(sc, npanel = Inf) {
  P <- sc$P; m <- sc$m
  if (is.infinite(npanel)) {
    r1 <- colMeans(P); r2 <- colMeans(1 - (1 - P) * (1 - P * sc$decay)); D <- crossprod(P) / nrow(P)
  } else {
    idx <- sample.int(nrow(P), npanel)
    X1 <- matrix(rbinom(npanel * m, 1, P[idx, ]), npanel, m)
    X2 <- matrix(rbinom(npanel * m, 1, P[idx, ] * sc$decay), npanel, m)
    r1 <- colMeans(X1); r2 <- colMeans(pmax(X1, X2)); D <- crossprod(X1) / npanel
  }
  diag(D) <- r1
  list(r1 = r1, r2 = r2, D = D)
}

# ---- Part B: hybrid plans (BBD + NBD) ----
gen_B <- function(dgp) {
  n1 <- sample(3:8, 1); r1_1 <- runif(1, .1, .35); th <- runif(1, .15, .6)
  n2 <- sample(c(30, 60, 100), 1); s_blk <- runif(1, .006, .03); shp <- runif(1, .4, 2)
  rho <- runif(1, .1, .7)
  if (dgp == "H1_gauss_beta_gamma") {
    Z <- matrix(rnorm(NPOP * 2), NPOP, 2) %*% chol(equi(2, rho))
    bp <- beta_par(r1_1, th)
    P1 <- qbeta(pnorm(Z[, 1]), bp$a, bp$b)
    lam <- qgamma(pnorm(Z[, 2]), shape = shp, rate = shp / s_blk)
  } else if (dgp == "H2_latent_class") {
    w <- c(.25, .35, .40); cls <- sample(1:3, NPOP, TRUE, prob = w)
    m1 <- c(2.5, 1, .25) * exp(.3 * rnorm(3)); m1 <- m1 / sum(w * m1)
    m2 <- c(3, 1, .2) * exp(.3 * rnorm(3)); m2 <- m2 / sum(w * m2)
    P1 <- pmin(r1_1 * m1[cls], .97)
    lam <- s_blk * m2[cls] * rgamma(NPOP, 3, 3)
  } else stop("dgp?")
  P1 <- pmin(pmax(P1, 1e-9), 1 - 1e-9); lam <- pmax(lam, 1e-12)
  list(n1 = n1, n2 = n2, P1 = P1, lam = lam, dgp = dgp)
}

truth_B <- function(sc) {
  K <- sc$n1 + sc$n2; n1 <- sc$n1
  bin <- matrix(dbinom(rep(0:n1, each = NPOP), n1, rep(sc$P1, n1 + 1L)), NPOP, n1 + 1L)
  poi <- matrix(dpois(rep(0:(K - 1L), each = NPOP), rep(sc$n2 * sc$lam, K)), NPOP, K)
  out <- matrix(0, NPOP, K)
  for (a in 0:n1) { len <- K - a; cols <- (a + 1L):K; out[, cols] <- out[, cols] + bin[, a + 1L] * poi[, 1:len, drop = FALSE] }
  pm <- colMeans(out); c(pm, max(0, 1 - sum(pm)))
}

inputs_B <- function(sc, npanel = Inf) {
  pb <- 1 - exp(-sc$lam); P1 <- sc$P1
  if (is.infinite(npanel)) {
    r1 <- c(mean(P1), mean(pb)); r2 <- c(mean(1 - (1 - P1)^2), mean(1 - exp(-2 * sc$lam)))
    d12 <- mean(P1 * pb)
  } else {
    idx <- sample.int(NPOP, npanel)
    x11 <- rbinom(npanel, 1, P1[idx]); x12 <- rbinom(npanel, 1, P1[idx])
    x21 <- rbinom(npanel, 1, pb[idx]); x22 <- rbinom(npanel, 1, pb[idx])
    r1 <- c(mean(x11), mean(x21)); r2 <- c(mean(pmax(x11, x12)), mean(pmax(x21, x22))); d12 <- mean(x11 * x21)
  }
  D <- matrix(c(r1[1], d12, d12, r1[2]), 2)
  list(r1 = r1, r2 = r2, D = D)
}

# ---- Model runner (returns pmf over 0..N, N = total insertions/blocks) ----
safe <- function(expr) {
  nw <- 0L
  res <- tryCatch(withCallingHandlers(expr, warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") }),
                  error = function(e) structure(list(msg = conditionMessage(e)), class = "failure"))
  list(res = res, nw = nw)
}

run_models <- function(inp, ni, types, Kopen = NULL) {
  N <- sum(ni)
  vd <- data.frame(insertions = ni, R1 = inp$r1, R2 = inp$r2)
  M <- inp$D; diag(M) <- 2 * inp$r1 - inp$r2
  ex <- function(f) { s <- safe(f()); if (inherits(s$res, "failure")) list(pmf = NULL, nw = s$nw) else list(pmf = as.numeric(s$res), nw = s$nw) }
  out <- list()
  cop_call <- function(dup) {
    function() distrib_hybrid_copula(inp$r1, inp$r2, ni, types, dup_matrix = dup, k_max = Kopen, n_qmc = 4096)$distribution$pmf
  }
  out$copula <- ex(cop_call(inp$D))
  out$indep  <- ex(cop_call(NULL))
  out$sainsbury <- ex(function() { x <- calc_sainsbury(inp$r1, 1, ni); c(1 - x$reach$percent / 100, x$distribution$percent / 100) })
  out$metheringham <- ex(function() { x <- calc_metheringham(inp$r1, ni, M, 1); c(1 - x$reach$percent / 100, x$distribution$percent / 100) })
  out$canex <- ex(function() calc_canex(vd, inp$D)$distribution$probability)
  out$msad  <- ex(function() calc_msad(vd, inp$D)$distribution$probability)
  out$cbd   <- ex(function() calc_cbd(vd, inp$D)$distribution$probability)
  out$mbd   <- ex(function() calc_mbd(vd, inp$D)$distribution$probability)
  out$csd   <- ex(function() calc_csd(vd, inp$D)$distribution$probability)
  out
}

score <- function(pmf, truth) {
  if (is.null(pmf) || length(pmf) != length(truth) || anyNA(pmf)) return(c(tv = NA, e1 = NA, e3 = NA))
  pmf <- pmax(pmf, 0); pmf <- pmf / sum(pmf)
  c(tv = 0.5 * sum(abs(pmf - truth)),
    e1 = 100 * ((1 - pmf[1]) - (1 - truth[1])),
    e3 = 100 * ((1 - sum(pmf[1:3])) - (1 - sum(truth[1:3]))))
}

one_scenario <- function(id, part, dgp, npanel_levels, seed0) {
  set.seed(seed0 + id)
  if (part == "A") { sc <- gen_A(dgp); tr <- truth_A(sc); types <- rep("bbd", sc$m); ni <- sc$ni; Kopen <- NULL; inf <- inputs_A }
  else { sc <- gen_B(dgp); tr <- truth_B(sc); types <- c("bbd", "nbd"); ni <- c(sc$n1, sc$n2); Kopen <- sc$n1 + sc$n2; inf <- inputs_B }
  rows <- list()
  for (npn in npanel_levels) {
    inp <- inf(sc, npn)
    t0 <- proc.time()[[3]]
    ms <- run_models(inp, ni, types, Kopen)
    for (nm in names(ms)) {
      sc_vec <- score(ms[[nm]]$pmf, tr)
      rows[[length(rows) + 1L]] <- data.frame(id = id, part = part, dgp = dgp, regime = if (is.infinite(npn)) "population" else paste0("panel", npn),
        model = nm, ok = !is.null(ms[[nm]]$pmf) && !is.na(sc_vec[["tv"]]), nwarn = ms[[nm]]$nw,
        tv = sc_vec[["tv"]], e1 = sc_vec[["e1"]], e3 = sc_vec[["e3"]], m = length(ni), N = sum(ni),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
