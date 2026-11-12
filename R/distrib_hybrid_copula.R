#' Exposure distribution of a mixed media plan with a Gaussian copula
#'
#' Computes the distribution of the total number of contacts of a media plan
#' whose vehicles are linked by a Gaussian copula. Each vehicle has either a
#' bounded number of insertions (Beta-Binomial marginal, `"bbd"`) or an
#' unbounded stream of impressions (Negative Binomial marginal, `"nbd"`).
#' `distrib_gcbbd()` is the special case in which every vehicle is `"bbd"`.
#'
#' @details
#' **Marginals.** For a `"bbd"` vehicle the number of contacts is
#' `X | P ~ Binomial(n, P)` with `P ~ Beta(alpha, beta)`; `alpha + beta` is
#' obtained exactly from `r1` and `r2`, because the probability of being
#' exposed to two insertions is `2 * r1 - r2 = E[P^2]`. For an `"nbd"` vehicle
#' the number of impressions is `X | L ~ Poisson(L)` with `L` Gamma distributed
#' over `n` unit blocks; the shape and the mean per block are obtained by
#' solving exactly the two equations `r1 = 1 - (1 + s / r)^(-r)` and
#' `r2 = 1 - (1 + 2 s / r)^(-r)`. Hence, for both families, the reach of one
#' and of two unit blocks is reproduced exactly, and the expected number of
#' contacts is `n * r1` for `"bbd"` and `n * s` for `"nbd"` (where `s >= r1`
#' because one block can deliver several impressions to the same person).
#' Contacts are therefore exposures for `"bbd"` vehicles and impressions for
#' `"nbd"` vehicles, and the gross rating points (GRPs) are `100` times the
#' expected number of contacts.
#'
#' **Dependence.** Given a latent Gaussian vector with correlation matrix
#' `Sigma`, vehicles are conditionally independent; the latent quantity of each
#' vehicle is mapped to its propensity to be exposed to one unit block
#' (`P` for `"bbd"`, `1 - exp(-L / n)` for `"nbd"`). The duplication of the
#' first unit blocks of vehicles `i` and `j` is then
#' `D_ij = E[pi_i * pi_j]`, which is a strictly increasing function of
#' `Sigma[i, j]`. Each off-diagonal element is calibrated, one pair at a time,
#' by solving `D_ij(Sigma[i, j]) = target` with quadrature on a normal grid.
#' Pairwise calibration does not guarantee a positive definite `Sigma` when
#' there are three or more vehicles; if it is not, the eigenvalues are floored
#' and a warning is issued.
#'
#' **Integration.** The distribution is the average, over a deterministic
#' antithetic Halton sample of the latent vector (moment-matched to a unit
#' covariance), of the exact convolution of the conditional marginals. No random
#' number generator is used, so results are reproducible and the global random
#' state is untouched. Increase `n_qmc` to assess the numerical error, and
#' check `diagnostics`, which reports the error of the fitted reach of one
#' block and of the fitted duplications. With more than six vehicles a larger
#' `n_qmc` is advisable; at most 12 vehicles are supported.
#'
#' **Truncation.** If any vehicle is `"nbd"`, the support is unbounded. The
#' distribution is reported for `0, ..., k_max`, and the last category is
#' "`k_max` or more" (`open_last_bin = TRUE`). The default `k_max` is the sum of
#' the number of insertions of the `"bbd"` vehicles and of the `1 - 1e-6`
#' quantiles of the `"nbd"` vehicles. For a plan with only `"bbd"` vehicles the
#' support is `0, ..., sum(insertions)` and is exact unless a smaller `k_max`
#' is requested.
#'
#' **Scale.** Probabilities are entered and returned in the unit given by
#' `scale`. The scale is never guessed: with `scale = "proportion"` any value
#' above one is an error.
#'
#' **Limitations.** The reach of one and of two insertions of each vehicle is
#' reproduced exactly, and the expected number of contacts follows from them.
#' The duplication of the first insertions is reproduced only approximately,
#' because the correlations are calibrated one pair at a time; when the
#' calibrated matrix is not positive definite it is adjusted, and the size of the
#' deviation is reported in `diagnostics$dup_error`. This happens more often as
#' the number of vehicles grows and when the duplications come from observed
#' data. The reach of one and two insertions and the duplication do not
#' determine the distribution of contacts: the shape of the distribution beyond
#' them comes from the assumptions of the model (Beta or Gamma propensities, a
#' Gaussian copula, and exchangeable insertions within a vehicle), and models
#' with the same inputs and other assumptions give other distributions. The
#' accuracy of the model relative to other multivehicle models depends on the
#' mechanism that generates the data: in simulations its advantage varied across
#' mechanisms, and on survey-derived data the differences were of a few
#' thousandths of total variation distance and, out of sample, not
#' distinguishable from sampling noise. The choice of model should therefore
#' rest on its assumptions and on the diagnostics, not on a presumption of
#' higher accuracy.
#'
#' @param r1 Numeric vector with the reach of one insertion (one unit block)
#'   of each vehicle.
#' @param r2 Numeric vector with the reach of two insertions (two unit blocks).
#'   For every vehicle it must satisfy `r1 < r2 < 2 * r1 - r1^2`.
#' @param insertions Positive numeric vector with the number of insertions (or
#'   unit blocks). It must be integer valued for `"bbd"` vehicles.
#' @param model_type `"bbd"` or `"nbd"`, one per vehicle or a single value
#'   recycled.
#' @param dup_matrix Optional symmetric `m x m` matrix. With
#'   `dup_type = "duplication"` the element `[i, j]` is the proportion exposed
#'   to the first insertion of both vehicles; with `dup_type = "reach"` it is
#'   the combined reach of the first insertion of both vehicles. The diagonal
#'   is ignored. If `NULL`, vehicles are independent.
#' @param dup_type `"duplication"` (default) or `"reach"`.
#' @param scale `"proportion"` (default) or `"percent"`; unit of `r1`, `r2`,
#'   `dup_matrix` and of the returned probabilities.
#' @param k_max Optional integer, the largest number of contacts reported (see
#'   Details).
#' @param n_qmc Number of quasi-Monte Carlo points (default 4096).
#'
#' @return An object of class `mp_copula_exposure`, a list with
#' \describe{
#'   \item{`distribution`}{Data frame with `k`, `pmf`, `cdf` and
#'     `reach_k_plus` (probability of `k` or more contacts), in the unit of
#'     `scale`.}
#'   \item{`metrics`}{List with `grp`, `reach_1plus` and `reach_3plus` (in the
#'     unit of `scale`, except `grp`), `ots` (average opportunities to see
#'     among those reached, `NA` if `reach_1plus` is zero), `mean_contacts`,
#'     `variance` (`NA` if the last bin is open), `max_contacts` and
#'     `open_last_bin`.}
#'   \item{`marginals`}{Data frame with the calibrated marginal of each vehicle
#'     (`alpha`, `beta`: shapes of the Beta for `"bbd"`; shape and rate of the
#'     Gamma intensity per block for `"nbd"`).}
#'   \item{`Sigma_latent`}{Calibrated latent correlation matrix.}
#'   \item{`dup_observed`, `dup_fitted`}{Target and fitted duplication matrices.}
#'   \item{`diagnostics`}{List with `pairs` (target, attained value and
#'     status of each calibration), `psd_adjusted`, `reach1_error` and
#'     `dup_error` (largest absolute errors, in proportions) and
#'     `quadrature_error` (relative error of the quadrature of `r1` for each
#'     vehicle; a warning is issued above 1e-3, which happens only for very
#'     polarized marginals with `r2` extremely close to `r1`).}
#'   \item{`scale`}{The scale used.}
#' }
#'
#' @references
#' Danaher, P. J. and Smith, M. S. (2011). Modeling multivariate distributions
#' using copulas: Applications in marketing. *Marketing Science*, 30(1), 4-21.
#' \doi{10.1287/mksc.1090.0491}
#'
#' Danaher, P. J. (1991). A canonical expansion model for multivariate media
#' exposure distributions: A generalization of the "duplication of viewing
#' law". *Journal of Marketing Research*, 28(3), 361-367.
#' \doi{10.1177/002224379102800311}
#'
#' @examples
#' # TV (Beta-Binomial, 5 insertions) and a digital vehicle (Negative Binomial,
#' # 50 unit blocks), reach in proportions.
#' fit <- distrib_hybrid_copula(
#'   r1 = c(TV = 0.20, Digital = 0.02),
#'   r2 = c(TV = 0.312, Digital = 0.0390),
#'   insertions = c(5, 50),
#'   model_type = c("bbd", "nbd"),
#'   dup_matrix = matrix(c(0.20, 0.006, 0.006, 0.02), 2),
#'   n_qmc = 1024
#' )
#' fit
#'
#' # Two Beta-Binomial vehicles, reach in percent.
#' distrib_gcbbd(
#'   r1 = c(A = 20, B = 15), r2 = c(A = 31.2, B = 23.925),
#'   insertions = c(5, 5), dup_matrix = matrix(c(20, 6, 6, 15), 2),
#'   scale = "percent", n_qmc = 1024
#' )
#' @export
distrib_hybrid_copula <- function(r1, r2, insertions, model_type = "bbd",
                                  dup_matrix = NULL,
                                  dup_type = c("duplication", "reach"),
                                  scale = c("proportion", "percent"),
                                  k_max = NULL, n_qmc = 4096L) {
  dup_type <- match.arg(dup_type)
  scale <- match.arg(scale)
  fac <- if (scale == "percent") 100 else 1

  inp <- .check_inputs(r1, r2, insertions, model_type, fac, scale, k_max, n_qmc)
  m <- inp$m
  nm <- inp$names
  margs <- inp$marginals
  p <- inp$p

  d_obs <- .prepare_dup(dup_matrix, dup_type, fac, p, nm)

  grid <- .normal_grid()
  sigma <- diag(m)
  dimnames(sigma) <- list(nm, nm)
  pair_tab <- NULL
  psd_adjusted <- FALSE
  quad_err <- NULL

  if (m > 1L) {
    curves <- lapply(margs, .propensity_curve)
    quad_err <- vapply(seq_len(m), function(i) {
      abs(sum(grid$w * curves[[i]](grid$x)) / p[i] - 1)
    }, 0)
    names(quad_err) <- nm
    for (i in which(quad_err > 1e-3)) {
      warning(sprintf(paste0("Vehicle '%s' is strongly polarized (r2 very close ",
                             "to r1): the quadrature of its latent propensity ",
                             "has a relative error of %.2g in r1, so the ",
                             "calibrated dependence may be imprecise."),
                      nm[i], quad_err[i]), call. = FALSE)
    }
    idx <- which(upper.tri(sigma), arr.ind = TRUE)
    pair_tab <- data.frame(
      i = nm[idx[, 1]], j = nm[idx[, 2]],
      target = d_obs[idx], attained = NA_real_, status = NA_character_,
      stringsAsFactors = FALSE
    )
    for (r in seq_len(nrow(idx))) {
      i <- idx[r, 1]
      j <- idx[r, 2]
      cal <- .calibrate_rho(curves[[i]], curves[[j]], d_obs[i, j], grid)
      sigma[i, j] <- sigma[j, i] <- cal$rho
      pair_tab$attained[r] <- cal$attained
      pair_tab$status[r] <- cal$status
      if (cal$status != "ok") {
        warning(sprintf(
          paste0("Duplication between '%s' and '%s' is not attainable by the ",
                 "Gaussian copula with these marginals (target %.6g, nearest ",
                 "attainable %.6g); the nearest bound is used."),
          nm[i], nm[j], d_obs[i, j], cal$attained), call. = FALSE)
      }
    }
    ev <- eigen(sigma, symmetric = TRUE)
    if (min(ev$values) < 1e-8) {
      psd_adjusted <- TRUE
      warning(paste0("The pairwise calibrated latent correlation matrix is not ",
                     "positive definite; eigenvalues were floored at 1e-8, so ",
                     "the fitted duplications differ from the targets."),
              call. = FALSE)
      sigma <- cov2cor(ev$vectors %*% diag(pmax(ev$values, 1e-8)) %*% t(ev$vectors))
      dimnames(sigma) <- list(nm, nm)
    }
  }

  open_last_bin <- any(inp$types == "nbd")
  n_total <- sum(vapply(margs, function(g) if (g$type == "bbd") g$n else 0, 0))
  if (is.null(k_max)) {
    k_max <- if (open_last_bin) .default_k_max(margs) else n_total
  }
  if (!open_last_bin && k_max >= n_total) {
    k_max <- n_total
  } else {
    open_last_bin <- TRUE
  }
  kc <- if (open_last_bin) k_max else k_max + 1L    # exactly computed columns

  if (m == 1L) {
    pmf <- .exact_pmf_single(margs[[1]], kc)
    d_fit <- d_obs
    reach1_err <- 0
    dup_err <- 0
  } else {
    sim <- .simulate_pmf(margs, sigma, kc, as.integer(n_qmc))
    pmf <- sim$pmf
    d_fit <- sim$d_fit
    dimnames(d_fit) <- list(nm, nm)
    reach1_err <- max(abs(sim$mean_prop - p))
    dup_err <- max(abs(d_fit - d_obs)[upper.tri(d_obs)])
  }
  if (open_last_bin) pmf <- c(pmf, max(0, 1 - sum(pmf)))
  pmf <- pmax(pmf, 0)
  pmf <- pmf / sum(pmf)

  .build_result(pmf, margs, nm, sigma, d_obs, d_fit, pair_tab, psd_adjusted,
                reach1_err, dup_err, open_last_bin, fac, scale, quad_err)
}

#' @rdname distrib_hybrid_copula
#' @export
distrib_gcbbd <- function(r1, r2, insertions, dup_matrix = NULL,
                          dup_type = c("duplication", "reach"),
                          scale = c("proportion", "percent"),
                          k_max = NULL, n_qmc = 4096L) {
  distrib_hybrid_copula(r1 = r1, r2 = r2, insertions = insertions,
                        model_type = "bbd", dup_matrix = dup_matrix,
                        dup_type = match.arg(dup_type), scale = match.arg(scale),
                        k_max = k_max, n_qmc = n_qmc)
}

# --- Input handling ----------------------------------------------------------

.check_inputs <- function(r1, r2, insertions, model_type, fac, scale, k_max, n_qmc) {
  if (!is.numeric(r1) || !is.numeric(r2) || !is.numeric(insertions)) {
    stop("'r1', 'r2' and 'insertions' must be numeric.", call. = FALSE)
  }
  m <- length(r1)
  if (m < 1L || length(r2) != m || length(insertions) != m) {
    stop("'r1', 'r2' and 'insertions' must be non-empty and of equal length.",
         call. = FALSE)
  }
  if (anyNA(c(r1, r2, insertions)) || !all(is.finite(c(r1, r2, insertions)))) {
    stop("'r1', 'r2' and 'insertions' must be finite.", call. = FALSE)
  }
  if (scale == "proportion" && any(c(r1, r2) > 1)) {
    stop("Values above 1 found with scale = \"proportion\"; use ",
         "scale = \"percent\" if the inputs are percentages.", call. = FALSE)
  }
  p <- as.numeric(r1) / fac
  q <- as.numeric(r2) / fac
  if (any(p <= 0 | p >= 1) || any(q <= 0 | q >= 1)) {
    stop("'r1' and 'r2' must lie strictly between 0 and 1 (0 and 100 with ",
         "scale = \"percent\").", call. = FALSE)
  }
  if (any(insertions <= 0)) stop("'insertions' must be positive.", call. = FALSE)

  types <- tolower(as.character(model_type))
  if (length(types) == 1L) types <- rep(types, m)
  if (length(types) != m || !all(types %in% c("bbd", "nbd"))) {
    stop("'model_type' must be \"bbd\" or \"nbd\" for every vehicle.",
         call. = FALSE)
  }
  if (any(types == "bbd" & abs(insertions - round(insertions)) > 1e-8)) {
    stop("'insertions' must be integer valued for \"bbd\" vehicles.",
         call. = FALSE)
  }
  if (!is.null(k_max) && (length(k_max) != 1L || !is.finite(k_max) || k_max < 1)) {
    stop("'k_max' must be a single integer of at least 1.", call. = FALSE)
  }
  if (length(n_qmc) != 1L || !is.finite(n_qmc) || n_qmc < 64) {
    stop("'n_qmc' must be a single number of at least 64.", call. = FALSE)
  }
  nm <- names(r1)
  if (is.null(nm)) nm <- paste0("V", seq_len(m))

  margs <- lapply(seq_len(m), function(i) {
    if (types[i] == "bbd") {
      .marginal_bbd(p[i], q[i], round(insertions[i]), nm[i])
    } else {
      .marginal_nbd(p[i], q[i], insertions[i], nm[i])
    }
  })
  list(m = m, p = p, names = nm, types = types, marginals = margs)
}

.prepare_dup <- function(dup_matrix, dup_type, fac, p, nm) {
  m <- length(p)
  if (is.null(dup_matrix) || m == 1L) {
    d <- outer(p, p)
  } else {
    d <- as.matrix(dup_matrix)
    if (!is.numeric(d) || any(dim(d) != m) || anyNA(d)) {
      stop("'dup_matrix' must be a numeric m x m matrix without missing values.",
           call. = FALSE)
    }
    if (max(abs(d - t(d))) > 1e-8 * max(1, max(abs(d)))) {
      stop("'dup_matrix' must be symmetric.", call. = FALSE)
    }
    d <- d / fac
    if (dup_type == "reach") d <- outer(p, p, "+") - d
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        lo <- max(0, p[i] + p[j] - 1)
        hi <- min(p[i], p[j])
        if (d[i, j] < lo - 1e-10 || d[i, j] > hi + 1e-10) {
          stop(sprintf(paste0("Duplication between '%s' and '%s' (%.6g) is ",
                              "impossible for their reaches: it must lie in ",
                              "[%.6g, %.6g]."), nm[i], nm[j], d[i, j], lo, hi),
               call. = FALSE)
        }
      }
    }
  }
  diag(d) <- p
  dimnames(d) <- list(nm, nm)
  d
}

.default_k_max <- function(margs) {
  k <- sum(vapply(margs, function(g) {
    if (g$type == "bbd") g$n else qnbinom(1 - 1e-6, size = g$alpha,
                                          mu = g$mean_contacts)
  }, 0))
  k <- as.integer(ceiling(k)) + 1L
  if (k > 3000L) {
    stop("The default 'k_max' exceeds 3000; set 'k_max' explicitly.",
         call. = FALSE)
  }
  k
}

# --- Distribution ------------------------------------------------------------

.exact_pmf_single <- function(mg, kc) {
  k <- 0:(kc - 1L)
  if (mg$type == "bbd") {
    out <- numeric(kc)
    ok <- k <= mg$n
    kk <- k[ok]
    out[ok] <- exp(lchoose(mg$n, kk) + lbeta(kk + mg$alpha, mg$n - kk + mg$beta) -
                     lbeta(mg$alpha, mg$beta))
    out
  } else {
    dnbinom(k, size = mg$alpha, mu = mg$mean_contacts)
  }
}

.simulate_pmf <- function(margs, sigma, kc, n_qmc) {
  m <- length(margs)
  half <- ceiling(n_qmc / 2)
  z0 <- qnorm(.halton(half, m))
  z <- rbind(z0, -z0)                           # antithetic: mean exactly zero
  n_sim <- nrow(z)
  s_hat <- crossprod(z) / (n_sim - 1L)
  z <- z %*% backsolve(chol(s_hat), diag(m))    # unit sample covariance
  z <- z %*% chol(sigma)                        # target latent correlation

  prop <- matrix(0, n_sim, m)
  acc <- NULL
  for (i in seq_len(m)) {
    mg <- margs[[i]]
    lat <- .latent(mg, z[, i])
    prop[, i] <- .propensity(mg, lat)
    ks <- rep(0:(kc - 1L), each = n_sim)
    pm <- if (mg$type == "bbd") {
      matrix(dbinom(ks, size = mg$n, prob = rep(lat, kc)), n_sim, kc)
    } else {
      matrix(dpois(ks, lambda = rep(mg$n * lat, kc)), n_sim, kc)
    }
    acc <- if (is.null(acc)) pm else .row_convolve(acc, pm, kc)
  }
  list(pmf = colMeans(acc), d_fit = crossprod(prop) / n_sim,
       mean_prop = colMeans(prop))
}

.build_result <- function(pmf, margs, nm, sigma, d_obs, d_fit, pair_tab,
                          psd_adjusted, reach1_err, dup_err, open_last_bin,
                          fac, scale, quad_err) {
  k <- seq_along(pmf) - 1L
  cdf <- pmin(cumsum(pmf), 1)
  reach_k <- c(1, 1 - cdf[-length(cdf)])
  dist <- data.frame(k = k, pmf = pmf * fac, cdf = cdf * fac,
                     reach_k_plus = reach_k * fac)

  mean_contacts <- sum(vapply(margs, function(g) g$mean_contacts, 0))
  reach1 <- 1 - pmf[1]
  reach3 <- if (length(pmf) >= 4L) 1 - sum(pmf[1:3]) else NA_real_
  variance <- if (open_last_bin) NA_real_ else sum((k - sum(k * pmf))^2 * pmf)
  metrics <- list(
    grp = 100 * mean_contacts,
    reach_1plus = reach1 * fac,
    reach_3plus = reach3 * fac,
    ots = if (reach1 > 0) mean_contacts / reach1 else NA_real_,
    mean_contacts = mean_contacts,
    variance = variance,
    max_contacts = max(k),
    open_last_bin = open_last_bin
  )
  marginals <- data.frame(
    vehicle = nm,
    family = vapply(margs, function(g) g$type, ""),
    insertions = vapply(margs, function(g) g$n, 0),
    r1 = vapply(margs, function(g) g$p, 0) * fac,
    r2 = vapply(margs, function(g) g$q, 0) * fac,
    mean_contacts = vapply(margs, function(g) g$mean_contacts, 0),
    alpha = vapply(margs, function(g) g$alpha, 0),
    beta = vapply(margs, function(g) g$beta, 0),
    stringsAsFactors = FALSE
  )
  structure(
    list(distribution = dist, metrics = metrics, marginals = marginals,
         Sigma_latent = sigma, dup_observed = d_obs * fac,
         dup_fitted = d_fit * fac,
         diagnostics = list(pairs = pair_tab, psd_adjusted = psd_adjusted,
                            reach1_error = reach1_err, dup_error = dup_err,
                            quadrature_error = quad_err),
         scale = scale),
    class = "mp_copula_exposure"
  )
}

#' @export
print.mp_copula_exposure <- function(x, digits = 2, max_rows = 16, ...) {
  mt <- x$metrics
  unit <- if (x$scale == "percent") "%" else ""
  cat("Copula exposure distribution (mediaPlanCopula)\n")
  cat("Vehicles: ",
      paste0(x$marginals$vehicle, " (", toupper(x$marginals$family), ")",
             collapse = ", "), "\n", sep = "")
  cat(sprintf(paste0("GRPs: %.*f | Reach 1+: %.*f%s | Reach 3+: %s | ",
                     "OTS: %.*f\n\n"),
              digits, mt$grp, digits, mt$reach_1plus, unit,
              if (is.na(mt$reach_3plus)) "NA" else
                sprintf("%.*f%s", digits, mt$reach_3plus, unit),
              digits, mt$ots))
  shown <- x$distribution[seq_len(min(nrow(x$distribution), max_rows)), ]
  print(round(shown, digits + 2L), row.names = FALSE)
  if (nrow(x$distribution) > max_rows) {
    cat(sprintf("... (%d more rows)\n", nrow(x$distribution) - max_rows))
  }
  if (mt$open_last_bin) {
    cat(sprintf("The last category (k = %d) is '%d or more'.\n",
                mt$max_contacts, mt$max_contacts))
  }
  invisible(x)
}
