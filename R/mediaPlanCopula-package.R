#' mediaPlanCopula: Copula-Based Exposure Distributions for Mixed Media Plans
#'
#' Exposure distributions for plans that combine vehicles with a bounded number
#' of insertions (Beta-Binomial, BBD) and vehicles that deliver an unbounded
#' stream of impressions (Negative Binomial, NBD), linked by a Gaussian copula.
#' See [distrib_hybrid_copula()] for the model and its assumptions.
#'
#' @keywords internal
#' @importFrom stats cov2cor dbinom dnbinom dnorm dpois pnorm qbeta qgamma
#'   qnbinom qnorm splinefun uniroot
"_PACKAGE"
