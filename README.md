# mediaPlanCopula

Exposure distribution of a media plan that mixes vehicles with a bounded number
of insertions (Beta-Binomial) and vehicles with an unbounded stream of impressions
(Negative Binomial), linked by a Gaussian copula.

## Installation

```r
# install.packages("remotes")
remotes::install_github("majesus/mediaPlanCopula")
```

## Example

```r
library(mediaPlanCopula)

fit <- distrib_hybrid_copula(
  r1 = c(TV = 0.20, Digital = 0.02),
  r2 = c(TV = 0.312, Digital = 0.0390),
  insertions = c(5, 50),
  model_type = c("bbd", "nbd"),
  dup_matrix = matrix(c(0.20, 0.006, 0.006, 0.02), 2)
)
fit
fit$diagnostics   # calibration errors and whether the correlation matrix was adjusted
```

## What it guarantees and what it does not

* The reach of one and of two insertions of each vehicle is reproduced exactly.
* The duplication of first insertions is calibrated pair by pair. When the
  calibrated correlation matrix is not positive definite it is adjusted, a warning
  is issued, and `diagnostics$dup_error` reports the largest deviation.
* The reach of one and two insertions and the duplication do not determine the
  distribution of contacts; its shape comes from the model's assumptions. Other
  multivehicle models with the same inputs give other distributions.
* In comparisons with other multivehicle models on survey-derived data, differences
  in accuracy were small and not stable out of sample. The function is offered as
  an implementation available for every plan, with explicit diagnostics, and not
  as a more accurate alternative.

The companion package [mediaPlanSim](https://github.com/majesus/mediaPlanSim)
generates synthetic exposure data with a known truth to study how models behave.

See `vignette("calibration-and-limits", package = "mediaPlanCopula")`.

## License

MIT.
