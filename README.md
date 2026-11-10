# mediaPlanCopula

Exposure (contact) distributions for media plans that mix vehicles with a
bounded number of insertions (Beta-Binomial) and vehicles with an unbounded
stream of impressions (Negative Binomial), linked by a Gaussian copula whose
latent correlations are calibrated to the observed duplication.

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
```

See `?distrib_hybrid_copula` for the model, its assumptions and diagnostics.
