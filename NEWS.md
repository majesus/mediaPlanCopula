# mediaPlanCopula 0.1.0

* First release. `distrib_hybrid_copula()` and `distrib_gcbbd()` compute the
  exposure distribution of media plans that mix Beta-Binomial and Negative
  Binomial vehicles linked by a Gaussian copula.
* The reach of one and of two insertions of each vehicle is reproduced exactly;
  the duplication of first insertions is calibrated pair by pair and the
  calibration error is reported in `diagnostics`.
* Documentation states the limits of the model: approximate duplication when the
  calibrated correlation matrix is not positive definite, dependence of the
  distribution on the model's assumptions, and no claim of higher accuracy than
  other multivehicle models.
* A vignette illustrates calibration, diagnostics and numerical accuracy.
