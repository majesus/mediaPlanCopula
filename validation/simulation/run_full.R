source("./harness.R")
library(parallel)
R <- 150L
specs <- list(
  list(part="A", dgp="G1_gauss_beta", seed=100000), list(part="A", dgp="G2_gauss_logitnormal", seed=200000),
  list(part="A", dgp="G3_t_copula", seed=300000), list(part="A", dgp="G4_latent_class", seed=400000),
  list(part="A", dgp="G5_decay", seed=500000), list(part="A", dgp="G6_negative", seed=600000),
  list(part="B", dgp="H1_gauss_beta_gamma", seed=700000), list(part="B", dgp="H2_latent_class", seed=800000))
for (sp in specs) {
  f <- sprintf("./res_%s.rds", sp$dgp)
  if (file.exists(f)) next
  t0 <- proc.time()[[3]]
  out <- mclapply(seq_len(R), function(i) tryCatch(one_scenario(i, sp$part, sp$dgp, c(Inf, 3000), sp$seed), error = function(e) NULL), mc.cores = 2, mc.preschedule = FALSE)
  saveRDS(do.call(rbind, out), f)
  cat(format(Sys.time(), "%H:%M:%S"), sp$dgp, "listo en", round(proc.time()[[3]] - t0), "s;", sum(vapply(out, is.null, TRUE)), "escenarios fallidos\n", file = "./progress.log", append = TRUE)
}
cat("FIN\n", file = "./progress.log", append = TRUE)
