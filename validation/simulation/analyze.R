# Statistical analysis of the simulation: paired contrasts copula vs. each competitor.
files <- list.files("/home/claude/sim", pattern = "^res_.*\\.rds$", full.names = TRUE)
d <- do.call(rbind, lapply(files, readRDS))
d$absE1 <- abs(d$e1); d$absE3 <- abs(d$e3)
DELTA <- c(tv = 0.005, absE1 = 0.5, absE3 = 0.5)      # margins of practical equivalence (TV; percentage points)
models <- setdiff(unique(d$model), "copula")

# ---- descriptive ----
desc <- do.call(rbind, lapply(split(d, list(d$part, d$dgp, d$regime, d$model), drop = TRUE), function(x) {
  ok <- x[x$ok, ]
  data.frame(part = x$part[1], dgp = x$dgp[1], regime = x$regime[1], model = x$model[1], n = nrow(x),
             fail_rate = mean(!x$ok), warn_rate = mean(x$nwarn > 0),
             tv = mean(ok$tv), tv_se = sd(ok$tv) / sqrt(nrow(ok)),
             bias_r1 = mean(ok$e1), mae_r1 = mean(ok$absE1), bias_r3 = mean(ok$e3), mae_r3 = mean(ok$absE3))
}))
rownames(desc) <- NULL

# ---- paired contrasts ----
verdict <- function(lo, hi, delta) {
  if (lo > delta) "copula mejor (relevante)" else if (lo > 0) "copula mejor (pequeno)"
  else if (hi < -delta) "copula peor (relevante)" else if (hi < 0) "copula peor (pequeno)"
  else if (lo > -delta && hi < delta) "equivalentes" else "no concluyente"
}
contr <- list()
for (key in unique(paste(d$part, d$dgp, d$regime))) {
  s <- d[paste(d$part, d$dgp, d$regime) == key, ]
  cop <- s[s$model == "copula", ]
  for (mo in models) {
    oth <- s[s$model == mo, ]
    mm <- merge(cop[, c("id", "ok", "tv", "absE1", "absE3")], oth[, c("id", "ok", "tv", "absE1", "absE3")], by = "id", suffixes = c("_c", "_o"))
    mm <- mm[mm$ok_c & mm$ok_o, ]
    for (met in c("tv", "absE1", "absE3")) {
      dif <- mm[[paste0(met, "_o")]] - mm[[paste0(met, "_c")]]   # > 0: copula smaller error
      n <- length(dif); if (n < 10) next
      tt <- t.test(dif); wp <- suppressWarnings(wilcox.test(dif, exact = FALSE))$p.value
      contr[[length(contr) + 1L]] <- data.frame(part = s$part[1], dgp = s$dgp[1], regime = s$regime[1], metric = met, model = mo,
        n_pairs = n, mean_diff = mean(dif), ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2], p_t = tt$p.value, p_wilcoxon = wp,
        win_rate = mean(dif > 0), dz = mean(dif) / sd(dif), mde80 = (qnorm(.975) + qnorm(.8)) * sd(dif) / sqrt(n),
        verdict = verdict(tt$conf.int[1], tt$conf.int[2], DELTA[[met]]), stringsAsFactors = FALSE)
    }
  }
}
contr <- do.call(rbind, contr)
contr$p_holm <- ave(contr$p_wilcoxon, contr$metric, FUN = function(p) p.adjust(p, "holm"))
rownames(contr) <- NULL
write.csv(desc, "./tab_descriptivo.csv", row.names = FALSE)
write.csv(contr, "./tab_contrastes.csv", row.names = FALSE)
options(width = 200)
cat("\n===== TV medio (exposicion: error de la distribucion completa), regimen poblacion =====\n")
print(reshape(desc[desc$regime == "population", c("dgp", "model", "tv")], idvar = "dgp", timevar = "model", direction = "wide"), digits = 3, row.names = FALSE)
cat("\n===== TV medio, regimen panel 3000 =====\n")
print(reshape(desc[desc$regime == "panel3000", c("dgp", "model", "tv")], idvar = "dgp", timevar = "model", direction = "wide"), digits = 3, row.names = FALSE)
cat("\n===== Tasa de fallo (panel 3000) =====\n")
print(reshape(desc[desc$regime == "panel3000", c("dgp", "model", "fail_rate")], idvar = "dgp", timevar = "model", direction = "wide"), digits = 2, row.names = FALSE)
cat("\n===== Veredictos TV: copula frente a cada modelo (regimen poblacion) =====\n")
v <- contr[contr$metric == "tv" & contr$regime == "population", ]
print(reshape(v[, c("dgp", "model", "verdict")], idvar = "dgp", timevar = "model", direction = "wide"), row.names = FALSE)
cat("\n===== Veredictos TV (panel 3000) =====\n")
v <- contr[contr$metric == "tv" & contr$regime == "panel3000", ]
print(reshape(v[, c("dgp", "model", "verdict")], idvar = "dgp", timevar = "model", direction = "wide"), row.names = FALSE)
cat("\nMDE80 mediano (TV):", median(contr$mde80[contr$metric == "tv"]), " (max:", max(contr$mde80[contr$metric == "tv"]), ")\n")
