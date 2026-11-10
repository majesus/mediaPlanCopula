res <- readRDS("./res_real_cheong.rds")
res$aer <- 100 * abs(res$e1) / res$obs_reach
res$grp <- cut(res$m, c(1, 4, 8, 12), labels = c("2-4 soportes", "5-8 soportes", "9-12 soportes"))
options(width = 190, digits = 3)
cat("\n=== Tasa de fallo por modelo (n = 440) ===\n"); print(round(tapply(!res$ok, res$model, mean), 3))
cat("\n=== TV medio (solo escenarios donde el modelo converge) ===\n")
print(round(sapply(split(res[res$ok, ], res$model[res$ok]), function(x) mean(x$tv)), 4))
cat("\n=== TV medio por tamano de programa ===\n")
print(round(with(res[res$ok, ], tapply(tv, list(grp, model), mean)), 4))
cat("\n=== AER (error relativo medio del alcance, %) ===\n")
print(round(sapply(split(res[res$ok, ], res$model[res$ok]), function(x) mean(x$aer)), 2))
cat("\n=== Sesgo del alcance 3+ (pp) ===\n")
print(round(sapply(split(res[res$ok, ], res$model[res$ok]), function(x) mean(x$e3)), 2))
# paired contrasts vs copula, restricted to schedules where both ok
cop <- res[res$model == "copula", ]
out <- do.call(rbind, lapply(setdiff(unique(res$model), "copula"), function(mo) {
  o <- res[res$model == mo, ]; mm <- merge(cop[, c("sid", "ok", "tv", "aer", "grp")], o[, c("sid", "ok", "tv", "aer")], by = "sid", suffixes = c("_c", "_o"))
  mm <- mm[mm$ok_c & mm$ok_o, ]; d <- mm$tv_o - mm$tv_c; n <- length(d); tt <- t.test(d)
  data.frame(model = mo, n_pairs = n, copula_TV = mean(mm$tv_c), rival_TV = mean(mm$tv_o), mean_diff = mean(d), ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
             win_rate = mean(d > 0), p_wilcoxon = suppressWarnings(wilcox.test(d, exact = FALSE)$p.value), dz = mean(d) / sd(d))
}))
out$p_holm <- p.adjust(out$p_wilcoxon, "holm"); cat("\n=== Contrastes pareados TV (positivo = la copula comete menos error) ===\n"); print(out, row.names = FALSE)
write.csv(out, "./tab_real_contrastes.csv", row.names = FALSE)
write.csv(res, "./tab_real_resultados.csv", row.names = FALSE)
cat("\nAdvertencias copula (schedules con >=1):", mean(cop$nwarn > 0), "\n")
