suppressMessages({library(mediaPlanR); library(mediaPlanCopula); library(parallel)})
source("./harness.R")   # safe(), run_models(), score()
dir <- "/mnt/user-data/uploads/_copula_mediaPlanR/validation_real"
veh <- read.csv(file.path(dir, "vehicles.csv")); prs <- read.csv(file.path(dir, "pairs.csv")); sch <- read.csv(file.path(dir, "schedules.csv"))
Dfull <- matrix(0, 25, 25); for (k in seq_len(nrow(prs))) { Dfull[prs$i[k], prs$j[k]] <- prs$dup_pct[k] / 100; Dfull[prs$j[k], prs$i[k]] <- prs$dup_pct[k] / 100 }
r1all <- veh$avg_pct / 100; r2all <- veh$cum_pct / 100
run_one <- function(k) {
  v <- as.integer(strsplit(sch$veh[k], " ")[[1]]); o <- as.numeric(strsplit(sch$obs[k], " ")[[1]]); o <- o / sum(o)
  m <- length(v); inp <- list(r1 = r1all[v], r2 = r2all[v], D = Dfull[v, v, drop = FALSE]); diag(inp$D) <- inp$r1
  ni <- rep(2L, m); types <- rep("bbd", m)
  t0 <- proc.time()[[3]]
  ms <- run_models(inp, ni, types, NULL)
  do.call(rbind, lapply(names(ms), function(nm) { s <- score(ms[[nm]]$pmf, o)
    data.frame(sid = k, m = m, model = nm, ok = !is.null(ms[[nm]]$pmf) && !is.na(s[["tv"]]), nwarn = ms[[nm]]$nw, tv = s[["tv"]], e1 = s[["e1"]], e3 = s[["e3"]],
               obs_reach = 100 * (1 - o[1]), sec = proc.time()[[3]] - t0, stringsAsFactors = FALSE) }))
}
args <- commandArgs(TRUE); ids <- if (length(args)) eval(parse(text = args[1])) else seq_len(nrow(sch))
res <- do.call(rbind, mclapply(ids, run_one, mc.cores = 2, mc.preschedule = FALSE))
saveRDS(res, if (length(args) > 1) args[2] else "./res_real_cheong.rds")
print(aggregate(cbind(tv, sec) ~ model, res, mean))
