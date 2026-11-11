#!/usr/bin/env Rscript
## Summarise results/replications.jsonl into performance and validation measures
## (writes results/summary.json and prints a readable summary).
## Usage: Rscript summarize.R results
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results"

recs <- lapply(readLines(file.path(dir, "replications.jsonl")), fromJSON, simplifyVector = TRUE)
TRUEV <- c(alpha = 0.70, b1 = 0.70, b2 = 0.003, b3 = -0.015)
COEFS <- names(TRUEV); EST <- c("FE", "IV", "DGMM", "SGMM")
LAMS <- sort(unique(sapply(recs, `[[`, "lam")))
z <- qnorm(0.975)
key <- function(l) sprintf("%.1f", l)
ok_fit <- function(x) is.list(x)

out <- list(perf = list(), valid = list(), oracle = list())
for (lam in LAMS) {
  R <- Filter(function(r) r$lam == lam, recs)
  Nret <- sapply(R, `[[`, "Nret"); nobs <- sapply(R, `[[`, "nobs")
  v <- list(reps = length(R), Nret_mean = mean(Nret), Nret_min = min(Nret), Nret_max = max(Nret),
            nobs_mean = mean(nobs), nobs_min = min(nobs), nobs_max = max(nobs))
  perf_l <- list()
  for (e in EST) {
    ok <- Filter(ok_fit, lapply(R, `[[`, e))
    ve <- list(success = length(ok), failed = length(R) - length(ok),
               nobs_mean = mean(sapply(ok, `[[`, "nobs")))
    g <- function(f) sapply(ok, `[[`, f)
    if (e %in% c("DGMM", "SGMM")) {
      a1 <- g("ar1"); a2 <- g("ar2")
      ve <- c(ve, list(ninst = unique(g("ninst")), Jdf = unique(g("Jdf")),
                       J_rej = mean(g("Jp") < 0.05), Jit_rej = mean(g("Jit_p") < 0.05),
                       ar1_rej = mean(abs(a1[is.finite(a1)]) > z),
                       ar2_rej = mean(abs(a2[is.finite(a2)]) > z)))
    }
    if (e == "IV") {
      Fm <- t(sapply(ok, `[[`, "F1"))
      ve <- c(ve, list(ninst = 2, F_med = apply(Fm, 2, median), F_lt10 = colMeans(Fm < 10)))
    }
    if (e == "FE") {
      ve <- c(ve, list(cd_rej = mean(g("cd_p") < 0.05), cd_med = median(g("cd")),
                       absrho_med = median(g("abs_rho")), wool_rej = mean(g("wool_p") < 0.05),
                       wool_rho_med = median(g("wool_rho"))))
    }
    v[[e]] <- ve
    B <- t(sapply(ok, `[[`, "b")); S <- t(sapply(ok, `[[`, "se"))
    pe <- list()
    for (j in seq_along(COEFS)) {
      th <- TRUEV[[j]]; x <- B[, j]; s <- S[, j]
      pe[[COEFS[j]]] <- list(bias = mean(x) - th, rb = 100 * (mean(x) - th) / th, sd = sd(x),
                             rmse = sqrt(mean((x - th)^2)), mcse = sd(x) / sqrt(length(x)),
                             medbias = median(x) - th, mae = median(abs(x - th)),
                             cover = mean(abs(x - th) <= z * s))
    }
    perf_l[[e]] <- pe
  }
  out$valid[[key(lam)]] <- v
  out$perf[[key(lam)]] <- perf_l
  O <- lapply(R, `[[`, "oracle")
  out$oracle[[key(lam)]] <- lapply(setNames(names(O[[1]]), names(O[[1]])), function(k) {
    x <- sapply(O, `[[`, k); list(mean = mean(x), p95abs = unname(quantile(abs(x), 0.95)))
  })
}
write_json(out, file.path(dir, "summary.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)

for (lam in LAMS) {
  k <- key(lam); v <- out$valid[[k]]
  cat(sprintf("\n===== lambda %s | countries %.1f [%d,%d] | obs %.0f\n", k, v$Nret_mean, v$Nret_min, v$Nret_max, v$nobs_mean))
  cat(sprintf("  fits: %s\n", paste(sapply(EST, function(e) paste0(e, "=", v[[e]]$success)), collapse = " ")))
  cat(sprintf("  DGMM inst %d J %.3f/%.3f AR2 %.3f | SGMM inst %d J %.3f/%.3f AR2 %.3f | CD %.3f Wool %.3f (rho %.2f) | IV F %.1f/%.1f\n",
              v$DGMM$ninst, v$DGMM$J_rej, v$DGMM$Jit_rej, v$DGMM$ar2_rej, v$SGMM$ninst, v$SGMM$J_rej,
              v$SGMM$Jit_rej, v$SGMM$ar2_rej, v$FE$cd_rej, v$FE$wool_rej, v$FE$wool_rho_med,
              v$IV$F_med[1], v$IV$F_med[2]))
  for (e in EST) for (c in COEFS) {
    p <- out$perf[[k]][[e]][[c]]
    cat(sprintf("  %-5s%-6s bias %+.4f rb %+7.1f sd %.4f rmse %.4f mcse %.4f medb %+.4f mae %.4f cov %.3f\n",
                e, c, p$bias, p$rb, p$sd, p$rmse, p$mcse, p$medbias, p$mae, p$cover))
  }
}
