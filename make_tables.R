#!/usr/bin/env Rscript
## Build the manuscript's tables and figure data from results/summary.json and
## results/illustrative.json, and substitute them into the LaTeX template
## (markers of the form %%NAME%%).
## Usage: Rscript make_tables.R results manuscript_template.tex manuscript.tex
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly = TRUE)
res <- args[1]; tmpl <- args[2]; outp <- args[3]
S <- fromJSON(file.path(res, "summary.json"), simplifyVector = FALSE)
I <- fromJSON(file.path(res, "illustrative.json"), simplifyVector = FALSE)
LAMS <- c("0.0", "0.3", "0.6", "0.9")
EST <- c("FE", "IV", "DGMM", "SGMM")
ENAME <- c(FE = "FE", IV = "IV", DGMM = "Diff-GMM", SGMM = "Sys-GMM")
CN <- c(alpha = "$\\alpha$", b1 = "$\\beta_1$", b2 = "$\\beta_2$", b3 = "$\\beta_3$")
pf <- S$perf; V <- S$valid
un <- function(x) unlist(x)

f <- function(x, d = 4) { x <- round(x, d); if (x == 0) x <- 0; s <- formatC(x, format = "f", digits = d); if (x < 0) sub("-", "$-$", s, fixed = TRUE) else s }
big <- function(x, d = 4) if (abs(x) < 10) f(x, d) else f(x, 2)
cm <- function(x, d = 0) gsub(",", "{,}", formatC(x, format = "f", digits = d, big.mark = ","), fixed = TRUE)
bold <- function(s, b) if (b) paste0("\\textbf{", s, "}") else s
argmin <- function(lam, c) EST[which.min(sapply(EST, function(e) pf[[lam]][[e]][[c]]$rmse))]

T <- list()
## --- performance table for alpha and beta1
rows <- character()
for (lam in LAMS) {
  rows <- c(rows, "\\addlinespace", sprintf("\\multicolumn{9}{@{}l}{\\textit{$\\lambda=%s$}} \\\\", lam))
  for (e in EST) for (c in c("alpha", "b1")) {
    p <- pf[[lam]][[e]][[c]]
    rows <- c(rows, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %.3f \\\\",
      if (c == "alpha") ENAME[[e]] else "", CN[[c]], big(p$bias), f(p$rb, 1), big(p$sd),
      bold(big(p$rmse), e == argmin(lam, c)), f(p$medbias), f(p$mae), p$cover))
  }
}
T$TAB_PERF <- paste(rows, collapse = "\n")

## --- climatic coefficients
rows <- character()
for (c in c("b2", "b3")) {
  for (e in EST) {
    cells <- unlist(lapply(LAMS, function(lam) {
      p <- pf[[lam]][[e]][[c]]
      c(f(p$bias), bold(f(p$rmse), e == argmin(lam, c)), sprintf("%.2f", p$cover))
    }))
    rows <- c(rows, paste0(if (e == "FE") CN[[c]] else "", " & ", ENAME[[e]], " & ",
                           paste(cells, collapse = " & "), " \\\\"))
  }
  if (c == "b2") rows <- c(rows, "\\addlinespace")
}
T$TAB_CLIM <- paste(rows, collapse = "\n")

## --- validation table
row <- function(label, fn) paste0(label, " & ", paste(sapply(LAMS, fn), collapse = " & "), " \\\\")
rows <- c(
  row("Retained countries, mean [min, max]", function(l) sprintf("%.1f [%d, %d]", V[[l]]$Nret_mean, V[[l]]$Nret_min, V[[l]]$Nret_max)),
  row("Observed country-years, mean", function(l) cm(V[[l]]$nobs_mean)),
  row("Usable rows: FE / IV, Diff-GMM / Sys-GMM", function(l) paste(cm(V[[l]]$FE$nobs_mean), cm(V[[l]]$IV$nobs_mean), cm(V[[l]]$SGMM$nobs_mean), sep = " / ")),
  row("Successful fits (each estimator)", function(l) {
    s <- sapply(EST, function(e) V[[l]][[e]]$success); if (all(s == 1000)) "1{,}000" else paste(s, collapse = "/") }),
  "\\addlinespace",
  row("IV: median first-stage $F$ ($\\Delta y_{t-1}$ / $\\Delta q_t$)", function(l) sprintf("%.1f / %.1f", un(V[[l]]$IV$F_med)[1], un(V[[l]]$IV$F_med)[2])),
  row("IV: share with $F<10$ ($\\Delta y_{t-1}$ / $\\Delta q_t$)", function(l) sprintf("%.2f / %.2f", un(V[[l]]$IV$F_lt10)[1], un(V[[l]]$IV$F_lt10)[2])))
for (e in c("DGMM", "SGMM")) rows <- c(rows, "\\addlinespace",
  row(paste0(ENAME[[e]], ": instruments (Hansen df)"), function(l) sprintf("%d (%d)", un(V[[l]][[e]]$ninst)[1], un(V[[l]][[e]]$Jdf)[1])),
  row(paste0(ENAME[[e]], ": Hansen rejection, two-step / iterated"), function(l) sprintf("%.3f / %.3f", V[[l]][[e]]$J_rej, V[[l]][[e]]$Jit_rej)),
  row(paste0(ENAME[[e]], ": AR(1) / AR(2) rejection"), function(l) sprintf("%.3f / %.3f", V[[l]][[e]]$ar1_rej, V[[l]][[e]]$ar2_rej)))
rows <- c(rows, "\\addlinespace",
  row("FE residuals: Pesaran CD rejection", function(l) sprintf("%.3f", V[[l]]$FE$cd_rej)),
  row("FE residuals: median CD statistic", function(l) f(V[[l]]$FE$cd_med, 2)),
  row("Wooldridge AR test rejection (median $\\hat\\rho$)", function(l) sprintf("%.3f (%s)", V[[l]]$FE$wool_rej, f(V[[l]]$FE$wool_rho_med, 2))))
T$TAB_VALID <- paste(rows, collapse = "\n")

## --- oracle moment table
O <- S$oracle
lab <- c(diff_y2 = "$y_{i,t-2}$ vs $\\Delta\\varepsilon_{it}$", diff_y3 = "$y_{i,t-3}$ vs $\\Delta\\varepsilon_{it}$",
         diff_q2 = "$q_{i,t-2}$ vs $\\Delta\\varepsilon_{it}$", diff_q3 = "$q_{i,t-3}$ vs $\\Delta\\varepsilon_{it}$",
         lev_dy1 = "$\\Delta y_{i,t-1}$ vs $\\mu_i+\\varepsilon_{it}$", lev_dq1 = "$\\Delta q_{i,t-1}$ vs $\\mu_i+\\varepsilon_{it}$",
         corr_q_eps = "$q_{it}$ vs $\\varepsilon_{it}$ (endogeneity, reference)")
rows <- character()
for (k in names(lab)) {
  rows <- c(rows, paste0(lab[[k]], " & ", paste(sapply(LAMS, function(l)
    sprintf("%s & %.3f", f(O[[l]][[k]]$mean), O[[l]][[k]]$p95abs)), collapse = " & "), " \\\\"))
  if (k == "lev_dq1") rows <- c(rows, "\\addlinespace")
}
T$TAB_ORACLE <- paste(rows, collapse = "\n")

## --- illustrative replication
de <- I$desc
drow <- function(lab, s) sprintf("%s & %s & %.2f & %.1f & %s \\\\", lab, cm(de$GDP[[s]]),
                                 de$CO2[[s]], de$Precip[[s]], f(de$Temp[[s]]))
T$TAB_DESC <- paste(c(drow("Minimum", "min"), drow("Median", "median"), drow("Mean", "mean"),
                      drow("Std.\\ deviation", "sd"), drow("Maximum", "max")), collapse = "\n")
a <- I$ar1
T$TAB_AR1 <- paste(mapply(function(k, nm) sprintf("%s & %.3f & %.4f \\\\", nm, a[[k]]$b, a[[k]]$se),
  c("GDP", "CO2", "Precip", "Temp"),
  c("Log GDP", "Log \\COtwo{}", "Log precipitation", "Temperature anomaly")), collapse = "\n")
E <- I$est
names_row <- c("Lag log GDP ($\\alpha=0.70$)", "Log \\COtwo{} ($\\beta_1=0.70$)",
               "Log precipitation ($\\beta_2=0.003$)", "Temp.\\ anomaly ($\\beta_3=-0.015$)")
rows <- character()
for (j in 1:4) {
  cells <- ses <- character()
  for (e in EST) {
    b <- un(E[[e]]$b)[j]; s <- un(E[[e]]$se)[j]
    pv <- 2 * pnorm(-abs(b / s))
    st <- if (pv < .01) "***" else if (pv < .05) "**" else if (pv < .1) "*" else ""
    cells <- c(cells, paste0(f(b), if (nzchar(st)) sprintf("$^{%s}$", st) else ""))
    ses <- c(ses, sprintf("(%.4f)", s))
  }
  rows <- c(rows, paste0(names_row[j], " & ", paste(cells, collapse = " & "), " \\\\"),
            paste0(" & ", paste(ses, collapse = " & "), " \\\\"))
}
T$TAB_SINGLE <- paste(rows, collapse = "\n")
pz <- function(z) { pv <- 2 * pnorm(-abs(z)); if (pv < 0.001) "$<0.001$" else sprintf("%.3f", pv) }
dg <- E$DGMM; sg <- E$SGMM; iv <- E$IV
cd <- un(I$cd); wo <- un(I$wooldridge); G1 <- I$Nret - 1
fF <- un(iv$F1)
T$TAB_DIAG <- paste(c(
  sprintf("FE & Pesaran CD & $z=%.2f$ & -- & %.3f \\\\", cd[1], cd[2]),
  sprintf("FE & Wooldridge AR & $F=%.1f$ & (1, %d) & %s \\\\", wo[1], G1, if (wo[2] < .001) "$<0.001$" else sprintf("%.3f", wo[2])),
  "\\midrule",
  sprintf("IV & First stage, $\\Delta y_{t-1}$ & $F=%.2f$ & (2, %d) & %.3f \\\\", fF[1], G1, pf(fF[1], 2, G1, lower.tail = FALSE)),
  sprintf("IV & First stage, $\\Delta q_t$ & $F=%.2f$ & (2, %d) & %.3f \\\\", fF[2], G1, pf(fF[2], 2, G1, lower.tail = FALSE)),
  "\\midrule",
  sprintf("Diff-GMM & Hansen (two-step) & $\\chi^2=%.2f$ & %d & %.3f \\\\", dg$J, dg$Jdf, dg$Jp),
  sprintf("Diff-GMM & AR(1) & $z=%.2f$ & -- & %s \\\\", dg$ar1, pz(dg$ar1)),
  sprintf("Diff-GMM & AR(2) & $z=%.2f$ & -- & %s \\\\", dg$ar2, pz(dg$ar2)),
  "\\midrule",
  sprintf("Sys-GMM & Hansen (two-step) & $\\chi^2=%.2f$ & %d & %.3f \\\\", sg$J, sg$Jdf, sg$Jp),
  sprintf("Sys-GMM & AR(1) & $z=%.2f$ & -- & %s \\\\", sg$ar1, pz(sg$ar1)),
  sprintf("Sys-GMM & AR(2) & $z=%.2f$ & -- & %s \\\\", sg$ar2, pz(sg$ar2))), collapse = "\n")

## --- figure coordinates
for (cc in list(c("alpha", "FIG_A"), c("b1", "FIG_B"))) {
  sfx <- c(FE = "_FE", IV = "_IV", DGMM = "_DG", SGMM = "_SG")
  for (e in EST) T[[paste0(cc[2], sfx[[e]])]] <- paste0("{", paste(sprintf("(%s,%.4f)", LAMS,
    sapply(LAMS, function(l) pf[[l]][[e]][[cc[1]]]$rmse)), collapse = ""), "}")
}

s <- paste(readLines(tmpl, warn = FALSE), collapse = "\n")
for (k in names(T)) s <- gsub(paste0("%%", k, "%%"), T[[k]], s, fixed = TRUE)
left <- regmatches(s, gregexpr("%%(TAB|FIG)_[A-Z_]+%%", s))[[1]]
if (length(left)) stop("unfilled markers: ", paste(left, collapse = ", "))
writeLines(s, outp)
cat("Wrote", outp, "\n")
