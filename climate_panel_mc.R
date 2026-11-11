#!/usr/bin/env Rscript
## =============================================================================
## Monte Carlo comparison of FE, Anderson-Hsiao IV, Difference-GMM and
## System-GMM in unbalanced dynamic climate-economy panels.
##
## Reproduces every simulation result in:
##   Opoku, Adu, Emmanuel & Owusu, "Evaluating Estimator Performance in
##   Unbalanced Dynamic Climate-Economy Panels", Journal of Science and
##   Technology (KNUST).
##
## Requirements: R >= 4.0 with the 'jsonlite' package ('parallel' is base R).
## Usage (from a terminal):
##   Rscript climate_panel_mc.R --reps 1000 --out results      # full study
##   Rscript climate_panel_mc.R --illustrative --out results   # one replication
## or, inside R:  source("climate_panel_mc.R"); run_study(1000, "results")
##
## Random numbers: R's default Mersenne-Twister generator (R >= 3.6 sampling).
## Replication r (0..999) of setting k (lambda index 0..3) is generated after
##   set.seed(20261030 + 100000 * k + r)
## so any replication can be regenerated on its own, independently of the
## number of cores used.
## All estimators and tests are coded from their published formulae in base R,
## so the results do not depend on package defaults.
## =============================================================================

suppressPackageStartupMessages(library(jsonlite))

## ---------------------------------------------------------------------------
## 1. Design constants (Table 1 of the manuscript)
## ---------------------------------------------------------------------------
MASTER_SEED <- 20261030
LAMBDAS <- c(0, 0.3, 0.6, 0.9)
N0 <- 135; TT <- 50; BURN <- 50
TRUE_THETA <- c(alpha = 0.70, b1 = 0.70, b2 = 0.003, b3 = -0.015)
P <- list(
  sigma_mu = 0.30,                                   # sd of country effect mu_i
  mq0 = 5.55, sigma_mq = 0.30, corr_mq_mu = 0.60,    # country mean of log CO2
  rho_q = 0.80, sigma_eta = 0.10,                    # CO2 AR(1), innovation sd
  mp0 = log(870), sigma_mp = 0.25,                   # country mean of log precipitation
  rho_p = 0.30, sigma_zeta = 0.25,                   # precipitation AR(1), innovation sd
  rho_h = 0.60, sigma_xi = 0.70,                     # temperature anomaly AR(1), idios. sd
  rho_g = 0.50, sigma_g = 0.30,                      # global temperature component
  sigma_u = 0.05,                                    # GDP shock sd
  max_entry = 30, p_delete = 0.08, min_obs = 20)

rep_seed <- function(k, r) MASTER_SEED + 100000 * k + r

## ---------------------------------------------------------------------------
## 2. Data-generating process
## ---------------------------------------------------------------------------
simulate_panel <- function(lam, return_truth = FALSE) {
  a <- TRUE_THETA[["alpha"]]; b1 <- TRUE_THETA[["b1"]]
  b2 <- TRUE_THETA[["b2"]];   b3 <- TRUE_THETA[["b3"]]
  TB <- BURN + TT
  z <- rnorm(N0); w <- rnorm(N0); v <- rnorm(N0)
  mu <- P$sigma_mu * z
  cc <- P$corr_mq_mu
  mq <- P$mq0 + P$sigma_mq * (cc * z + sqrt(1 - cc^2) * w)
  mp <- P$mp0 + P$sigma_mp * v
  eta  <- matrix(rnorm(N0 * TB, sd = P$sigma_eta),  N0, TB)
  zeta <- matrix(rnorm(N0 * TB, sd = P$sigma_zeta), N0, TB)
  xi   <- matrix(rnorm(N0 * TB, sd = P$sigma_xi),   N0, TB)
  u    <- matrix(rnorm(N0 * TB, sd = P$sigma_u),    N0, TB)
  gsh  <- rnorm(TB, sd = P$sigma_g)
  q <- p <- h <- y <- matrix(NA_real_, N0, TB)
  ## initial values: country-specific steady states (mean-stationary start)
  q_prev <- mq; p_prev <- mp; h_prev <- rep(0, N0); g_prev <- 0
  y_prev <- (mu + b1 * mq + b2 * mp) / (1 - a)
  eps <- u + lam * eta
  for (t in seq_len(TB)) {
    g_t <- P$rho_g * g_prev + gsh[t]
    q[, t] <- mq + P$rho_q * (q_prev - mq) + eta[, t]
    p[, t] <- mp + P$rho_p * (p_prev - mp) + zeta[, t]
    h[, t] <- P$rho_h * h_prev + g_t + xi[, t]
    y[, t] <- a * y_prev + b1 * q[, t] + b2 * p[, t] + b3 * h[, t] + mu + eps[, t]
    q_prev <- q[, t]; p_prev <- p[, t]; h_prev <- h[, t]; g_prev <- g_t; y_prev <- y[, t]
  }
  keep_t <- (BURN + 1):TB
  y <- y[, keep_t]; q <- q[, keep_t]; p <- p[, keep_t]; h <- h[, keep_t]
  eps <- eps[, keep_t]
  ## imbalance: staggered entry and random intermittent deletion
  e <- sample.int(P$max_entry + 1, N0, replace = TRUE) - 1
  obs <- outer(e, seq_len(TT), function(ei, t) t > ei) &
         matrix(runif(N0 * TT) > P$p_delete, N0, TT)
  keep <- rowSums(obs) >= P$min_obs
  mask <- function(A) { A[!obs] <- NA; A[keep, , drop = FALSE] }
  out <- list(y = mask(y), q = mask(q), p = mask(p), h = mask(h))
  if (return_truth) { out$eps <- eps[keep, , drop = FALSE]; out$mu <- mu[keep] }
  out
}

## ---------------------------------------------------------------------------
## 3. Helpers
## ---------------------------------------------------------------------------
lagm <- function(A, k = 1) {
  if (k >= ncol(A)) return(matrix(NA_real_, nrow(A), ncol(A)))
  cbind(matrix(NA_real_, nrow(A), k), A[, 1:(ncol(A) - k), drop = FALSE])
}
nz <- function(x) { x[is.na(x)] <- 0; x }
## group sums returning an N-row matrix (zero rows for empty groups)
gsum <- function(X, gid, N) {
  X <- as.matrix(X)
  out <- matrix(0, N, ncol(X))
  r <- rowsum(X, gid)
  out[as.integer(rownames(r)), ] <- r
  out
}
pinv <- function(A, tol = 1e-12) {
  s <- svd(A)
  d <- ifelse(s$d > tol * max(s$d), 1 / s$d, 0)
  s$v %*% (d * t(s$u))
}
cluster_vcov <- function(X, e, gid, N, bread) {
  s <- gsum(X * e, gid, N)
  bread %*% crossprod(s) %*% bread
}
idx <- function(M) {          # row (country) and column (year) of TRUE cells, column-major
  w <- which(M, arr.ind = TRUE)
  list(g = w[, 1], t = w[, 2])
}

## ---------------------------------------------------------------------------
## 4. Estimators
## ---------------------------------------------------------------------------
est_fe <- function(d) {
  Y <- d$y; YL <- lagm(Y)
  M <- !is.na(Y) & !is.na(YL)
  ix <- idx(M); N <- nrow(Y)
  X <- cbind(YL[M], d$q[M], d$p[M], d$h[M]); yv <- Y[M]
  n <- tabulate(ix$g, N)
  Xd <- X - (gsum(X, ix$g, N) / pmax(n, 1))[ix$g, , drop = FALSE]
  yd <- yv - (gsum(yv, ix$g, N)[, 1] / pmax(n, 1))[ix$g]
  XtXi <- solve(crossprod(Xd))
  b <- drop(XtXi %*% crossprod(Xd, yd))
  e <- drop(yd - Xd %*% b)
  V <- cluster_vcov(Xd, e, ix$g, N, XtXi)
  R <- matrix(NA_real_, nrow(Y), ncol(Y)); R[cbind(ix$g, ix$t)] <- e
  list(b = b, se = sqrt(diag(V)), resid = R, nobs = length(yv))
}

first_diff_rows <- function(d) {
  Y <- d$y; Q <- d$q; Pp <- d$p; H <- d$h
  DY <- Y - lagm(Y); DYL <- lagm(Y) - lagm(Y, 2)
  DQ <- Q - lagm(Q); DP <- Pp - lagm(Pp); DH <- H - lagm(H)
  M <- !is.na(DY) & !is.na(DYL)   # y_t, y_{t-1}, y_{t-2} observed (regressors share the mask)
  ix <- idx(M)
  list(M = M, gid = ix$g, tid = ix$t, dy = DY[M],
       X = cbind(DYL[M], DQ[M], DP[M], DH[M]),
       y2 = lagm(Y, 2)[M], q2 = lagm(Q, 2)[M],
       y3 = nz(lagm(Y, 3)[M]), q3 = nz(lagm(Q, 3)[M]),
       dp = DP[M], dh = DH[M])
}

## Anderson-Hsiao: first differences; y_{t-2}, q_{t-2} instrument dy_{t-1}, dq_t
est_iv <- function(d) {
  r <- first_diff_rows(d); N <- nrow(d$y)
  Z <- cbind(r$y2, r$q2, r$dp, r$dh); X <- r$X
  ZZi <- solve(crossprod(Z))
  PzX <- Z %*% (ZZi %*% crossprod(Z, X))
  A <- solve(crossprod(PzX, X))
  b <- drop(A %*% crossprod(PzX, r$dy))
  e <- drop(r$dy - X %*% b)
  V <- cluster_vcov(PzX, e, r$gid, N, A)
  ## first-stage cluster-robust F for the two excluded instruments
  F1 <- sapply(1:2, function(j) {
    pi <- drop(ZZi %*% crossprod(Z, X[, j]))
    res <- X[, j] - drop(Z %*% pi)
    Vj <- cluster_vcov(Z, res, r$gid, N, ZZi)[1:2, 1:2]
    drop(t(pi[1:2]) %*% solve(Vj, pi[1:2])) / 2
  })
  list(b = b, se = sqrt(diag(V)), nobs = length(r$dy), F1 = F1, ninst = 2)
}

build_gmm <- function(d, system) {
  r <- first_diff_rows(d)
  Zd <- cbind(r$y2, r$y3, r$q2, r$q3, r$dp, r$dh)
  Xd <- r$X
  if (system) Xd <- cbind(Xd, 0)
  B <- list(Xd = Xd, yd = r$dy, Zd = Zd, gd = r$gid, td = r$tid)
  if (system) {
    Y <- d$y; Q <- d$q
    DYL <- lagm(Y) - lagm(Y, 2); DQL <- lagm(Q) - lagm(Q, 2)
    M <- !is.na(Y) & !is.na(lagm(Y))
    ix <- idx(M); n <- sum(M)
    B$Xl <- cbind(lagm(Y)[M], Q[M], d$p[M], d$h[M], rep(1, n))
    B$yl <- Y[M]
    B$Zl <- cbind(nz(DYL[M]), nz(DQL[M]), d$p[M], d$h[M], rep(1, n))
    B$gl <- ix$g
  }
  B
}

gmm_fit <- function(B, N, system) {
  Xd <- B$Xd; yd <- B$yd; Zd <- B$Zd; gd <- B$gd; td <- B$td
  Ld <- ncol(Zd); nd <- length(yd)
  if (system) {
    Ll <- ncol(B$Zl)
    X <- rbind(Xd, B$Xl); yv <- c(yd, B$yl)
    Z <- matrix(0, nd + length(B$yl), Ld + Ll)
    Z[1:nd, 1:Ld] <- Zd; Z[(nd + 1):nrow(Z), (Ld + 1):(Ld + Ll)] <- B$Zl
    gid <- c(gd, B$gl)
  } else {
    X <- Xd; yv <- yd; Z <- Zd; gid <- gd
  }
  ## one-step weight: H (MA(1) structure of differenced errors) for diff rows, I for levels
  Zdd <- Z[1:nd, , drop = FALSE]
  ZHZ <- 2 * crossprod(Zdd)
  o <- order(gd, td); a_ <- o[-length(o)]; b_ <- o[-1]
  adj <- gd[a_] == gd[b_] & (td[b_] - td[a_]) == 1
  a_ <- a_[adj]; b_ <- b_[adj]
  C <- crossprod(Zdd[a_, , drop = FALSE], Zdd[b_, , drop = FALSE])
  ZHZ <- ZHZ - C - t(C)
  if (system) { Zll <- Z[(nd + 1):nrow(Z), , drop = FALSE]; ZHZ <- ZHZ + crossprod(Zll) }
  A1 <- pinv(ZHZ)
  ZX <- crossprod(Z, X); Zy <- crossprod(Z, yv)
  M1 <- solve(t(ZX) %*% A1 %*% ZX)
  b1 <- drop(M1 %*% t(ZX) %*% A1 %*% Zy)
  e1 <- drop(yv - X %*% b1)
  g1 <- gsum(Z * e1, gid, N)
  S1 <- crossprod(g1)
  A2 <- pinv(S1)
  M2 <- solve(t(ZX) %*% A2 %*% ZX)
  b2 <- drop(M2 %*% t(ZX) %*% A2 %*% Zy)
  e2 <- drop(yv - X %*% b2)
  ## one-step robust variance and Windmeijer (2005) corrected two-step variance
  V1 <- M1 %*% t(ZX) %*% A1 %*% S1 %*% A1 %*% ZX %*% M1
  Ze2 <- crossprod(Z, e2)
  K <- ncol(X)
  D <- matrix(0, K, K)
  for (k in seq_len(K)) {
    ak <- gsum(Z * X[, k], gid, N)          # Z_i' x_ik
    dS <- crossprod(ak, g1) + crossprod(g1, ak)
    D[, k] <- M2 %*% t(ZX) %*% A2 %*% dS %*% A2 %*% Ze2
  }
  V2 <- M2 + D %*% M2 + M2 %*% t(D) + D %*% V1 %*% t(D)
  J <- drop(t(Ze2) %*% A2 %*% Ze2)
  df <- ncol(Z) - K
  ## iterated GMM (Hansen, Heaton & Yaron 1996)
  bi <- b2
  for (it in 1:50) {
    gi <- gsum(Z * drop(yv - X %*% bi), gid, N)
    Ai <- pinv(crossprod(gi))
    bn <- drop(solve(t(ZX) %*% Ai %*% ZX, t(ZX) %*% Ai %*% Zy))
    conv <- max(abs(bn - bi)) < 1e-10
    bi <- bn
    if (conv) break
  }
  gi <- gsum(Z * drop(yv - X %*% bi), gid, N); gs <- colSums(gi)
  Jit <- drop(t(gs) %*% pinv(crossprod(gi)) %*% gs)
  ## Arellano-Bond AR(1), AR(2) tests on differenced residuals
  ed <- drop(yd - Xd %*% b2)
  Tn <- max(td)
  E <- matrix(NA_real_, N, TT); E[cbind(gd, td)] <- ed
  XS <- lapply(seq_len(K), function(k) { A <- matrix(NA_real_, N, TT); A[cbind(gd, td)] <- Xd[, k]; A })
  g2 <- gsum(Z * e2, gid, N)
  ar <- sapply(1:2, function(j) {
    EL <- lagm(E, j)
    m <- !is.na(E) & !is.na(EL)
    ix <- idx(m)
    es <- E[m]; ws <- EL[m]
    Xs <- sapply(XS, function(A) A[m])
    s_i <- gsum(ws * es, ix$g, N)[, 1]
    num <- sum(s_i)
    wX <- drop(ws %*% Xs)
    term2 <- drop(t(wX) %*% M2 %*% t(ZX) %*% A2 %*% crossprod(g2, s_i))
    vv <- sum(s_i^2) - 2 * term2 + drop(t(wX) %*% V2 %*% wX)
    if (vv > 0) num / sqrt(vv) else NA_real_
  })
  list(b = b2[1:4], se = sqrt(diag(V2))[1:4], J = J, Jdf = df,
       Jp = pchisq(J, df, lower.tail = FALSE), Jit = Jit,
       Jit_p = pchisq(Jit, df, lower.tail = FALSE), ninst = ncol(Z),
       ar1 = ar[1], ar2 = ar[2], nobs = length(yv), nobs_d = nd)
}

est_dgmm <- function(d) gmm_fit(build_gmm(d, FALSE), nrow(d$y), FALSE)
est_sgmm <- function(d) gmm_fit(build_gmm(d, TRUE), nrow(d$y), TRUE)

## ---------------------------------------------------------------------------
## 5. Diagnostics
## ---------------------------------------------------------------------------
pesaran_cd <- function(E, min_overlap = 3) {
  M <- (!is.na(E)) * 1; E0 <- nz(E)
  n <- M %*% t(M)
  sx <- E0 %*% t(M); sy <- M %*% t(E0); sxy <- E0 %*% t(E0)
  sxx <- (E0^2) %*% t(M); syy <- M %*% t(E0^2)
  rho <- suppressWarnings((sxy - sx * sy / n) / sqrt((sxx - sx^2 / n) * (syy - sy^2 / n)))
  iu <- upper.tri(n)
  ok <- iu & n >= min_overlap & is.finite(rho)
  Nn <- nrow(E)
  cd <- sqrt(2 / (Nn * (Nn - 1))) * sum(sqrt(n[ok]) * rho[ok])
  c(cd = cd, p = 2 * pnorm(-abs(cd)), abs_rho = mean(abs(rho[ok])))
}

## Wooldridge (2010) test: regress FD residuals on their lag, H0: coefficient = -0.5
wooldridge_ar <- function(d) {
  r <- first_diff_rows(d); N <- nrow(d$y)
  b <- qr.coef(qr(r$X), r$dy)
  e <- drop(r$dy - r$X %*% b)
  E <- matrix(NA_real_, N, ncol(d$y)); E[cbind(r$gid, r$tid)] <- e
  EL <- lagm(E); m <- !is.na(E) & !is.na(EL)
  ix <- idx(m); x <- EL[m]; yv <- E[m]
  rho <- sum(x * yv) / sum(x * x)
  s <- gsum(x * (yv - rho * x), ix$g, N)[, 1]
  se <- sqrt(sum(s^2)) / sum(x * x)
  Fs <- ((rho + 0.5) / se)^2
  c(F = Fs, p = pf(Fs, 1, N - 1, lower.tail = FALSE), rho = rho)
}

## Sample correlation between each instrument and the TRUE transformed disturbance
oracle_moments <- function(d) {
  eps <- d$eps; mu <- d$mu
  r <- first_diff_rows(d)
  deps <- (eps - lagm(eps))[r$M]
  out <- list()
  for (nm in c("y2", "y3", "q2", "q3")) {
    zz <- r[[nm]]; m <- zz != 0
    out[[paste0("diff_", nm)]] <- cor(zz[m], deps[m])
  }
  Y <- d$y; Q <- d$q
  M <- !is.na(Y) & !is.na(lagm(Y))
  lev <- (matrix(mu, nrow(Y), ncol(Y)) + eps)[M]
  for (nm in c("lev_dy1", "lev_dq1")) {
    A <- if (nm == "lev_dy1") lagm(Y) - lagm(Y, 2) else lagm(Q) - lagm(Q, 2)
    zz <- A[M]; m <- !is.na(zz)
    out[[nm]] <- cor(zz[m], lev[m])
  }
  Mq <- !is.na(Q)
  Qd <- Q - rowMeans(Q, na.rm = TRUE)
  out$corr_q_eps <- cor(Qd[Mq], eps[Mq])     # endogeneity of q (reference)
  out
}

ESTIMATORS <- list(FE = est_fe, IV = est_iv, DGMM = est_dgmm, SGMM = est_sgmm)

one_rep <- function(k, r) {
  lam <- LAMBDAS[k + 1]
  set.seed(rep_seed(k, r))
  d <- simulate_panel(lam, return_truth = TRUE)
  rec <- list(lam = lam, rep = r, Nret = nrow(d$y), nobs = sum(!is.na(d$y)))
  for (nm in names(ESTIMATORS)) {
    res <- tryCatch(ESTIMATORS[[nm]](d), error = function(e) NULL)
    if (is.null(res) || !all(is.finite(res$b)) || !all(is.finite(res$se))) {
      rec[[nm]] <- NA; next
    }
    if (nm == "FE") {
      cd <- pesaran_cd(res$resid); wo <- wooldridge_ar(d)
      res$cd <- cd[["cd"]]; res$cd_p <- cd[["p"]]; res$abs_rho <- cd[["abs_rho"]]
      res$wool_F <- wo[["F"]]; res$wool_p <- wo[["p"]]; res$wool_rho <- wo[["rho"]]
    }
    res$resid <- NULL
    rec[[nm]] <- res
  }
  rec$oracle <- oracle_moments(d)
  rec
}

## ---------------------------------------------------------------------------
## 6. Illustrative replication (Section 3.2 of the manuscript)
## ---------------------------------------------------------------------------
illustrative <- function(outdir, lam = 0.9, r = 0) {
  k <- match(lam, LAMBDAS) - 1
  set.seed(rep_seed(k, r))
  d <- simulate_panel(lam, return_truth = TRUE)
  res <- lapply(ESTIMATORS, function(f) f(d))
  lv <- list(GDP = exp(d$y), CO2 = exp(d$q), Precip = exp(d$p), Temp = d$h)
  desc <- lapply(lv, function(v) list(min = min(v, na.rm = TRUE), median = median(v, na.rm = TRUE),
                                      mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
                                      max = max(v, na.rm = TRUE)))
  ar1 <- lapply(list(GDP = d$y, CO2 = d$q, Precip = d$p, Temp = d$h), function(A) {
    AL <- lagm(A); m <- !is.na(A) & !is.na(AL); ix <- idx(m)
    X <- cbind(1, AL[m]); yv <- A[m]
    XtXi <- solve(crossprod(X)); b <- drop(XtXi %*% crossprod(X, yv))
    se <- sqrt(diag(cluster_vcov(X, drop(yv - X %*% b), ix$g, nrow(A), XtXi)))
    list(b = b[2], se = se[2])
  })
  cd <- pesaran_cd(res$FE$resid); wo <- wooldridge_ar(d)
  est <- lapply(res, function(x) { x$resid <- NULL; x })
  out <- list(lam = lam, rep = r, Nret = nrow(d$y), nobs = sum(!is.na(d$y)),
              desc = desc, ar1 = ar1, cd = unname(cd), wooldridge = unname(wo), est = est)
  write_json(out, file.path(outdir, "illustrative.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
  invisible(out)
}

## ---------------------------------------------------------------------------
## 7. Driver
## ---------------------------------------------------------------------------
run_study <- function(reps = 1000, outdir = "results",
                      cores = max(1, parallel::detectCores())) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  jobs <- expand.grid(r = 0:(reps - 1), k = seq_along(LAMBDAS) - 1)
  t0 <- Sys.time()
  run <- function(i) one_rep(jobs$k[i], jobs$r[i])
  recs <- if (.Platform$OS.type == "unix" && cores > 1) {
    parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores, mc.preschedule = TRUE)
  } else lapply(seq_len(nrow(jobs)), run)
  con <- file(file.path(outdir, "replications.jsonl"), "w")
  for (rec in recs) writeLines(toJSON(rec, auto_unbox = TRUE, digits = NA, na = "null"), con)
  close(con)
  write_json(list(master_seed = MASTER_SEED, seed_rule = "set.seed(master_seed + 100000*k + r)",
                  rng = RNGkind(), lambdas = LAMBDAS, N0 = N0, T = TT, burn = BURN, reps = reps,
                  true = as.list(TRUE_THETA), params = P, R_version = R.version.string),
             file.path(outdir, "design.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
  message(sprintf("%d replications in %.0f s", nrow(jobs), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(recs)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  getarg <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1] }
  outdir <- getarg("--out", "results")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if ("--illustrative" %in% args) {
    illustrative(outdir)
  } else {
    run_study(as.integer(getarg("--reps", "1000")), outdir,
              as.integer(getarg("--cores", parallel::detectCores())))
  }
}
