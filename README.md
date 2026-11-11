# Supplementary Material S1: R replication code

Reproduces every table and figure in Opoku, Adu, Emmanuel & Owusu,
"Evaluating Estimator Performance in Unbalanced Dynamic Climate–Economy Panels".

Requirements: R >= 4.0 and the `jsonlite` package (`install.packages("jsonlite")`).
Everything else is base R. Tested with R 4.3.3.

From a terminal, in this folder:

    Rscript climate_panel_mc.R --reps 1000 --out results       # 4 x 1,000 replications (~2 min on 2 cores)
    Rscript climate_panel_mc.R --illustrative --out results    # illustrative replication (lambda = 0.9, r = 0)
    Rscript summarize.R results                                # performance + validation summary
    Rscript make_tables.R results manuscript_template.tex manuscript.tex

Or inside R / RStudio:

    source("climate_panel_mc.R")
    run_study(reps = 1000, outdir = "results")   # on Windows this runs on one core (slower)
    illustrative("results")

Seeds: replication r (0-999) of setting k (lambda = 0, 0.3, 0.6, 0.9 -> k = 0..3) is generated
after `set.seed(20261030 + 100000*k + r)` with R's default Mersenne-Twister generator, so results
do not depend on the number of cores.

Files:
- `climate_panel_mc.R`  data-generating process, FE, Anderson-Hsiao IV, Difference- and
  System-GMM (collapsed lags 2-3, two-step, Windmeijer SEs, iterated Hansen), AR(1)/AR(2),
  Pesaran CD, Wooldridge tests and direct moment checks.
- `summarize.R`         bias, RB, SD, RMSE, MCSE, median bias, MAE, coverage, validation rates.
- `make_tables.R`       writes all LaTeX tables and figure data into the manuscript template.
- `results/replications.jsonl`  one record per replication (estimates, SEs, instrument counts,
  test statistics, moment checks); `results/summary.json`, `results/illustrative.json`.
