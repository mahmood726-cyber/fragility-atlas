# Extended validation: DL + REML + FE across 10 reviews
# Run: Rscript tests/validate_all_estimators.R

library(metafor)

specs <- read.csv("C:/FragilityAtlas/data/output/fragility_atlas_specifications.csv",
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8")
results <- read.csv("C:/FragilityAtlas/data/output/fragility_atlas_results.csv",
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8")
pairwise_dir <- "C:/Models/Pairwise70/data"

# Same 10 reviews as before
set.seed(42)
eligible <- results[results$k >= 5 & results$scale == "ratio", ]
sample_ids <- eligible$review_id[sample(nrow(eligible), min(10, nrow(eligible)))]

cat("============================================================\n")
cat("Extended R Validation: DL + REML + FE (7 estimators)\n")
cat("============================================================\n\n")

total_checks <- 0
pass_checks <- 0

for (rid in sample_ids) {
  py_result <- results[results$review_id == rid, ]
  py_analysis <- py_result$analysis_name[1]
  py_k <- py_result$k[1]

  rda_files <- list.files(pairwise_dir, pattern = paste0("^", rid, "_"), full.names = TRUE)
  if (length(rda_files) == 0) next

  env <- new.env()
  load(rda_files[1], envir = env)
  df <- get(ls(env)[1], envir = env)

  if (!"Analysis.name" %in% names(df)) next
  matching <- df[df$Analysis.name == py_analysis, ]
  if (nrow(matching) == 0) next

  candidates <- unique(matching[, c("Analysis.group", "Analysis.number")])
  best_diff <- Inf; primary <- NULL
  for (i in seq_len(nrow(candidates))) {
    cand <- matching[matching$Analysis.group == candidates$Analysis.group[i] &
                     matching$Analysis.number == candidates$Analysis.number[i], ]
    if (abs(nrow(cand) - py_k) < best_diff) {
      best_diff <- abs(nrow(cand) - py_k); primary <- cand
    }
  }
  if (is.null(primary)) next

  primary <- primary[!is.na(primary$Mean) & !is.na(primary$CI.start) & !is.na(primary$CI.end) &
                     primary$Mean > 0 & primary$CI.start > 0 & primary$CI.end > 0, ]
  if (nrow(primary) < 3) next
  yi <- log(primary$Mean)
  sei <- (log(primary$CI.end) - log(primary$CI.start)) / (2 * qnorm(0.975))
  valid <- sei > 0 & is.finite(yi) & is.finite(sei)
  yi <- yi[valid]; sei <- sei[valid]
  if (length(yi) < 3) next

  cat(sprintf("  %s (k=%d):\n", rid, length(yi)))

  # Test DL, REML, FE, PM, HS, HE, SJ
  methods <- c("DL", "REML", "FE", "PM", "HS", "HE", "SJ")
  for (m in methods) {
    # Get Python estimate
    py_row <- specs[specs$review_id == rid & specs$estimator == m &
                    specs$ci_method == "Wald" & specs$bias_correction == "none" &
                    specs$leave_out == "", ]
    if (nrow(py_row) == 0) next

    tryCatch({
      fit <- rma(yi = yi, sei = sei, method = m)
      r_theta <- as.numeric(coef(fit))
      py_theta <- py_row$theta[1]
      d <- abs(r_theta - py_theta)
      total_checks <<- total_checks + 1
      status <- ifelse(d < 1e-3, "PASS", "FAIL")
      if (status == "PASS") pass_checks <<- pass_checks + 1
      cat(sprintf("    %4s: Py=%.6f R=%.6f diff=%.2e [%s]\n", m, py_theta, r_theta, d, status))
    }, error = function(e) {
      cat(sprintf("    %4s: R error - %s\n", m, conditionMessage(e)))
    })
  }
}

cat(sprintf("\n  TOTAL: %d/%d checks passed (%.1f%%)\n", pass_checks, total_checks,
            100 * pass_checks / max(1, total_checks)))
