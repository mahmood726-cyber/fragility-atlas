# Proper R cross-validation: match Python's primary outcome selection
# Compares Python DL estimates from specifications CSV against metafor
# Run: Rscript tests/validate_proper.R

library(metafor)

pairwise_dir <- "C:/Models/Pairwise70/data"

# Read Python DL results (DL/Wald/none/no-LOO)
specs <- read.csv("C:/FragilityAtlas/data/output/fragility_atlas_specifications.csv",
                   stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Also read results for analysis selection info
results <- read.csv("C:/FragilityAtlas/data/output/fragility_atlas_results.csv",
                     stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Get Python DL/Wald/none/no-LOO specs
dl_specs <- specs[specs$estimator == "DL" &
                  specs$ci_method == "Wald" &
                  specs$bias_correction == "none" &
                  specs$leave_out == "", ]

# Select 10 reviews: seeded, k>=5, ratio scale
set.seed(42)
eligible <- results[results$k >= 5 & results$scale == "ratio", ]
sample_ids <- eligible$review_id[sample(nrow(eligible), min(10, nrow(eligible)))]

cat("============================================================\n")
cat("R Cross-Validation: Python DL vs metafor (same analysis)\n")
cat("============================================================\n\n")

pass_count <- 0
fail_count <- 0
validation_rows <- list()

for (rid in sample_ids) {
  # Get Python's chosen analysis name
  py_result <- results[results$review_id == rid, ]
  py_analysis <- py_result$analysis_name[1]
  py_k <- py_result$k[1]

  # Get Python's DL theta and tau2
  py_dl <- dl_specs[dl_specs$review_id == rid, ]
  if (nrow(py_dl) == 0) {
    cat(sprintf("  SKIP %s: no Python DL spec found\n", rid))
    next
  }
  py_theta <- py_dl$theta[1]
  py_tau2 <- py_dl$tau2[1]
  py_scale <- py_result$scale[1]

  # Load RDA
  rda_files <- list.files(pairwise_dir, pattern = paste0("^", rid, "_"), full.names = TRUE)
  if (length(rda_files) == 0) {
    cat(sprintf("  SKIP %s: no RDA file\n", rid))
    next
  }

  env <- new.env()
  load(rda_files[1], envir = env)
  df <- get(ls(env)[1], envir = env)

  # Match Python's analysis: find the (group, number) pair whose name matches
  # and whose filtered k is closest to Python's k
  if (!"Analysis.name" %in% names(df)) {
    cat(sprintf("  SKIP %s: no Analysis.name column\n", rid))
    next
  }

  matching <- df[df$Analysis.name == py_analysis, ]
  if (nrow(matching) == 0) {
    cat(sprintf("  SKIP %s: analysis '%s' not found in RDA\n", rid, py_analysis))
    next
  }

  # If multiple (group, number) pairs share the same name, pick the one closest to py_k
  candidates <- unique(matching[, c("Analysis.group", "Analysis.number")])
  best_diff <- Inf
  primary <- NULL
  for (i in seq_len(nrow(candidates))) {
    cand <- matching[matching$Analysis.group == candidates$Analysis.group[i] &
                     matching$Analysis.number == candidates$Analysis.number[i], ]
    this_diff <- abs(nrow(cand) - py_k)
    if (this_diff < best_diff) {
      best_diff <- this_diff
      primary <- cand
    }
  }

  if (is.null(primary) || nrow(primary) == 0) {
    cat(sprintf("  SKIP %s: no matching analysis found\n", rid))
    next
  }

  # Compute yi/sei following Python's logic
  if (py_scale == "ratio") {
    # Log-transform Mean, CI
    primary <- primary[!is.na(primary$Mean) & !is.na(primary$CI.start) &
                       !is.na(primary$CI.end) &
                       primary$Mean > 0 & primary$CI.start > 0 & primary$CI.end > 0, ]
    if (nrow(primary) < 3) {
      cat(sprintf("  SKIP %s: too few valid ratio studies (%d)\n", rid, nrow(primary)))
      next
    }
    yi <- log(primary$Mean)
    sei <- (log(primary$CI.end) - log(primary$CI.start)) / (2 * qnorm(0.975))
  } else {
    # Difference scale: yi = Mean, sei from CI width
    primary <- primary[!is.na(primary$Mean) & !is.na(primary$CI.start) &
                       !is.na(primary$CI.end), ]
    if (nrow(primary) < 3) {
      cat(sprintf("  SKIP %s: too few valid diff studies (%d)\n", rid, nrow(primary)))
      next
    }
    yi <- primary$Mean
    sei <- (primary$CI.end - primary$CI.start) / (2 * qnorm(0.975))
  }

  # Filter invalid SE
  valid <- sei > 0 & is.finite(yi) & is.finite(sei)
  yi <- yi[valid]
  sei <- sei[valid]

  if (length(yi) < 3) {
    cat(sprintf("  SKIP %s: too few after SE filter (%d)\n", rid, length(yi)))
    next
  }

  # Run metafor DL
  tryCatch({
    fit_dl <- rma(yi = yi, sei = sei, method = "DL")
    r_theta <- as.numeric(coef(fit_dl))
    r_tau2 <- fit_dl$tau2

    diff_theta <- abs(py_theta - r_theta)
    diff_tau2 <- abs(py_tau2 - r_tau2)
    status <- ifelse(diff_theta < 1e-4, "PASS", ifelse(diff_theta < 0.01, "WARN", "FAIL"))
    if (status == "PASS") pass_count <<- pass_count + 1 else fail_count <<- fail_count + 1

    cat(sprintf("  %s (k=%d/%d): Py=%.6f R=%.6f diff=%.2e tau2_diff=%.2e [%s]\n",
                rid, length(yi), py_k, py_theta, r_theta, diff_theta, diff_tau2, status))

    validation_rows[[length(validation_rows) + 1]] <- data.frame(
      review_id = rid,
      k_r = length(yi),
      k_py = py_k,
      py_dl_theta = py_theta,
      r_dl_theta = r_theta,
      theta_diff = diff_theta,
      py_dl_tau2 = py_tau2,
      r_dl_tau2 = r_tau2,
      tau2_diff = diff_tau2,
      status = status,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  ERROR %s: %s\n", rid, conditionMessage(e)))
  })
}

cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")
cat(sprintf("  Passed: %d, Failed: %d, Total: %d\n", pass_count, fail_count, pass_count + fail_count))
cat(sprintf("  Pass rate: %.1f%%\n", 100 * pass_count / max(1, pass_count + fail_count)))

if (length(validation_rows) > 0) {
  out <- do.call(rbind, validation_rows)
  write.csv(out, "C:/FragilityAtlas/data/output/r_crossvalidation_proper.csv", row.names = FALSE)
  cat("  Saved to data/output/r_crossvalidation_proper.csv\n")
}
