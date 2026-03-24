# R cross-validation: compare Python Fragility Atlas results against metafor
# Run: Rscript tests/validate_vs_metafor.R
#
# Reads 10 randomly selected reviews from Pairwise70 RDA files,
# computes DL/REML/FE pooled estimates using metafor, and compares
# against Python pipeline results.

library(metafor)

pairwise_dir <- "C:/Models/Pairwise70/data"
python_results <- read.csv("C:/FragilityAtlas/data/output/fragility_atlas_results.csv",
                           stringsAsFactors = FALSE)

# Select 10 reviews with varying k (seeded for reproducibility)
set.seed(42)
eligible <- python_results[python_results$k >= 5 & python_results$scale == "ratio", ]
sample_ids <- eligible$review_id[sample(nrow(eligible), min(10, nrow(eligible)))]

cat("=" , rep("=", 58), "\n", sep = "")
cat("R Cross-Validation: Python vs metafor\n")
cat("=" , rep("=", 58), "\n\n", sep = "")

results <- data.frame()

for (rid in sample_ids) {
  # Find the RDA file
  rda_files <- list.files(pairwise_dir, pattern = paste0("^", rid, "_"), full.names = TRUE)
  if (length(rda_files) == 0) {
    cat("  SKIP:", rid, "- RDA file not found\n")
    next
  }

  # Load RDA
  load(rda_files[1])
  df_name <- ls()[!ls() %in% c("pairwise_dir", "python_results", "eligible",
                                "sample_ids", "results", "rid", "rda_files")]
  # The RDA contains one dataframe - get it
  env <- new.env()
  load(rda_files[1], envir = env)
  df <- get(ls(env)[1], envir = env)

  # Select primary analysis (largest k among binary)
  analyses <- aggregate(Study ~ Analysis.group + Analysis.number, data = df, FUN = length)
  names(analyses)[3] <- "k"
  best <- analyses[which.max(analyses$k), ]
  primary <- df[df$Analysis.group == best$Analysis.group &
                df$Analysis.number == best$Analysis.number, ]

  # Back-calculate yi and sei (log-transform for ratio scale)
  primary <- primary[!is.na(primary$Mean) & !is.na(primary$CI.start) &
                     !is.na(primary$CI.end) & primary$Mean > 0 &
                     primary$CI.start > 0 & primary$CI.end > 0, ]

  if (nrow(primary) < 3) {
    cat("  SKIP:", rid, "- too few valid studies\n")
    next
  }

  yi <- log(primary$Mean)
  sei <- (log(primary$CI.end) - log(primary$CI.start)) / (2 * 1.96)
  sei[sei <= 0] <- NA
  valid <- !is.na(yi) & !is.na(sei) & sei > 0
  yi <- yi[valid]
  sei <- sei[valid]

  if (length(yi) < 3) next

  # Get Python results for this review
  py_row <- python_results[python_results$review_id == rid, ]
  if (nrow(py_row) == 0) next

  # Run metafor
  tryCatch({
    # DL
    fit_dl <- rma(yi = yi, sei = sei, method = "DL")
    # REML
    fit_reml <- rma(yi = yi, sei = sei, method = "REML")
    # FE
    fit_fe <- rma(yi = yi, sei = sei, method = "FE")

    cat(sprintf("  %s (k=%d):\n", rid, length(yi)))
    cat(sprintf("    DL   : R theta=%.4f, tau2=%.4f\n", coef(fit_dl), fit_dl$tau2))
    cat(sprintf("    REML : R theta=%.4f, tau2=%.4f\n", coef(fit_reml), fit_reml$tau2))
    cat(sprintf("    FE   : R theta=%.4f\n", coef(fit_fe)))
    cat(sprintf("    Py   : median_theta=%.4f (from pipeline CSV)\n", py_row$median_theta))

    results <- rbind(results, data.frame(
      review_id = rid,
      k = length(yi),
      r_dl_theta = as.numeric(coef(fit_dl)),
      r_dl_tau2 = fit_dl$tau2,
      r_reml_theta = as.numeric(coef(fit_reml)),
      r_reml_tau2 = fit_reml$tau2,
      r_fe_theta = as.numeric(coef(fit_fe)),
      py_median_theta = py_row$median_theta,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat("  ERROR:", rid, "-", conditionMessage(e), "\n")
  })
}

cat("\n")
cat("=" , rep("=", 58), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 58), "\n\n", sep = "")

if (nrow(results) > 0) {
  # Compare Python median theta (which is the median across all specs) to R DL theta
  diffs <- abs(results$r_dl_theta - results$py_median_theta)
  cat(sprintf("  Reviews validated: %d\n", nrow(results)))
  cat(sprintf("  Mean |R_DL - Py_median|: %.4f\n", mean(diffs)))
  cat(sprintf("  Max  |R_DL - Py_median|: %.4f\n", max(diffs)))
  cat(sprintf("  All within 0.1: %s\n", ifelse(all(diffs < 0.1), "YES", "NO")))

  # Save
  write.csv(results, "C:/FragilityAtlas/data/output/r_crossvalidation.csv", row.names = FALSE)
  cat("\n  Saved to data/output/r_crossvalidation.csv\n")
} else {
  cat("  No reviews validated.\n")
}
