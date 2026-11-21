# run_all.R
# ---------------------------------------------------------------------------
# Run the whole pipeline end to end on the synthetic cohort:
#
#   Rscript run_all.R              # default n_boot from config.yml (fast)
#   Rscript run_all.R --n-boot 200 # full bootstrap used for the README figures
#
# Each step is an ordinary script that can also be run on its own.
# ---------------------------------------------------------------------------

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- if (length(this_file)) dirname(normalizePath(this_file)) else getwd()
setwd(root)

steps <- c("R/00_simulate_cohort.R", "R/01_linkage_checks.R", "R/02_features.R",
           "R/03_models.R", "R/04_report.R")
extra_args <- commandArgs(trailingOnly = TRUE)

t_all <- Sys.time()
timings <- data.frame(step = steps, seconds = NA_real_)
for (i in seq_along(steps)) {
  cat("\n==================== ", steps[i], " ====================\n", sep = "")
  t0 <- Sys.time()
  # Each step runs in a fresh R process so that nothing leaks between steps.
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c(steps[i], if (grepl("03_models", steps[i])) extra_args),
                    stdout = "", stderr = "")
  if (status != 0) stop("Step failed: ", steps[i])
  timings$seconds[i] <- round(as.numeric(Sys.time() - t0, units = "secs"), 1)
}

cat("\n==================== timings ====================\n")
print(timings, row.names = FALSE)
cat(sprintf("Total: %.1f s\n", as.numeric(Sys.time() - t_all, units = "secs")))
write.csv(timings, "output/run_timings.csv", row.names = FALSE)
