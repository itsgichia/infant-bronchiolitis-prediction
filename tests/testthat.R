# Run the test suite:  Rscript tests/testthat.R
# The tests read the outputs of a completed pipeline run (Rscript run_all.R).
root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config.yml")) && dirname(d) != d) d <- dirname(d); d })
setwd(root)
library(testthat)
res <- test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
