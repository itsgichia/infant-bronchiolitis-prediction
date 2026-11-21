# Shared setup for the tests: locate the project root and load the helpers.
root <- local({ d <- getwd(); while (!file.exists(file.path(d, "config.yml")) && dirname(d) != d) d <- dirname(d); d })
source(file.path(root, "R", "utils.R"))
cfg <- load_config(root)
out <- cfg$paths$output_dir
syn <- cfg$paths$synthetic_dir
rscript <- file.path(R.home("bin"), "Rscript")
