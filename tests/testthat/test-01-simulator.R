test_that("the simulator is reproducible from its seed", {
  # Re-run the simulator in a throw-away copy of the project and compare file hashes.
  tmp <- file.path(tempdir(), "sim_repro")
  unlink(tmp, recursive = TRUE)
  dir.create(file.path(tmp, "R"), recursive = TRUE)
  file.copy(file.path(root, "config.yml"), tmp)
  file.copy(file.path(root, "R", c("utils.R", "00_simulate_cohort.R")), file.path(tmp, "R"))
  # The script finds its project root from the working directory, so run it from the copy.
  status <- system(sprintf("cd '%s' && '%s' R/00_simulate_cohort.R", tmp, rscript),
                   ignore.stdout = TRUE, ignore.stderr = TRUE)
  expect_equal(status, 0L)
  for (f in c("births.csv", "gp_visits.csv", "ed_presentations.csv", "injected_faults.csv")) {
    expect_identical(unname(tools::md5sum(file.path(tmp, "data", "synthetic", f))),
                     unname(tools::md5sum(file.path(syn, f))), info = f)
  }
})

test_that("the synthetic cohort has the intended shape", {
  births <- fread(file.path(syn, "births.csv"))
  key <- fread(file.path(syn, "injected_faults.csv"))
  expect_gte(length(unique(births$child_id)), cfg$simulation$n_births)
  expect_true(all(c("maternal_asthma", "maternal_smoking_pregnancy", "seifa_quintile", "season_of_birth") %in% names(births)))
  expect_gte(length(unique(key$fault)), 8)
  feat <- fread(file.path(out, "feature_matrix.csv"))
  prev <- mean(feat$outcome)
  expect_gt(prev, 0.12)
  expect_lt(prev, 0.22)
})
