# Every injected fault must be caught by at least one rule, record by record.
fault_to_issue <- list(
  duplicate_birth_record       = "duplicate_birth_record",
  implausible_birth_weight     = "implausible_birth_weight",
  missing_gestational_age      = "missing_gestational_age",
  implausible_gestational_age  = "implausible_gestational_age",
  implausible_maternal_age     = "implausible_maternal_age",
  implausible_sibling_interval = "implausible_sibling_interval",
  orphan_visit                 = c("gp_orphan_visit", "ed_orphan_visit"),
  visit_before_birth           = c("gp_visit_before_birth", "ed_visit_before_birth"),
  visit_after_extract          = c("gp_visit_after_extract", "ed_visit_after_extract"),
  duplicate_visit_row          = c("gp_duplicate_visit_row", "ed_duplicate_visit_row")
)

key <- fread(file.path(syn, "injected_faults.csv"))
detail <- fread(file.path(out, "data_quality_detail.csv"))
log <- fread(file.path(out, "data_quality_log.csv"))

test_that("the log has the required columns", {
  expect_true(all(c("issue", "source_table", "n_rows", "rule", "action") %in% names(log)))
})

for (fault_class in unique(key$fault)) {
  test_that(sprintf("fault class '%s' is fully detected", fault_class), {
    expect_true(fault_class %in% names(fault_to_issue), info = "fault class has a mapped rule")
    injected <- unique(key[fault == fault_class, record_id])
    caught <- unique(detail[issue %in% fault_to_issue[[fault_class]], record_id])
    missed <- setdiff(injected, caught)
    expect_length(missed, 0)
  })
}

test_that("the clean tables satisfy the rules they were built from", {
  con <- open_db(cfg); on.exit(dbDisconnect(con))
  n_distinct <- dbGetQuery(con, "SELECT COUNT(DISTINCT child_id) AS n FROM births")$n
  expect_equal(dbGetQuery(con, "SELECT COUNT(*) AS n FROM births_clean")$n, n_distinct)
  expect_equal(dbGetQuery(con, "SELECT COUNT(*) AS n FROM births_clean WHERE birth_weight_g < 300 OR birth_weight_g > 6500")$n, 0)
  expect_equal(nrow(dbGetQuery(con, "SELECT gp_record_id FROM gp_visits_clean GROUP BY gp_record_id HAVING COUNT(*) > 1")), 0)
  expect_equal(dbGetQuery(con, "SELECT COUNT(*) AS n FROM ed_presentations_clean e LEFT JOIN births_clean b USING (child_id) WHERE b.child_id IS NULL")$n, 0)
})
