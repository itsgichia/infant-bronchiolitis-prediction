feat <- fread(file.path(out, "feature_matrix.csv"))
dict <- fread(file.path(out, "feature_dictionary.csv"))
predictors <- dict[role == "predictor", feature]

test_that("feature matrix and dictionary agree", {
  expect_setequal(names(feat), dict$feature)
  expect_gte(length(predictors), 40)
  expect_false(anyNA(feat))
})

test_that("no predictor leaks the outcome", {
  # Nothing derived from outcome codes or from after the landmark may be a predictor.
  leaky <- grepl("outcome|first_outcome|early|resp|J21|R06|bronch|wheeze", predictors, ignore.case = TRUE)
  expect_false(any(leaky), info = paste(predictors[leaky], collapse = ", "))
  # No single predictor should separate the outcome almost perfectly.
  uni_auc <- vapply(predictors, function(f) {
    a <- fast_auc(feat$outcome, feat[[f]]); max(a, 1 - a)
  }, numeric(1))
  expect_lt(max(uni_auc), 0.85)
})

test_that("utilisation features only use encounters up to the landmark", {
  # Recompute one feature from the raw clean tables and compare.
  con <- open_db(cfg); on.exit(dbDisconnect(con))
  chk <- as.data.table(dbGetQuery(con, sprintf(
    "SELECT b.child_id, COUNT(g.gp_record_id) AS n
     FROM births_clean b LEFT JOIN gp_visits_clean g
       ON g.child_id = b.child_id
      AND julianday(g.visit_date) - julianday(b.dob) <= %d
      AND NOT (g.reason_code LIKE 'J21%%' OR g.reason_code = 'R06.2')
     GROUP BY b.child_id", cfg$design$landmark_days)))
  m <- merge(feat[, .(child_id, gp_visits_lm)], chk, by = "child_id")
  expect_equal(m$gp_visits_lm, m$n)
})
