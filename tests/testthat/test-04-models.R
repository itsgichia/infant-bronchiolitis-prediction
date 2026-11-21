cmp <- fread(file.path(out, "model_comparison.csv"))

test_that("all four model families are reported and discriminate on synthetic data", {
  expect_setequal(cmp$model, c("logistic", "elastic_net", "random_forest", "xgboost"))
  expect_true(all(cmp$auc_test > 0.6))
  expect_true(all(cmp$auc_optimism_corrected > 0.6))
  expect_true(all(cmp$auc_optimism >= 0))                   # optimism cannot be negative on average
  expect_true(all(cmp$auc_optimism_corrected <= cmp$auc_apparent))
  expect_true(all(cmp$brier > 0 & cmp$brier < 0.25))
})

test_that("net benefit matches a reference implementation", {
  set.seed(1)
  y <- rbinom(500, 1, 0.2); p <- plogis(rnorm(500, -1.5 + 1.5 * y, 1))
  thr <- c(0.05, 0.1, 0.2, 0.3, 0.5)
  nb <- net_benefit(y, p, thr)
  # Hand-computed definition
  nb_ref <- vapply(thr, function(t) {
    tp <- sum(p >= t & y == 1); fp <- sum(p >= t & y == 0)
    tp / 500 - fp / 500 * t / (1 - t)
  }, numeric(1))
  expect_equal(nb, nb_ref)
  if (requireNamespace("dcurves", quietly = TRUE) && requireNamespace("tibble", quietly = TRUE)) {
    d <- data.frame(y = y, p = p)
    dc <- dcurves::dca(y ~ p, data = d, thresholds = thr)
    dc_tbl <- as.data.frame(tibble::as_tibble(dc))
    ref <- dc_tbl[dc_tbl$variable == "p", c("threshold", "net_benefit")]
    expect_equal(ref$net_benefit[order(ref$threshold)], nb, tolerance = 1e-8)
    all_ref <- dc_tbl[dc_tbl$variable == "all", c("threshold", "net_benefit")]
    expect_equal(all_ref$net_benefit[order(all_ref$threshold)], net_benefit_all(y, thr), tolerance = 1e-8)
  }
})

test_that("figures were written", {
  for (f in c("roc_curves.png", "calibration.png", "decision_curves.png",
              "shap_beeswarm.png", "feature_importance.png", "bootstrap_optimism.png")) {
    expect_true(file.exists(file.path(cfg$paths$figures_dir, f)), info = f)
  }
})
