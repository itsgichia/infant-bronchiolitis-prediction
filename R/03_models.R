# 03_models.R
# ---------------------------------------------------------------------------
# Benchmark four model families for predicting wheeze/bronchiolitis in the
# first year of life, with internal validation designed around the questions
# a clinician would ask: does it discriminate, is it calibrated, and would
# acting on it do more good than harm?
#
#   1. Logistic regression on 13 pre-specified clinical predictors (baseline)
#   2. Elastic-net logistic regression on all predictors (glmnet)
#   3. Random forest (ranger)
#   4. Gradient-boosted trees (xgboost)
#
# Validation
#   * stratified 70/30 hold-out
#   * Harrell's bootstrap optimism correction of AUC and calibration slope
#     on the training set (n_boot resamples, see config.yml)
#   * on the hold-out: AUC with percentile-bootstrap CI, Brier score and its
#     scaled version, calibration intercept/slope, calibration plot
#   * decision curve analysis: net benefit across thresholds vs treat-all
#   * SHAP values for the XGBoost model
#
# Outputs
#   output/model_comparison.csv, output/bootstrap_optimism.csv,
#   output/calibration_data.csv, output/decision_curve_data.csv,
#   output/test_predictions.csv, output/xgb_importance.csv
#   figures/roc_curves.png, calibration.png, decision_curves.png,
#   shap_beeswarm.png, feature_importance.png, bootstrap_optimism.png
# ---------------------------------------------------------------------------

source(local({ d <- getwd(); while (!file.exists(file.path(d, "config.yml")) && dirname(d) != d) d <- dirname(d); file.path(d, "R", "utils.R") }))
suppressPackageStartupMessages({
  library(glmnet); library(ranger); library(xgboost); library(pROC); library(ggplot2)
})
cfg <- load_config()
mod <- cfg$modelling
set.seed(mod$seed)
t0 <- Sys.time()

# Allow `Rscript R/03_models.R --n-boot 200` for the full run.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2 && args[1] == "--n-boot") mod$n_boot <- as.integer(args[2])

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
feat <- fread(file.path(cfg$paths$output_dir, "feature_matrix.csv"))
dict <- fread(file.path(cfg$paths$output_dir, "feature_dictionary.csv"))
predictors <- dict[role == "predictor", feature]
clinical   <- dict[clinical_model == 1, feature]
stopifnot(all(predictors %in% names(feat)), !"outcome" %in% predictors)

# Stratified hold-out split.
idx1 <- which(feat$outcome == 1); idx0 <- which(feat$outcome == 0)
test_idx <- c(sample(idx1, round(mod$test_fraction * length(idx1))),
              sample(idx0, round(mod$test_fraction * length(idx0))))
train <- feat[-test_idx]; test <- feat[test_idx]
msg("Train n = %d (%.1f%% outcome), test n = %d (%.1f%% outcome), %d predictors",
    nrow(train), 100 * mean(train$outcome), nrow(test), 100 * mean(test$outcome), length(predictors))

model_keys   <- c("logistic", "elastic_net", "random_forest", "xgboost")
model_labels <- c(logistic = "Logistic regression (clinical)", elastic_net = "Elastic net",
                  random_forest = "Random forest", xgboost = "XGBoost")
# Fixed colour per model (identity never changes across figures).
model_cols <- c(logistic = "#2a78d6", elastic_net = "#eb6834",
                random_forest = "#1baf7a", xgboost = "#eda100")

as_matrix <- function(dt) as.matrix(dt[, ..predictors])
xgb_params <- list(objective = "binary:logistic", eval_metric = "logloss",
                   eta = 0.05, max_depth = 3, min_child_weight = 5,
                   subsample = 0.8, colsample_bytree = 0.8, lambda = 1)

# ---------------------------------------------------------------------------
# Hyper-parameters tuned once on the training set, then held fixed in the
# bootstrap (a pragmatic shortcut; the fully honest version re-tunes inside
# every resample, which multiplies run time without changing the picture).
# ---------------------------------------------------------------------------
msg("Tuning elastic-net lambda (5-fold CV) and XGBoost rounds (5-fold CV)")
cv_enet <- cv.glmnet(as_matrix(train), train$outcome, family = "binomial", alpha = 0.5, nfolds = 5)
cv_xgb  <- xgb.cv(params = xgb_params, data = xgb.DMatrix(as_matrix(train), label = train$outcome),
                  nrounds = mod$xgb_max_rounds, nfold = 5, early_stopping_rounds = 30, verbose = 0)
best_rounds <- tryCatch(cv_xgb$best_iteration, error = function(e) NULL)
if (is.null(best_rounds)) {
  log <- as.data.frame(cv_xgb$evaluation_log)
  best_rounds <- log$iter[which.min(log$test_logloss_mean)]
}
hp <- list(lambda = cv_enet$lambda.min, nrounds = max(50L, as.integer(best_rounds)))
msg("  lambda.min = %.4f, xgboost rounds = %d", hp$lambda, hp$nrounds)

# ---------------------------------------------------------------------------
# Fit / predict wrappers
# ---------------------------------------------------------------------------
fit_models <- function(tr, hp, importance = FALSE) {
  X <- as_matrix(tr)
  list(
    # suppressWarnings: resamples occasionally produce a cell with fitted
    # probabilities of 0/1 in the glm, which is harmless for prediction.
    logistic = suppressWarnings(glm(reformulate(clinical, "outcome"), data = tr, family = binomial())),
    elastic_net = glmnet(X, tr$outcome, family = "binomial", alpha = 0.5, lambda = hp$lambda),
    # Leaves of at least 40 children and depth <= 8 stop the forest from
    # memorising the training set, which keeps the bootstrap optimism meaningful.
    random_forest = ranger(x = tr[, ..predictors], y = factor(tr$outcome, levels = c(0, 1)),
                           probability = TRUE, num.trees = mod$rf_num_trees,
                           min.node.size = 40, max.depth = 8,
                           importance = if (importance) "impurity" else "none",
                           seed = sample.int(1e6, 1), verbose = FALSE),
    xgboost = xgb.train(params = xgb_params, data = xgb.DMatrix(X, label = tr$outcome),
                        nrounds = hp$nrounds, verbose = 0)
  )
}

predict_models <- function(fits, newdata) {
  X <- as_matrix(newdata)
  data.table(
    logistic = as.numeric(predict(fits$logistic, newdata = newdata, type = "response")),
    elastic_net = as.numeric(predict(fits$elastic_net, newx = X, type = "response")),
    random_forest = predict(fits$random_forest, data = newdata[, ..predictors])$predictions[, "1"],
    xgboost = as.numeric(predict(fits$xgboost, newdata = X))
  )
}

# ---------------------------------------------------------------------------
# Fit on the full training set
# ---------------------------------------------------------------------------
msg("Fitting the four models on the training set")
fits <- fit_models(train, hp, importance = TRUE)
p_train <- predict_models(fits, train)
p_test  <- predict_models(fits, test)

apparent <- rbindlist(lapply(model_keys, function(m) {
  cs <- calibration_stats(train$outcome, p_train[[m]])
  data.table(model = m, auc_apparent = fast_auc(train$outcome, p_train[[m]]), slope_apparent = cs[["slope"]])
}))

# ---------------------------------------------------------------------------
# Bootstrap optimism (Harrell): refit on a resample, measure the gap between
# performance on the resample and on the original training data.
# ---------------------------------------------------------------------------
msg("Bootstrap optimism with %d resamples", mod$n_boot)
opt <- vector("list", mod$n_boot)
for (b in seq_len(mod$n_boot)) {
  boot <- train[sample(.N, replace = TRUE)]
  fb <- fit_models(boot, hp)
  pb <- predict_models(fb, boot); po <- predict_models(fb, train)
  opt[[b]] <- rbindlist(lapply(model_keys, function(m) data.table(
    boot = b, model = m,
    auc_optimism = fast_auc(boot$outcome, pb[[m]]) - fast_auc(train$outcome, po[[m]]),
    slope_optimism = calibration_stats(boot$outcome, pb[[m]])[["slope"]] -
      calibration_stats(train$outcome, po[[m]])[["slope"]]
  )))
  if (b %% 10 == 0) msg("  resample %d/%d", b, mod$n_boot)
}
opt <- rbindlist(opt)
opt_summary <- opt[, .(auc_optimism = mean(auc_optimism), slope_optimism = mean(slope_optimism, na.rm = TRUE)), by = model]

# ---------------------------------------------------------------------------
# Hold-out performance
# ---------------------------------------------------------------------------
msg("Evaluating on the hold-out set")
thr <- seq(mod$dca_thresholds[[1]], mod$dca_thresholds[[2]], by = mod$dca_thresholds[[3]])
boot_ci <- function(y, p, B) {
  n <- length(y)
  a <- replicate(B, { i <- sample.int(n, replace = TRUE); fast_auc(y[i], p[i]) })
  quantile(a, c(0.025, 0.975), names = FALSE)
}
brier_null <- brier(test$outcome, mean(train$outcome))

test_perf <- rbindlist(lapply(model_keys, function(m) {
  p <- p_test[[m]]; y <- test$outcome
  ci <- boot_ci(y, p, mod$n_boot_ci)
  cs <- calibration_stats(y, p)
  nb <- net_benefit(y, p, mod$clinical_threshold)
  data.table(model = m,
             auc_test = fast_auc(y, p), auc_test_ci_low = ci[1], auc_test_ci_high = ci[2],
             brier = brier(y, p), brier_scaled = 1 - brier(y, p) / brier_null,
             calibration_intercept = cs[["intercept"]], calibration_slope = cs[["slope"]],
             net_benefit_at_threshold = nb)
}))

comparison <- merge(merge(apparent, opt_summary, by = "model"), test_perf, by = "model")
comparison[, auc_optimism_corrected := auc_apparent - auc_optimism]
comparison[, slope_optimism_corrected := slope_apparent - slope_optimism]
comparison[, n_predictors := ifelse(model == "logistic", length(clinical), length(predictors))]
comparison[, model_label := model_labels[model]]
comparison <- comparison[match(model_keys, model)]
setcolorder(comparison, c("model", "model_label", "n_predictors", "auc_apparent", "auc_optimism",
                          "auc_optimism_corrected", "auc_test", "auc_test_ci_low", "auc_test_ci_high",
                          "brier", "brier_scaled", "calibration_intercept", "calibration_slope",
                          "slope_optimism_corrected", "net_benefit_at_threshold"))
comparison[, net_benefit_treat_all := net_benefit_all(test$outcome, mod$clinical_threshold)]
print(comparison[, .(model_label, auc_apparent = round(auc_apparent, 3),
                     auc_corrected = round(auc_optimism_corrected, 3), auc_test = round(auc_test, 3),
                     brier = round(brier, 4), cal_slope = round(calibration_slope, 2),
                     nb = round(net_benefit_at_threshold, 4))])

# Decision curves
dca <- rbindlist(c(
  lapply(model_keys, function(m) data.table(model = m, threshold = thr, net_benefit = net_benefit(test$outcome, p_test[[m]], thr))),
  list(data.table(model = "treat_all", threshold = thr, net_benefit = net_benefit_all(test$outcome, thr)),
       data.table(model = "treat_none", threshold = thr, net_benefit = 0))
))

# Calibration by decile of predicted risk
calib <- rbindlist(lapply(model_keys, function(m) {
  p <- p_test[[m]]
  g <- cut(p, breaks = unique(quantile(p, probs = seq(0, 1, 0.1))), include.lowest = TRUE)
  d <- data.table(p = p, y = test$outcome, g = g)[, .(predicted = mean(p), observed = mean(y), n = .N), by = g]
  d[, se := sqrt(observed * (1 - observed) / n)]
  d[, model := m][]
}))

# ---------------------------------------------------------------------------
# Save tables
# ---------------------------------------------------------------------------
out <- cfg$paths$output_dir
fwrite(comparison, file.path(out, "model_comparison.csv"))
fwrite(opt, file.path(out, "bootstrap_optimism.csv"))
fwrite(calib, file.path(out, "calibration_data.csv"))
fwrite(dca, file.path(out, "decision_curve_data.csv"))
fwrite(cbind(test[, .(child_id, outcome)], p_test), file.path(out, "test_predictions.csv"))
xgb_imp <- as.data.table(xgb.importance(model = fits$xgboost))
fwrite(xgb_imp, file.path(out, "xgb_importance.csv"))
saveRDS(list(fits = fits, hp = hp, predictors = predictors, clinical = clinical),
        file.path(out, "models.rds"))

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
fig <- cfg$paths$figures_dir
theme_lab <- function() {
  theme_minimal(base_size = 12) +
    theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
          axis.line = element_line(colour = "#c3c2b7", linewidth = 0.3),
          text = element_text(colour = "#0b0b0b"),
          axis.text = element_text(colour = "#52514e"),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "#52514e"),
          plot.caption = element_text(colour = "#898781", hjust = 0),
          legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(colour = "#0b0b0b", face = "bold"))
}
synthetic_caption <- "Synthetic cohort generated by R/00_simulate_cohort.R. Not study results."
save_fig <- function(p, name, w = 8, h = 5.5) ggsave(file.path(fig, name), p, width = w, height = h, dpi = 150, bg = "#fcfcfb")

# ROC overlay
roc_dt <- rbindlist(lapply(model_keys, function(m) {
  r <- pROC::roc(test$outcome, p_test[[m]], quiet = TRUE)
  data.table(model = m, fpr = 1 - r$specificities, tpr = r$sensitivities)
}))
auc_lab <- comparison[, sprintf("%s  AUC %.2f (%.2f-%.2f)", model_labels[model], auc_test, auc_test_ci_low, auc_test_ci_high)]
names(auc_lab) <- comparison$model
p_roc <- ggplot(roc_dt, aes(fpr, tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_path(linewidth = 0.8) +
  scale_colour_manual(values = model_cols, labels = auc_lab, breaks = model_keys) +
  coord_equal() +
  labs(title = "Discrimination on the hold-out set", subtitle = "ROC curves, synthetic cohort",
       x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)",
       caption = synthetic_caption) +
  guides(colour = guide_legend(ncol = 2)) + theme_lab()
save_fig(p_roc, "roc_curves.png", 7, 7.5)

# Calibration
calib[, model_label := factor(model_labels[model], levels = model_labels)]
cal_dt <- rbindlist(lapply(model_keys, function(m) data.table(model_label = model_labels[m], p = p_test[[m]], y = test$outcome)))
cal_dt[, model_label := factor(model_label, levels = model_labels)]
p_cal <- ggplot(calib, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_smooth(data = cal_dt, aes(p, y), method = "loess", formula = y ~ x, se = FALSE,
              colour = "#898781", linewidth = 0.6, span = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, observed - 1.96 * se), ymax = pmin(1, observed + 1.96 * se)),
                width = 0, colour = "#52514e", linewidth = 0.4) +
  geom_point(aes(colour = model), size = 2.4) +
  scale_colour_manual(values = model_cols, guide = "none") +
  facet_wrap(~ model_label, ncol = 2) +
  coord_equal(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
  labs(title = "Calibration on the hold-out set",
       subtitle = "Points: deciles of predicted risk with 95% CI. Grey line: loess smooth.\nDiagonal: perfect calibration.",
       x = "Predicted probability", y = "Observed proportion", caption = synthetic_caption) +
  theme_lab()
save_fig(p_cal, "calibration.png", 8, 8.5)

# Decision curves
dca[, kind := fifelse(model %in% model_keys, "model", "reference")]
dca_models <- dca[kind == "model"]; dca_ref <- dca[kind == "reference"]
p_dca <- ggplot() +
  geom_line(data = dca_ref, aes(threshold, net_benefit, linetype = model), colour = "#898781", linewidth = 0.6) +
  geom_line(data = dca_models, aes(threshold, net_benefit, colour = model), linewidth = 0.8) +
  scale_colour_manual(values = model_cols, labels = model_labels, breaks = model_keys) +
  scale_linetype_manual(values = c(treat_all = "dashed", treat_none = "solid"),
                        labels = c(treat_all = "Treat all", treat_none = "Treat none")) +
  coord_cartesian(ylim = c(-0.01, max(dca$net_benefit) * 1.05)) +
  labs(title = "Decision curve analysis on the hold-out set",
       subtitle = sprintf("Net benefit of acting on predicted risk above a threshold (prevalence %.1f%%)", 100 * mean(test$outcome)),
       x = "Risk threshold", y = "Net benefit", caption = synthetic_caption) +
  theme_lab() + theme(legend.box = "vertical")
save_fig(p_dca, "decision_curves.png", 8, 6)

# Feature importance: XGBoost gain and random-forest impurity, top 15 each
rf_imp <- data.table(Feature = names(fits$random_forest$variable.importance),
                     value = fits$random_forest$variable.importance)
rf_imp[, value := value / sum(value)][, method := "Random forest (impurity, share)"]
xg <- xgb_imp[, .(Feature, value = Gain / sum(Gain), method = "XGBoost (gain, share)")]
imp <- rbind(head(xg[order(-value)], 15), head(rf_imp[order(-value)], 15))
# Order bars within each facet: encode facet in the level name, strip it from the label.
imp[, key := paste(Feature, method, sep = " | ")]
imp[, key := factor(key, levels = imp[order(value), key])]
p_imp <- ggplot(imp, aes(value, key)) +
  geom_col(fill = "#2a78d6", width = 0.6) +
  scale_y_discrete(labels = function(x) sub(" \\| .*$", "", x)) +
  facet_wrap(~ method, scales = "free") +
  labs(title = "Feature importance (top 15 per model)", x = "Share of total importance", y = NULL,
       caption = synthetic_caption) +
  theme_lab() + theme(panel.grid.major.y = element_blank())
save_fig(p_imp, "feature_importance.png", 10, 6)

# Bootstrap optimism distribution
opt[, model_label := factor(model_labels[model], levels = model_labels)]
p_opt <- ggplot(opt, aes(model_label, auc_optimism, colour = model)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.6) +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 10, colour = "#0b0b0b") +
  scale_colour_manual(values = model_cols, guide = "none") +
  labs(title = sprintf("AUC optimism across %d bootstrap resamples", mod$n_boot),
       subtitle = "AUC on the resample minus AUC of the same fit on the original training data.\nThe black bar is the mean optimism, which is subtracted from the apparent AUC.",
       x = NULL, y = "AUC optimism", caption = synthetic_caption) +
  theme_lab()
save_fig(p_opt, "bootstrap_optimism.png", 8, 5)

# SHAP beeswarm for XGBoost (shapviz; falls back to predcontrib + own plot)
shap_plot <- tryCatch({
  library(shapviz)
  sv <- shapviz(fits$xgboost, X_pred = as_matrix(test))
  p <- sv_importance(sv, kind = "bee", max_display = 20) +
    labs(title = "SHAP values for the XGBoost model (hold-out set)",
         subtitle = "Each point is one child; colour is the feature value (low to high)",
         caption = synthetic_caption)
  p <- suppressMessages(p + scale_colour_gradient(low = "#9ec5f4", high = "#0d366b", breaks = NULL))
  p + theme_lab() + theme(legend.position = "right", panel.grid.major.y = element_blank())
}, error = function(e) {
  msg("shapviz failed (%s); using predcontrib fallback", conditionMessage(e))
  X <- as_matrix(test)
  S <- predict(fits$xgboost, newdata = X, predcontrib = TRUE)
  S <- S[, colnames(S) != "BIAS", drop = FALSE]
  top <- names(sort(colMeans(abs(S)), decreasing = TRUE))[1:20]
  long <- rbindlist(lapply(top, function(f) data.table(feature = f, shap = S[, f], value = scale(X[, f])[, 1])))
  long[, feature := factor(feature, levels = rev(top))]
  ggplot(long, aes(shap, feature, colour = pmin(pmax(value, -2), 2))) +
    geom_jitter(height = 0.2, size = 0.8, alpha = 0.6) +
    scale_colour_gradient(low = "#9ec5f4", high = "#0d366b", name = "Feature value (z)") +
    labs(title = "SHAP values for the XGBoost model (hold-out set)", x = "SHAP value (log-odds)", y = NULL,
         caption = synthetic_caption) + theme_lab() + theme(legend.position = "right")
})
save_fig(shap_plot, "shap_beeswarm.png", 9, 8)

msg("03_models.R finished in %.1f s", as.numeric(Sys.time() - t0, units = "secs"))
