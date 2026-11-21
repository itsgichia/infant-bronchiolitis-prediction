# 02_features.R
# ---------------------------------------------------------------------------
# Feature engineering.
#
# Step 1 (SQL, sql/02_features.sql): joins and window aggregations build
#         `features_base` in SQLite: one row per child with the birth record,
#         utilisation up to the landmark age, and the outcome.
# Step 2 (R): derived clinical features, missingness handling, cyclic season
#         encoding, interactions, and the landmark exclusion.
#
# Outputs
#   output/feature_matrix.csv       analysis cohort (one row per child)
#   output/feature_dictionary.csv   name, group, description, type, role
#   output/cohort_flow.csv          counts at each exclusion step
#   SQLite tables features_base and features
# ---------------------------------------------------------------------------

source(local({ d <- getwd(); while (!file.exists(file.path(d, "config.yml")) && dirname(d) != d) d <- dirname(d); file.path(d, "R", "utils.R") }))
cfg <- load_config()
des <- cfg$design
t0 <- Sys.time()

con <- open_db(cfg)
sql <- read_sql_file(file.path(cfg$paths$sql_dir, "02_features.sql"))

# Outcome-code expression built from config so SQL and R agree on the definition.
code_expr <- function(col, codes) {
  parts <- vapply(codes, function(cd) {
    if (grepl("\\.", cd)) sprintf("%s = '%s'", col, cd) else sprintf("%s LIKE '%s%%'", col, cd)
  }, character(1))
  paste0("(", paste(parts, collapse = " OR "), ")")
}
params <- list(
  landmark_days = des$landmark_days,
  horizon_days = des$horizon_days,
  extract_end = cfg$simulation$extract_end,
  outcome_codes_gp = code_expr("g.reason_code", des$outcome_codes),
  outcome_codes_ed = code_expr("e.diagnosis_code", des$outcome_codes)
)

# ---------------------------------------------------------------------------
# Step 1: SQL aggregation
# ---------------------------------------------------------------------------
msg("Building features_base in SQL (landmark %d days, horizon %d days)", des$landmark_days, des$horizon_days)
invisible(dbExecute(con, "DROP TABLE IF EXISTS features_base"))
invisible(dbExecute(con, fill_sql(sql[["create_features_base"]], params)))
print(dbGetQuery(con, fill_sql(sql[["outcome_summary"]], params)))
fb <- as.data.table(dbGetQuery(con, "SELECT * FROM features_base"))
fb[, dob := as.Date(dob)]

# ---------------------------------------------------------------------------
# Step 2: cohort flow and landmark exclusion
# ---------------------------------------------------------------------------
flow <- data.table(step = "Children with a cleaned birth record", n = nrow(fb))
fb <- fb[follow_up_days >= des$horizon_days]
flow <- rbind(flow, data.table(step = sprintf("Follow-up of at least %d days", des$horizon_days), n = nrow(fb)))
fb <- fb[early_outcome == 0]
flow <- rbind(flow, data.table(step = sprintf("No wheeze/bronchiolitis code by day %d (landmark)", des$landmark_days), n = nrow(fb)))
flow <- rbind(flow, data.table(step = "  of which outcome in window", n = sum(fb$outcome)))
print(flow)

# ---------------------------------------------------------------------------
# Step 3: derived features
# ---------------------------------------------------------------------------
lm_days <- des$landmark_days
X <- fb[, .(child_id, outcome)]

# Demographics and perinatal
X[, sex_male := as.integer(fb$sex == "M")]
X[, ga_missing := as.integer(is.na(fb$gestational_age_weeks))]
X[, gestational_age_weeks := fifelse(is.na(fb$gestational_age_weeks),
                                     median(fb$gestational_age_weeks, na.rm = TRUE),
                                     as.numeric(fb$gestational_age_weeks))]
X[, preterm := as.integer(gestational_age_weeks < 37)]
X[, very_preterm := as.integer(gestational_age_weeks < 32)]
X[, bw_missing := as.integer(is.na(fb$birth_weight_g))]
X[, birth_weight_g := fifelse(is.na(fb$birth_weight_g), median(fb$birth_weight_g, na.rm = TRUE), as.numeric(fb$birth_weight_g))]
X[, low_birth_weight := as.integer(birth_weight_g < 2500)]
# Birth-weight z-score within gestational week and sex (cohort-internal reference)
X[, birth_weight_z := {
  z <- (birth_weight_g - mean(birth_weight_g)) / sd(birth_weight_g)
  fifelse(is.na(z), 0, z)
}, by = .(ga_bin = pmin(pmax(round(gestational_age_weeks), 30), 42), sex_male)]
X[, small_for_gestational_age := as.integer(birth_weight_z < -1.28)]   # bottom decile
X[, caesarean := as.integer(fb$delivery_mode == "Caesarean")]

# Maternal and family
X[, parity := fb$parity]
X[, has_older_siblings := as.integer(fb$parity > 0)]
X[, siblings_in_extract := fb$siblings_in_extract]
X[, maternal_age := fifelse(is.na(fb$maternal_age), median(fb$maternal_age, na.rm = TRUE), as.numeric(fb$maternal_age))]
X[, maternal_age_lt20 := as.integer(maternal_age < 20)]
X[, maternal_age_ge35 := as.integer(maternal_age >= 35)]
X[, maternal_asthma := fb$maternal_asthma]
X[, maternal_smoking_pregnancy := fb$maternal_smoking_pregnancy]
X[, seifa_quintile := fb$seifa_quintile]
X[, seifa_low := as.integer(fb$seifa_quintile <= 2)]

# Seasonality (Southern hemisphere: winter = Jun-Aug)
month <- as.integer(format(fb$dob, "%m"))
X[, born_summer := as.integer(fb$season_of_birth == "Summer")]
X[, born_autumn := as.integer(fb$season_of_birth == "Autumn")]
X[, born_winter := as.integer(fb$season_of_birth == "Winter")]
X[, season_sin := sin(2 * pi * (month - 1) / 12)]
X[, season_cos := cos(2 * pi * (month - 1) / 12)]
# Days from birth to the start of the next winter (0 if born in winter):
# the younger a baby is when RSV season starts, the higher the exposure risk.
yr <- as.integer(format(fb$dob, "%Y"))
next_winter <- as.Date(sprintf("%d-06-01", yr))
next_winter[fb$dob > as.Date(sprintf("%d-08-31", yr))] <- as.Date(sprintf("%d-06-01", yr[fb$dob > as.Date(sprintf("%d-08-31", yr))] + 1L))
X[, days_to_first_winter := pmax(0, as.integer(next_winter - fb$dob))]

# Health-service utilisation up to the landmark (outcome codes excluded in SQL)
X[, gp_visits_lm := fb$gp_visits_lm]
X[, gp_visits_0_30 := fb$gp_visits_0_30]
X[, gp_routine_lm := fb$gp_routine_lm]
X[, gp_problem_lm := fb$gp_problem_lm]
X[, gp_urti_lm := fb$gp_urti_lm]
X[, gp_distinct_codes_lm := fb$gp_distinct_codes_lm]
X[, any_gp_lm := as.integer(fb$gp_visits_lm > 0)]
X[, days_to_first_gp := fifelse(is.na(fb$days_to_first_gp), as.numeric(lm_days), fb$days_to_first_gp)]  # censored at landmark
X[, ed_visits_lm := fb$ed_visits_lm]
X[, ed_high_acuity_lm := fb$ed_high_acuity_lm]
X[, any_ed_lm := as.integer(fb$ed_visits_lm > 0)]
X[, total_contacts_lm := fb$gp_visits_lm + fb$ed_visits_lm]

# Interactions and a simple risk count (useful for the clinical baseline model)
X[, asthma_x_smoking := maternal_asthma * maternal_smoking_pregnancy]
X[, preterm_x_winter := preterm * born_winter]
X[, risk_factor_count := maternal_asthma + maternal_smoking_pregnancy + preterm +
      low_birth_weight + born_winter + has_older_siblings]

# ---------------------------------------------------------------------------
# Feature dictionary
# ---------------------------------------------------------------------------
# clinical_model = 1 marks the pre-specified predictors of the baseline
# logistic regression (the ML models use every predictor).
dict <- rbindlist(list(
  list("child_id", "id", "Pseudonymised child identifier", "id", "id", 0),
  list("outcome", "outcome", sprintf("Wheeze (R06.2) or bronchiolitis (J21.x) diagnosis in GP or ED on days %d-%d", lm_days + 1, des$horizon_days), "binary", "outcome", 0),
  list("sex_male", "demographic", "Male sex", "binary", "predictor", 1),
  list("gestational_age_weeks", "perinatal", "Completed weeks of gestation (median-imputed if missing)", "numeric", "predictor", 0),
  list("ga_missing", "perinatal", "Gestational age missing or implausible in source", "binary", "predictor", 0),
  list("preterm", "perinatal", "Born before 37 weeks", "binary", "predictor", 1),
  list("very_preterm", "perinatal", "Born before 32 weeks", "binary", "predictor", 0),
  list("birth_weight_g", "perinatal", "Birth weight in grams (median-imputed if missing)", "numeric", "predictor", 0),
  list("bw_missing", "perinatal", "Birth weight missing or implausible in source", "binary", "predictor", 0),
  list("low_birth_weight", "perinatal", "Birth weight below 2500 g", "binary", "predictor", 1),
  list("birth_weight_z", "perinatal", "Birth-weight z-score within gestational week and sex", "numeric", "predictor", 0),
  list("small_for_gestational_age", "perinatal", "Birth-weight z-score below -1.28 (bottom decile)", "binary", "predictor", 0),
  list("caesarean", "perinatal", "Caesarean delivery", "binary", "predictor", 1),
  list("parity", "family", "Number of previous births to the mother", "count", "predictor", 0),
  list("has_older_siblings", "family", "Parity greater than zero", "binary", "predictor", 1),
  list("siblings_in_extract", "family", "Other children of the same mother registered in the extract", "count", "predictor", 0),
  list("maternal_age", "maternal", "Maternal age at birth in years (median-imputed if missing)", "numeric", "predictor", 1),
  list("maternal_age_lt20", "maternal", "Maternal age under 20", "binary", "predictor", 0),
  list("maternal_age_ge35", "maternal", "Maternal age 35 or over", "binary", "predictor", 0),
  list("maternal_asthma", "maternal", "Mother has a diagnosis of asthma", "binary", "predictor", 1),
  list("maternal_smoking_pregnancy", "maternal", "Mother smoked during pregnancy", "binary", "predictor", 1),
  list("seifa_quintile", "socioeconomic", "SEIFA disadvantage quintile of residence (1 = most disadvantaged)", "ordinal", "predictor", 0),
  list("seifa_low", "socioeconomic", "SEIFA quintile 1 or 2", "binary", "predictor", 1),
  list("born_summer", "seasonal", "Born December-February", "binary", "predictor", 0),
  list("born_autumn", "seasonal", "Born March-May", "binary", "predictor", 1),
  list("born_winter", "seasonal", "Born June-August", "binary", "predictor", 1),
  list("season_sin", "seasonal", "Sine encoding of birth month", "numeric", "predictor", 0),
  list("season_cos", "seasonal", "Cosine encoding of birth month", "numeric", "predictor", 0),
  list("days_to_first_winter", "seasonal", "Days from birth to the next 1 June (0 if born in winter)", "numeric", "predictor", 0),
  list("gp_visits_lm", "utilisation", sprintf("GP visits in days 0-%d, excluding outcome codes", lm_days), "count", "predictor", 1),
  list("gp_visits_0_30", "utilisation", "GP visits in days 0-30, excluding outcome codes", "count", "predictor", 0),
  list("gp_routine_lm", "utilisation", "Routine GP visits (child health check, immunisation) up to landmark", "count", "predictor", 0),
  list("gp_problem_lm", "utilisation", "Problem-oriented GP visits up to landmark", "count", "predictor", 0),
  list("gp_urti_lm", "utilisation", "GP visits coded as upper respiratory tract infection (J06.9) up to landmark", "count", "predictor", 0),
  list("gp_distinct_codes_lm", "utilisation", "Distinct GP reason codes up to landmark", "count", "predictor", 0),
  list("any_gp_lm", "utilisation", "Any GP visit up to landmark", "binary", "predictor", 0),
  list("days_to_first_gp", "utilisation", sprintf("Age in days at first GP visit (censored at %d)", lm_days), "numeric", "predictor", 0),
  list("ed_visits_lm", "utilisation", "ED presentations up to landmark, excluding outcome codes", "count", "predictor", 1),
  list("ed_high_acuity_lm", "utilisation", "ED presentations with triage category 1-3 up to landmark", "count", "predictor", 0),
  list("any_ed_lm", "utilisation", "Any ED presentation up to landmark", "binary", "predictor", 0),
  list("total_contacts_lm", "utilisation", "GP visits plus ED presentations up to landmark", "count", "predictor", 0),
  list("asthma_x_smoking", "interaction", "Maternal asthma and smoking in pregnancy", "binary", "predictor", 0),
  list("preterm_x_winter", "interaction", "Preterm birth in winter", "binary", "predictor", 0),
  list("risk_factor_count", "composite", "Count of six established risk factors", "count", "predictor", 0)
))
setnames(dict, c("feature", "group", "description", "type", "role", "clinical_model"))

# Dictionary and matrix must agree exactly.
stopifnot(setequal(dict$feature, names(X)))
setcolorder(X, dict$feature)
stopifnot(!anyNA(X))

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
fwrite(X, file.path(cfg$paths$output_dir, "feature_matrix.csv"))
fwrite(dict, file.path(cfg$paths$output_dir, "feature_dictionary.csv"))
fwrite(flow, file.path(cfg$paths$output_dir, "cohort_flow.csv"))
dbWriteTable(con, "features", X, overwrite = TRUE)
dbDisconnect(con)

n_pred <- sum(dict$role == "predictor")
msg("Feature matrix: %d children, %d predictors (%d in the clinical baseline), prevalence %.1f%%",
    nrow(X), n_pred, sum(dict$clinical_model == 1), 100 * mean(X$outcome))
msg("02_features.R finished in %.1f s", as.numeric(Sys.time() - t0, units = "secs"))
