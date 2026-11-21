# 01_linkage_checks.R
# ---------------------------------------------------------------------------
# Linkage validation and plausibility checks.
#
# Runs every detection query in sql/01_linkage_checks.sql against the raw
# linked tables, logs how many records each rule catches and what was done
# about them, then builds the *_clean tables used downstream.
#
# Outputs
#   output/data_quality_log.csv      one row per rule (issue, table, n, rule, action)
#   output/data_quality_detail.csv   one row per flagged record (for auditing/tests)
#   SQLite tables births_clean, gp_visits_clean, ed_presentations_clean
#
# In the real study the equivalent log traced 200+ issues back to the
# supplying datasets; this script shows the mechanism on synthetic data.
# ---------------------------------------------------------------------------

source(local({ d <- getwd(); while (!file.exists(file.path(d, "config.yml")) && dirname(d) != d) d <- dirname(d); file.path(d, "R", "utils.R") }))
cfg <- load_config()
t0 <- Sys.time()

con <- open_db(cfg)
sql <- read_sql_file(file.path(cfg$paths$sql_dir, "01_linkage_checks.sql"))
params <- list(extract_end = cfg$simulation$extract_end)

# ---------------------------------------------------------------------------
# Rule book: plain-English rule and the action taken for each detection query.
# `record_col` names the column that identifies the flagged record.
# ---------------------------------------------------------------------------
rules <- rbindlist(list(
  list("duplicate_birth_record", "births", "child_id",
       "child_id appears more than once in the birth extract",
       "Keep the most complete record (then first supplied); conflicting versions listed for manual review"),
  list("implausible_birth_weight", "births", "child_id",
       "birth_weight_g outside 300-6500 g",
       "Values 0.3-6.5 treated as kilograms and rescaled to grams; other values set to NULL"),
  list("missing_gestational_age", "births", "child_id",
       "gestational_age_weeks is NULL",
       "Retain; median-imputed in feature step with a missingness indicator"),
  list("implausible_gestational_age", "births", "child_id",
       "gestational_age_weeks outside 20-44 weeks",
       "Set to NULL (then handled as missing)"),
  list("implausible_maternal_age", "births", "child_id",
       "maternal_age outside 12-55 years",
       "Set to NULL (then handled as missing)"),
  list("invalid_sex_code", "births", "child_id",
       "sex not in {M, F}",
       "Set to NULL"),
  list("birth_weight_gestation_mismatch", "births", "child_id",
       "term birth < 1000 g or <= 28 weeks and > 2500 g",
       "Flag for review; retained unchanged"),
  list("implausible_sibling_interval", "births", "mother_id",
       "two children of the same mother born 1-200 days apart",
       "Retain children; the flagged mother_id links are ignored when counting siblings in the extract"),
  list("gp_orphan_visit", "gp_visits", "gp_record_id",
       "child_id has no birth record",
       "Drop encounter (cannot be linked); counts reported to data custodian"),
  list("ed_orphan_visit", "ed_presentations", "ed_record_id",
       "child_id has no birth record",
       "Drop encounter (cannot be linked); counts reported to data custodian"),
  list("gp_visit_before_birth", "gp_visits", "gp_record_id",
       "visit_date earlier than child's date of birth",
       "Drop encounter (mis-link or DOB error)"),
  list("ed_visit_before_birth", "ed_presentations", "ed_record_id",
       "presentation_date earlier than child's date of birth",
       "Drop encounter (mis-link or DOB error)"),
  list("gp_visit_after_extract", "gp_visits", "gp_record_id",
       "visit_date after the extract date",
       "Drop encounter"),
  list("ed_visit_after_extract", "ed_presentations", "ed_record_id",
       "presentation_date after the extract date",
       "Drop encounter"),
  list("gp_duplicate_visit_row", "gp_visits", "gp_record_id",
       "same gp_record_id supplied more than once",
       "Keep one copy (SELECT DISTINCT)"),
  list("ed_duplicate_visit_row", "ed_presentations", "ed_record_id",
       "same ed_record_id supplied more than once",
       "Keep one copy (SELECT DISTINCT)"),
  list("ed_invalid_triage", "ed_presentations", "ed_record_id",
       "triage_category outside 1-5",
       "Drop encounter"),
  list("gp_missing_reason_code", "gp_visits", "gp_record_id",
       "reason_code missing",
       "Retain; counted as a visit with unknown reason")
))
setnames(rules, c("issue", "source_table", "record_col", "rule", "action"))

# ---------------------------------------------------------------------------
# Run the detection queries
# ---------------------------------------------------------------------------
msg("Running %d detection rules", nrow(rules))
log_rows <- list(); detail_rows <- list()
for (i in seq_len(nrow(rules))) {
  r <- rules[i]
  if (!r$issue %in% names(sql)) stop("No SQL statement named ", r$issue)
  hits <- as.data.table(dbGetQuery(con, fill_sql(sql[[r$issue]], params)))
  n <- nrow(hits)
  log_rows[[i]] <- data.table(issue = r$issue, source_table = r$source_table, n_rows = n,
                              rule = r$rule, action = r$action)
  if (n > 0) {
    detail_rows[[i]] <- data.table(issue = r$issue, source_table = r$source_table,
                                   record_id = as.character(hits[[r$record_col]]))
  }
  msg("  %-34s %-18s %5d rows", r$issue, r$source_table, n)
}
dq_log <- rbindlist(log_rows)
dq_detail <- rbindlist(detail_rows)

# Conflicting duplicates deserve their own line in the log.
dups <- as.data.table(dbGetQuery(con, sql[["duplicate_birth_record"]]))
dq_log <- rbind(dq_log, data.table(
  issue = "duplicate_birth_record_conflicting", source_table = "births",
  n_rows = sum(dups$n_versions > 1),
  rule = "duplicated child_id whose copies disagree on a clinical field",
  action = "Listed in output/data_quality_detail.csv for manual adjudication"))
if (any(dups$n_versions > 1)) {
  dq_detail <- rbind(dq_detail, data.table(issue = "duplicate_birth_record_conflicting",
                                           source_table = "births",
                                           record_id = dups[n_versions > 1, child_id]))
}

# ---------------------------------------------------------------------------
# Build the clean tables
# ---------------------------------------------------------------------------
msg("Building clean tables")
for (tbl in c("births_clean", "gp_visits_clean", "ed_presentations_clean")) {
  invisible(dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", tbl)))
}
invisible(dbExecute(con, fill_sql(sql[["clean_births"]], params)))
invisible(dbExecute(con, fill_sql(sql[["clean_gp_visits"]], params)))
invisible(dbExecute(con, fill_sql(sql[["clean_ed_presentations"]], params)))
for (tbl in c("births_clean", "gp_visits_clean", "ed_presentations_clean")) {
  invisible(dbExecute(con, sprintf("CREATE INDEX IF NOT EXISTS idx_%s ON %s(child_id)", tbl, tbl)))
}

count <- function(tbl) dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", tbl))$n
flow <- data.table(
  table = c("births", "gp_visits", "ed_presentations"),
  raw_rows = c(count("births"), count("gp_visits"), count("ed_presentations")),
  clean_rows = c(count("births_clean"), count("gp_visits_clean"), count("ed_presentations_clean"))
)
flow[, rows_removed := raw_rows - clean_rows]
print(flow)

# Post-condition checks: the clean tables must satisfy the rules.
stopifnot(
  count("births_clean") == dbGetQuery(con, "SELECT COUNT(DISTINCT child_id) AS n FROM births")$n,
  dbGetQuery(con, "SELECT COUNT(*) AS n FROM gp_visits_clean g LEFT JOIN births_clean b USING (child_id) WHERE b.child_id IS NULL")$n == 0,
  dbGetQuery(con, "SELECT COUNT(*) AS n FROM gp_visits_clean g JOIN births_clean b USING (child_id) WHERE julianday(g.visit_date) < julianday(b.dob)")$n == 0,
  dbGetQuery(con, "SELECT COUNT(*) AS n FROM births_clean WHERE birth_weight_g < 300 OR birth_weight_g > 6500")$n == 0
)

fwrite(dq_log, file.path(cfg$paths$output_dir, "data_quality_log.csv"))
fwrite(dq_detail, file.path(cfg$paths$output_dir, "data_quality_detail.csv"))
fwrite(flow, file.path(cfg$paths$output_dir, "cleaning_row_counts.csv"))
dbDisconnect(con)

msg("Logged %d issues across %d rules (%d rules fired)", sum(dq_log$n_rows), nrow(dq_log), sum(dq_log$n_rows > 0))
msg("01_linkage_checks.R finished in %.1f s", as.numeric(Sys.time() - t0, units = "secs"))
