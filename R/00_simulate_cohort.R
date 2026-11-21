# 00_simulate_cohort.R
# ---------------------------------------------------------------------------
# Generate a SYNTHETIC linked cohort that mimics the *structure* of the data
# used in the real study (birth registration, general practice and emergency
# department records for infants). No real records are used or reproduced.
#
# Outputs
#   data/synthetic/births.csv            one row per registered birth
#   data/synthetic/gp_visits.csv         general practice encounters
#   data/synthetic/ed_presentations.csv  emergency department presentations
#   data/synthetic/injected_faults.csv   answer key for the injected faults
#   data/synthetic/linked.sqlite         the same three tables in SQLite
#
# Design of the generator (assumptions are listed in README.md):
#   * Outcome = any J21.x (bronchiolitis) or R06.2 (wheeze) diagnosis in GP or
#     ED records in the first 365 days of life, target prevalence ~18%.
#   * True risk factors with modest log-odds: maternal asthma, maternal
#     smoking in pregnancy, prematurity, low birth weight, autumn/winter birth,
#     older siblings, male sex, socio-economic disadvantage.
#   * Respiratory events cluster in the Southern-hemisphere winter and in the
#     first six months of life.
#   * Health-service use is higher among cases and among families with more
#     children, so early utilisation carries some (weak) signal.
#   * Linkage faults are injected on top of the clean data and recorded in an
#     answer key so the checks in 01_linkage_checks.R can be tested.
# ---------------------------------------------------------------------------

# Locate R/utils.R by walking up from the working directory to the project root.
source(local({ d <- getwd(); while (!file.exists(file.path(d, "config.yml")) && dirname(d) != d) d <- dirname(d); file.path(d, "R", "utils.R") }))
cfg <- load_config()
sim <- cfg$simulation
set.seed(sim$seed)
t0 <- Sys.time()

n <- sim$n_births
extract_end <- as.Date(sim$extract_end)

# ---------------------------------------------------------------------------
# 1. Birth registration records
# ---------------------------------------------------------------------------
msg("Simulating %d births", n)

child_id  <- sprintf("C%06d", seq_len(n))
dob <- as.Date("2022-01-01") + sample(0:729, n, replace = TRUE)

# ~3% of mothers have a second child in the window (a sibling, 10-23 months
# later) or a twin (same day). Intervals are kept physiologically plausible
# here; implausible intervals are injected later as a linkage fault.
mother_id <- sprintf("M%06d", seq_len(n))
first_born <- sample(seq_len(n), size = round(0.03 * n))
second_born <- sample(setdiff(seq_len(n), first_born), length(first_born))
mother_id[second_born] <- mother_id[first_born]
gap <- ifelse(runif(length(first_born)) < 0.03, 0L, sample(300:700, length(first_born), TRUE))
dob2 <- dob[first_born] + gap
flip <- dob2 > as.Date("2023-12-31")               # reflect so the sibling stays in the window
dob2[flip] <- dob[first_born][flip] - gap[flip]
dob[second_born] <- dob2
month <- as.integer(format(dob, "%m"))
season <- c("Summer", "Summer", "Autumn", "Autumn", "Autumn", "Winter",
            "Winter", "Winter", "Spring", "Spring", "Spring", "Summer")[month]

sex <- sample(c("M", "F"), n, replace = TRUE, prob = c(0.512, 0.488))
seifa_quintile <- sample(1:5, n, replace = TRUE, prob = c(0.25, 0.22, 0.20, 0.18, 0.15))
seifa_low <- as.integer(seifa_quintile <= 2)

maternal_age <- pmin(pmax(round(rnorm(n, 30.3, 5.4)), 16), 48)
parity <- sample(0:4, n, replace = TRUE, prob = c(0.42, 0.33, 0.15, 0.07, 0.03))
older_sib <- as.integer(parity > 0)

# Smoking is more common with disadvantage and younger age; asthma ~15%.
p_smoke <- plogis(-2.6 + 0.9 * seifa_low - 0.05 * (maternal_age - 30))
maternal_smoking <- rbinom(n, 1, p_smoke)
maternal_asthma  <- rbinom(n, 1, 0.15)

# Gestational age: 8% preterm with a tail down to 24 weeks.
is_preterm <- rbinom(n, 1, plogis(-2.5 + 0.35 * maternal_smoking + 0.25 * maternal_asthma + 0.2 * seifa_low))
ga <- ifelse(is_preterm == 1,
             pmax(24, round(36 - rexp(n, rate = 0.45))),
             pmin(42, pmax(37, round(rnorm(n, 39.3, 1.1)))))

# Birth weight conditional on gestation, sex, smoking and parity.
bw_mean <- 3450 + 165 * (ga - 40) - 9 * (ga - 40)^2 + 115 * (sex == "M") -
  190 * maternal_smoking + 40 * older_sib
birth_weight <- pmax(450, round(rnorm(n, bw_mean, 400)))

p_cs <- plogis(-0.75 + 0.6 * is_preterm + 0.03 * (maternal_age - 30) - 0.25 * older_sib)
delivery_mode <- ifelse(rbinom(n, 1, p_cs) == 1, "Caesarean", "Vaginal")

births <- data.table(
  child_id, mother_id, dob, sex,
  gestational_age_weeks = ga,
  birth_weight_g = birth_weight,
  delivery_mode, parity, maternal_age,
  maternal_asthma, maternal_smoking_pregnancy = maternal_smoking,
  seifa_quintile, season_of_birth = season
)

# ---------------------------------------------------------------------------
# 2. Latent outcome: wheeze/bronchiolitis diagnosis in the first year
# ---------------------------------------------------------------------------
lbw <- as.integer(birth_weight < 2500)
# Log-odds contributions. These are moderate, literature-plausible effects
# (odds ratios roughly 1.3-2.5) chosen so the synthetic signal is detectable
# but far from deterministic; see README "Simulator assumptions".
lp <- 0.80 * maternal_asthma + 0.60 * maternal_smoking +
  0.70 * is_preterm + 0.12 * pmax(0, 37 - ga) +          # dose-response below term
  0.45 * lbw + 0.70 * (season == "Winter") + 0.40 * (season == "Autumn") +
  0.50 * older_sib + 0.25 * (sex == "M") + 0.30 * seifa_low +
  0.15 * (delivery_mode == "Caesarean") - 0.20 * (parity >= 3) +
  0.30 * maternal_asthma * maternal_smoking +
  0.35 * (season == "Winter") * is_preterm                 # interaction for the tree models
# Choose the intercept so the marginal prevalence hits the target exactly.
b0 <- uniroot(function(b) mean(plogis(b + lp)) - sim$target_prevalence, c(-10, 5))$root
y <- rbinom(n, 1, plogis(b0 + lp))
msg("Outcome prevalence (truth, before faults): %.1f%%", 100 * mean(y))

# ---------------------------------------------------------------------------
# 3. Encounter records
# ---------------------------------------------------------------------------
# Non-outcome reason codes (ICD-10-AM style) and sampling weights.
gp_codes <- c("Z00.1" = 20, "Z23" = 18, "J06.9" = 12, "R50.9" = 7, "L20.9" = 6,
              "R63.3" = 6, "K59.0" = 4, "H66.9" = 5, "B34.9" = 5, "R11" = 4,
              "A09" = 4, "L22" = 4, "Z76.2" = 5)
ed_codes <- c("R50.9" = 22, "A09" = 14, "J06.9" = 14, "S09.9" = 10, "R11" = 8,
              "H66.9" = 7, "T14.9" = 6, "Z03.8" = 8, "R56.0" = 3, "L03.9" = 3,
              "N39.0" = 5)
resp_gp  <- c("J21.9" = 45, "R06.2" = 40, "J21.0" = 15)
resp_ed  <- c("J21.0" = 35, "J21.9" = 45, "R06.2" = 20)
triage_p <- c("1" = 0.01, "2" = 0.08, "3" = 0.38, "4" = 0.43, "5" = 0.10)
triage_resp_p <- c("1" = 0.02, "2" = 0.18, "3" = 0.48, "4" = 0.29, "5" = 0.03)

is_winter_month <- function(d) as.integer(format(d, "%m")) %in% 6:8
is_shoulder_month <- function(d) as.integer(format(d, "%m")) %in% c(5, 9)

#' Sample event ages (days) weighted toward winter and the 2-6 month age band.
sample_event_days <- function(dob_i, k, lo = 14, hi = 365) {
  days <- lo:hi
  dates <- dob_i + days
  w <- (1 + 2.5 * is_winter_month(dates) + 1.0 * is_shoulder_month(dates)) *
    (0.25 + dnorm(days, 135, 80) / dnorm(135, 135, 80))
  sample(days, k, replace = FALSE, prob = w)
}

# 3a. Routine / non-outcome GP visits over two years of follow-up (front-loaded).
mu_gp <- exp(log(4.2) + 0.30 * y + 0.15 * older_sib + 0.30 * is_preterm +
             0.10 * seifa_low - 0.05 * (seifa_quintile == 5) + rnorm(n, 0, 0.30))
n_gp <- rpois(n, mu_gp)
gp_rows <- data.table(
  child_id   = rep(child_id, n_gp),
  age_days   = floor(730 * rbeta(sum(n_gp), 1, 1.7)),
  reason_code = sample(names(gp_codes), sum(n_gp), replace = TRUE, prob = gp_codes)
)
# Early upper-respiratory infections (J06.9, not an outcome code) in weeks
# 2-8 are more common among children who go on to wheeze.
early_urti <- rbinom(n, 1, ifelse(y == 1, 0.32, 0.11))
gp_rows <- rbind(gp_rows, data.table(
  child_id = child_id[early_urti == 1],
  age_days = sample(14:56, sum(early_urti), TRUE),
  reason_code = "J06.9"))
dob_by_child <- setNames(dob, child_id)
gp_rows[, visit_date := dob_by_child[child_id] + age_days]

# 3b. Non-outcome ED presentations.
mu_ed <- 0.38 + 0.15 * y + 0.10 * seifa_low + 0.25 * is_preterm
n_ed <- rpois(n, mu_ed)
ed_rows <- data.table(
  child_id = rep(child_id, n_ed),
  age_days = floor(730 * rbeta(sum(n_ed), 1, 1.4)),
  diagnosis_code = sample(names(ed_codes), sum(n_ed), replace = TRUE, prob = ed_codes),
  triage_category = as.integer(sample(names(triage_p), sum(n_ed), replace = TRUE, prob = triage_p))
)
ed_rows[, presentation_date := rep(dob, n_ed) + age_days]

# 3c. Outcome events for cases: 1 + Poisson(0.5) respiratory encounters in year 1.
case_idx <- which(y == 1)
k_events <- 1 + rpois(length(case_idx), 0.5)
ev <- rbindlist(lapply(seq_along(case_idx), function(j) {
  i <- case_idx[j]
  d <- sample_event_days(dob[i], k_events[j])
  data.table(child_id = child_id[i], age_days = d, event_date = dob[i] + d,
             setting = ifelse(runif(k_events[j]) < 0.30, "ED", "GP"))
}))
ev_gp <- ev[setting == "GP", .(child_id, age_days, visit_date = event_date,
                               reason_code = sample(names(resp_gp), .N, replace = TRUE, prob = resp_gp))]
ev_ed <- ev[setting == "ED", .(child_id, age_days, presentation_date = event_date,
                               diagnosis_code = sample(names(resp_ed), .N, replace = TRUE, prob = resp_ed),
                               triage_category = as.integer(sample(names(triage_resp_p), .N, replace = TRUE, prob = triage_resp_p)))]

# 3d. Some controls have wheeze/bronchiolitis codes AFTER the first year, so the
#     outcome window genuinely matters.
late_idx <- sample(which(y == 0), round(0.06 * sum(y == 0)))
late_gp <- data.table(child_id = child_id[late_idx], age_days = sample(366:730, length(late_idx), TRUE))
late_gp[, visit_date := dob[late_idx] + age_days]
late_gp[, reason_code := sample(names(resp_gp), .N, replace = TRUE, prob = resp_gp)]

gp_visits <- rbindlist(list(gp_rows, ev_gp, late_gp), use.names = TRUE)[, .(child_id, visit_date, reason_code)]
ed_presentations <- rbindlist(list(ed_rows, ev_ed), use.names = TRUE)[, .(child_id, presentation_date, diagnosis_code, triage_category)]

# Censor at the extract date and give every encounter a record id.
gp_visits <- gp_visits[visit_date <= extract_end]
ed_presentations <- ed_presentations[presentation_date <= extract_end]
setorder(gp_visits, child_id, visit_date)
setorder(ed_presentations, child_id, presentation_date)
gp_visits[, gp_record_id := sprintf("GP%07d", .I)]
ed_presentations[, ed_record_id := sprintf("ED%06d", .I)]
setcolorder(gp_visits, "gp_record_id")
setcolorder(ed_presentations, "ed_record_id")

msg("Clean records: %d GP visits, %d ED presentations", nrow(gp_visits), nrow(ed_presentations))

# ---------------------------------------------------------------------------
# 4. Inject linkage and data-quality faults (and keep an answer key)
# ---------------------------------------------------------------------------
faults <- list()
dob_lookup <- setNames(births$dob, births$child_id)   # clean births, before duplicates
log_fault <- function(fault, table, ids) {
  faults[[length(faults) + 1]] <<- data.table(fault = fault, source_table = table, record_id = ids)
}

# F1 Duplicate child records: 20 exact copies, 10 conflicting copies.
dup_exact <- sample(child_id, 20)
dup_conf  <- sample(setdiff(child_id, dup_exact), 10)
dup_rows_exact <- births[child_id %in% dup_exact]
dup_rows_conf  <- births[child_id %in% dup_conf]
dup_rows_conf[, birth_weight_g := birth_weight_g + sample(c(-60, -40, 40, 60), .N, TRUE)]
dup_rows_conf[, delivery_mode := ifelse(delivery_mode == "Vaginal", "Caesarean", "Vaginal")]
births <- rbindlist(list(births, dup_rows_exact, dup_rows_conf))
log_fault("duplicate_birth_record", "births", c(dup_exact, dup_conf))

# F2 Encounters dated before the child's birth (classic mis-link).
pre_gp <- sample(gp_visits$gp_record_id, 45)
gp_visits[gp_record_id %in% pre_gp,
          visit_date := dob_lookup[child_id] - sample(1:400, .N, TRUE)]
pre_ed <- sample(ed_presentations$ed_record_id, 12)
ed_presentations[ed_record_id %in% pre_ed,
                 presentation_date := dob_lookup[child_id] - sample(1:400, .N, TRUE)]
log_fault("visit_before_birth", "gp_visits", pre_gp)
log_fault("visit_before_birth", "ed_presentations", pre_ed)

# F3 Implausible birth weights: kilograms entered as grams, impossible lows, extra digits.
bw_ids <- sample(setdiff(child_id, c(dup_exact, dup_conf)), 18)
births[child_id %in% bw_ids[1:6],   birth_weight_g := round(birth_weight_g / 1000, 2)]
births[child_id %in% bw_ids[7:12],  birth_weight_g := sample(c(90, 120, 150, 210, 260, 280), .N)]
births[child_id %in% bw_ids[13:18], birth_weight_g := birth_weight_g * 10]
log_fault("implausible_birth_weight", "births", bw_ids)

# F4 Missing and implausible gestational age.
ga_na <- sample(setdiff(child_id, bw_ids), round(0.035 * n))
births[child_id %in% ga_na, gestational_age_weeks := NA_integer_]
ga_bad <- sample(setdiff(child_id, c(ga_na, bw_ids)), 8)
births[child_id %in% ga_bad, gestational_age_weeks := sample(c(12L, 18L, 47L, 52L), .N, TRUE)]
log_fault("missing_gestational_age", "births", ga_na)
log_fault("implausible_gestational_age", "births", ga_bad)

# F5 Orphan encounters: child ids with no birth record (failed linkage).
orphan_ids <- sprintf("C9%05d", sample(1:99999, 40))
orphan_gp <- data.table(
  gp_record_id = sprintf("GP%07d", nrow(gp_visits) + 1:70),
  child_id = sample(orphan_ids, 70, TRUE),
  visit_date = as.Date("2022-03-01") + sample(0:900, 70, TRUE),
  reason_code = sample(names(gp_codes), 70, TRUE, prob = gp_codes))
orphan_ed <- data.table(
  ed_record_id = sprintf("ED%06d", nrow(ed_presentations) + 1:20),
  child_id = sample(orphan_ids, 20, TRUE),
  presentation_date = as.Date("2022-03-01") + sample(0:900, 20, TRUE),
  diagnosis_code = sample(names(ed_codes), 20, TRUE, prob = ed_codes),
  triage_category = as.integer(sample(names(triage_p), 20, TRUE, prob = triage_p)))
gp_visits <- rbind(gp_visits, orphan_gp)
ed_presentations <- rbind(ed_presentations, orphan_ed)
log_fault("orphan_visit", "gp_visits", orphan_gp$gp_record_id)
log_fault("orphan_visit", "ed_presentations", orphan_ed$ed_record_id)

# F6 Exact duplicate encounter rows (same record id repeated in the extract).
dup_gp <- sample(gp_visits[!gp_record_id %in% orphan_gp$gp_record_id, gp_record_id], 60)
gp_visits <- rbind(gp_visits, gp_visits[gp_record_id %in% dup_gp])
log_fault("duplicate_visit_row", "gp_visits", dup_gp)

# F7 Encounters dated after the extract date.
fut_gp <- sample(gp_visits[!gp_record_id %in% c(pre_gp, dup_gp, orphan_gp$gp_record_id), gp_record_id], 10)
gp_visits[gp_record_id %in% fut_gp, visit_date := extract_end + sample(30:400, .N, TRUE)]
log_fault("visit_after_extract", "gp_visits", fut_gp)

# F8 Implausible maternal age.
ma_ids <- sample(setdiff(child_id, c(bw_ids, ga_bad, dup_exact, dup_conf)), 6)
births[child_id %in% ma_ids, maternal_age := sample(c(9L, 11L, 64L, 71L), .N, TRUE)]
log_fault("implausible_maternal_age", "births", ma_ids)

# F9 Mother-id mis-links: two unrelated children (born 1-200 days apart) given
#    the same mother_id. Recorded against the mother_id, as the check reports it.
singleton_m <- names(which(table(births$mother_id) == 1))
pair_a <- sample(births[mother_id %in% singleton_m, child_id], 12)
pair_b <- vapply(pair_a, function(a) {
  d <- dob_lookup[[a]]
  cands <- births[mother_id %in% singleton_m & !child_id %in% pair_a &
                  abs(as.integer(dob - d)) %between% c(1, 200), child_id]
  sample(cands, 1)
}, character(1))
mother_a <- births$mother_id[match(pair_a, births$child_id)]
births[child_id %in% pair_b, mother_id := mother_a[match(child_id, pair_b)]]
log_fault("implausible_sibling_interval", "births", mother_a)

# Shuffle row order so faults are not trivially positional.
births <- births[sample(.N)]
gp_visits <- gp_visits[sample(.N)]
ed_presentations <- ed_presentations[sample(.N)]
fault_key <- rbindlist(faults)

# ---------------------------------------------------------------------------
# 5. Write CSVs and the SQLite database
# ---------------------------------------------------------------------------
out <- cfg$paths$synthetic_dir
fwrite(births, file.path(out, "births.csv"))
fwrite(gp_visits, file.path(out, "gp_visits.csv"))
fwrite(ed_presentations, file.path(out, "ed_presentations.csv"))
fwrite(fault_key, file.path(out, "injected_faults.csv"))

if (file.exists(cfg$paths$sqlite)) invisible(file.remove(cfg$paths$sqlite))
con <- open_db(cfg)
# SQLite has no DATE type; store ISO-8601 text so julianday() works in SQL.
to_sql <- function(dt) {
  dt <- copy(dt)
  for (col in names(dt)) if (inherits(dt[[col]], "Date")) dt[[col]] <- format(dt[[col]], "%Y-%m-%d")
  dt
}
dbWriteTable(con, "births", to_sql(births))
dbWriteTable(con, "gp_visits", to_sql(gp_visits))
dbWriteTable(con, "ed_presentations", to_sql(ed_presentations))
dbWriteTable(con, "injected_faults", fault_key)
invisible(dbExecute(con, "CREATE INDEX idx_births_child ON births(child_id)"))
invisible(dbExecute(con, "CREATE INDEX idx_gp_child ON gp_visits(child_id)"))
invisible(dbExecute(con, "CREATE INDEX idx_ed_child ON ed_presentations(child_id)"))
dbDisconnect(con)

msg("Wrote %d birth rows, %d GP rows, %d ED rows, %d injected faults",
    nrow(births), nrow(gp_visits), nrow(ed_presentations), nrow(fault_key))
msg("00_simulate_cohort.R finished in %.1f s", as.numeric(Sys.time() - t0, units = "secs"))
