-- ===========================================================================
-- 02_features.sql
-- SQL-first feature engineering on the cleaned linked tables.
--
-- Builds `features_base`: one row per child with the birth record, the
-- outcome, and health-service utilisation summarised up to the landmark age.
-- Derived transformations (z-scores, cyclic season encoding, imputation,
-- interactions) are added in R/02_features.R.
--
-- Placeholders filled by R from config.yml (defaults in brackets):
--   {landmark_days}     age (days) at which the prediction is made     [60]
--   {horizon_days}      end of the outcome window                       [365]
--   {extract_end}       extract date, for follow-up calculation         [2024-12-31]
--   {outcome_codes_gp}  boolean expression over gp_visits_clean.reason_code
--                       e.g. (g.reason_code LIKE 'J21%' OR g.reason_code = 'R06.2')
--   {outcome_codes_ed}  same expression over ed_presentations_clean.diagnosis_code
--
-- Leakage rule: no utilisation feature may use an encounter after the
-- landmark, and no feature may use an outcome-coded encounter. Children with
-- an outcome code on or before the landmark are flagged (early_outcome = 1)
-- and excluded from the analysis cohort in R.
-- ===========================================================================

-- name: create_features_base
CREATE TABLE features_base AS
WITH gp AS (
  -- every GP encounter with age at visit and code groupings
  SELECT g.child_id,
         julianday(g.visit_date) - julianday(b.dob)              AS age_days,
         g.reason_code,
         CASE WHEN {outcome_codes_gp} THEN 1 ELSE 0 END          AS is_outcome_code,
         CASE WHEN g.reason_code = 'J06.9' THEN 1 ELSE 0 END     AS is_urti,
         CASE WHEN g.reason_code IN ('Z00.1', 'Z23', 'Z76.2') THEN 1 ELSE 0 END AS is_routine
  FROM gp_visits_clean g
  JOIN births_clean b ON b.child_id = g.child_id
),
ed AS (
  SELECT e.child_id,
         julianday(e.presentation_date) - julianday(b.dob)       AS age_days,
         e.diagnosis_code, e.triage_category,
         CASE WHEN {outcome_codes_ed} THEN 1 ELSE 0 END          AS is_outcome_code
  FROM ed_presentations_clean e
  JOIN births_clean b ON b.child_id = e.child_id
),
first_outcome AS (
  -- first wheeze/bronchiolitis-coded encounter in either setting
  SELECT child_id, MIN(age_days) AS first_outcome_day
  FROM (
    SELECT child_id, age_days FROM gp WHERE is_outcome_code = 1
    UNION ALL
    SELECT child_id, age_days FROM ed WHERE is_outcome_code = 1
  )
  GROUP BY child_id
),
gp_lm AS (
  -- GP utilisation up to the landmark, excluding outcome-coded encounters
  SELECT child_id,
         COUNT(*)                                                 AS gp_visits_lm,
         SUM(CASE WHEN age_days <= 30 THEN 1 ELSE 0 END)          AS gp_visits_0_30,
         SUM(is_routine)                                          AS gp_routine_lm,
         SUM(1 - is_routine)                                      AS gp_problem_lm,
         SUM(is_urti)                                             AS gp_urti_lm,
         COUNT(DISTINCT reason_code)                              AS gp_distinct_codes_lm,
         MIN(age_days)                                            AS days_to_first_gp
  FROM gp
  WHERE age_days <= {landmark_days} AND is_outcome_code = 0
  GROUP BY child_id
),
ed_lm AS (
  -- ED utilisation up to the landmark, excluding outcome-coded encounters
  SELECT child_id,
         COUNT(*)                                                 AS ed_visits_lm,
         SUM(CASE WHEN triage_category <= 3 THEN 1 ELSE 0 END)    AS ed_high_acuity_lm
  FROM ed
  WHERE age_days <= {landmark_days} AND is_outcome_code = 0
  GROUP BY child_id
),
sibs AS (
  -- other births to the same mother within this extract. Pairs born 1-200
  -- days apart fail the plausibility check in 01_linkage_checks.sql (the
  -- mother_id link is wrong), so they are not counted; same-day births (twins)
  -- and intervals over 200 days are.
  SELECT a.child_id, COUNT(b.child_id) AS siblings_in_extract
  FROM births_clean a
  LEFT JOIN births_clean b
    ON b.mother_id = a.mother_id AND b.child_id <> a.child_id
   AND (ABS(julianday(a.dob) - julianday(b.dob)) = 0
        OR ABS(julianday(a.dob) - julianday(b.dob)) > 200)
  GROUP BY a.child_id
)
SELECT b.child_id, b.mother_id, b.dob, b.sex, b.gestational_age_weeks, b.birth_weight_g,
       b.delivery_mode, b.parity, b.maternal_age, b.maternal_asthma,
       b.maternal_smoking_pregnancy, b.seifa_quintile, b.season_of_birth,
       COALESCE(s.siblings_in_extract, 0)         AS siblings_in_extract,
       COALESCE(gp_lm.gp_visits_lm, 0)            AS gp_visits_lm,
       COALESCE(gp_lm.gp_visits_0_30, 0)          AS gp_visits_0_30,
       COALESCE(gp_lm.gp_routine_lm, 0)           AS gp_routine_lm,
       COALESCE(gp_lm.gp_problem_lm, 0)           AS gp_problem_lm,
       COALESCE(gp_lm.gp_urti_lm, 0)              AS gp_urti_lm,
       COALESCE(gp_lm.gp_distinct_codes_lm, 0)    AS gp_distinct_codes_lm,
       gp_lm.days_to_first_gp,
       COALESCE(ed_lm.ed_visits_lm, 0)            AS ed_visits_lm,
       COALESCE(ed_lm.ed_high_acuity_lm, 0)       AS ed_high_acuity_lm,
       fo.first_outcome_day,
       CASE WHEN fo.first_outcome_day <= {landmark_days} THEN 1 ELSE 0 END AS early_outcome,
       CASE WHEN fo.first_outcome_day >  {landmark_days}
             AND fo.first_outcome_day <= {horizon_days} THEN 1 ELSE 0 END  AS outcome,
       julianday('{extract_end}') - julianday(b.dob) AS follow_up_days
FROM births_clean b
LEFT JOIN sibs          s     ON s.child_id = b.child_id
LEFT JOIN gp_lm               ON gp_lm.child_id = b.child_id
LEFT JOIN ed_lm               ON ed_lm.child_id = b.child_id
LEFT JOIN first_outcome fo    ON fo.child_id = b.child_id;

-- name: outcome_summary
-- Quick check of the cohort flow after the landmark exclusion.
SELECT COUNT(*)                                                    AS children,
       SUM(early_outcome)                                          AS excluded_early_outcome,
       SUM(CASE WHEN follow_up_days < {horizon_days} THEN 1 ELSE 0 END) AS excluded_short_follow_up,
       SUM(CASE WHEN early_outcome = 0 AND follow_up_days >= {horizon_days} THEN outcome ELSE 0 END) AS outcomes,
       ROUND(100.0 * SUM(CASE WHEN early_outcome = 0 AND follow_up_days >= {horizon_days} THEN outcome ELSE 0 END)
             / SUM(CASE WHEN early_outcome = 0 AND follow_up_days >= {horizon_days} THEN 1 ELSE 0 END), 1) AS prevalence_pct
FROM features_base;
