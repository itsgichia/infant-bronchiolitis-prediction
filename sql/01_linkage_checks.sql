-- ===========================================================================
-- 01_linkage_checks.sql
-- Linkage validation and plausibility checks on the raw linked extract.
--
-- Each statement is headed by `-- name: <id>` so R/01_linkage_checks.R can
-- run them one at a time and log the number of rows each one returns.
-- Detection statements return the offending records; cleaning statements
-- build *_clean tables and never modify the raw extract.
--
-- Placeholders filled by R from config.yml:
--   {extract_end}   ISO date of the data extract (default 2024-12-31)
--
-- SQLite dialect (julianday() for date arithmetic, window functions).
-- ===========================================================================


-- ---------------------------------------------------------------------------
-- A. Birth registration record checks
-- ---------------------------------------------------------------------------

-- name: duplicate_birth_record
-- A child_id should appear exactly once. n_versions > 1 means the duplicates
-- disagree on at least one clinical field (conflicting duplicate).
SELECT child_id,
       COUNT(*) AS n_rows,
       COUNT(DISTINCT COALESCE(gestational_age_weeks, -1) || '|' ||
                      COALESCE(birth_weight_g, -1)        || '|' ||
                      delivery_mode) AS n_versions
FROM births
GROUP BY child_id
HAVING COUNT(*) > 1;

-- name: implausible_birth_weight
-- Viable live births lie roughly between 300 g and 6500 g. Values below 10
-- are almost always kilograms keyed into a grams field (recoverable);
-- values above 6500 are usually an extra digit.
SELECT child_id, birth_weight_g,
       CASE WHEN birth_weight_g BETWEEN 0.3 AND 6.5 THEN 'kg_entered_as_g'
            WHEN birth_weight_g < 300 THEN 'too_low'
            ELSE 'too_high' END AS sub_class
FROM births
WHERE birth_weight_g < 300 OR birth_weight_g > 6500;

-- name: missing_gestational_age
SELECT child_id
FROM births
WHERE gestational_age_weeks IS NULL;

-- name: implausible_gestational_age
-- Registered live births are 20-44 completed weeks.
SELECT child_id, gestational_age_weeks
FROM births
WHERE gestational_age_weeks < 20 OR gestational_age_weeks > 44;

-- name: implausible_maternal_age
SELECT child_id, maternal_age
FROM births
WHERE maternal_age < 12 OR maternal_age > 55;

-- name: invalid_sex_code
SELECT child_id, sex
FROM births
WHERE sex NOT IN ('M', 'F');

-- name: birth_weight_gestation_mismatch
-- Plausibility cross-check between two fields: a term baby under 1000 g or
-- an extremely preterm baby over 2500 g points at a keying or linkage error.
SELECT child_id, gestational_age_weeks, birth_weight_g
FROM births
WHERE (gestational_age_weeks >= 37 AND birth_weight_g BETWEEN 300 AND 1000)
   OR (gestational_age_weeks <= 28 AND birth_weight_g BETWEEN 2500 AND 6500);

-- name: implausible_sibling_interval
-- Two different children registered to the same mother less than 200 days
-- apart cannot both be singleton births; the mother_id link is unreliable.
SELECT a.mother_id, a.child_id AS child_a, b.child_id AS child_b,
       ABS(julianday(a.dob) - julianday(b.dob)) AS days_apart
FROM births a
JOIN births b
  ON a.mother_id = b.mother_id AND a.child_id < b.child_id
WHERE ABS(julianday(a.dob) - julianday(b.dob)) BETWEEN 1 AND 200;


-- ---------------------------------------------------------------------------
-- B. Encounter record checks (general practice and emergency department)
-- ---------------------------------------------------------------------------

-- name: gp_orphan_visit
-- Encounters whose child_id has no birth record: the linkage key failed.
SELECT g.gp_record_id, g.child_id
FROM gp_visits g
LEFT JOIN births b ON b.child_id = g.child_id
WHERE b.child_id IS NULL;

-- name: ed_orphan_visit
SELECT e.ed_record_id, e.child_id
FROM ed_presentations e
LEFT JOIN births b ON b.child_id = e.child_id
WHERE b.child_id IS NULL;

-- name: gp_visit_before_birth
-- An encounter before the date of birth is a mis-link (or a wrong DOB).
SELECT g.gp_record_id, g.child_id, g.visit_date, b.dob
FROM gp_visits g
JOIN (SELECT child_id, MIN(dob) AS dob FROM births GROUP BY child_id) b
  ON b.child_id = g.child_id
WHERE julianday(g.visit_date) < julianday(b.dob);

-- name: ed_visit_before_birth
SELECT e.ed_record_id, e.child_id, e.presentation_date, b.dob
FROM ed_presentations e
JOIN (SELECT child_id, MIN(dob) AS dob FROM births GROUP BY child_id) b
  ON b.child_id = e.child_id
WHERE julianday(e.presentation_date) < julianday(b.dob);

-- name: gp_visit_after_extract
-- Dates after the extract was cut cannot be real.
SELECT gp_record_id, child_id, visit_date
FROM gp_visits
WHERE julianday(visit_date) > julianday('{extract_end}');

-- name: ed_visit_after_extract
SELECT ed_record_id, child_id, presentation_date
FROM ed_presentations
WHERE julianday(presentation_date) > julianday('{extract_end}');

-- name: gp_duplicate_visit_row
-- The same encounter supplied more than once in the extract.
SELECT gp_record_id, COUNT(*) AS n_rows
FROM gp_visits
GROUP BY gp_record_id
HAVING COUNT(*) > 1;

-- name: ed_duplicate_visit_row
SELECT ed_record_id, COUNT(*) AS n_rows
FROM ed_presentations
GROUP BY ed_record_id
HAVING COUNT(*) > 1;

-- name: ed_invalid_triage
-- Australasian Triage Scale categories are 1-5.
SELECT ed_record_id, triage_category
FROM ed_presentations
WHERE triage_category IS NULL OR triage_category NOT BETWEEN 1 AND 5;

-- name: gp_missing_reason_code
SELECT gp_record_id
FROM gp_visits
WHERE reason_code IS NULL OR TRIM(reason_code) = '';


-- ---------------------------------------------------------------------------
-- C. Cleaning: build analysis-ready tables (raw tables are left untouched)
-- ---------------------------------------------------------------------------

-- name: clean_births
-- One row per child. Among duplicates keep the most complete record, then
-- the first supplied. Implausible values are set to NULL except kilogram
-- birth weights, which are rescaled to grams.
CREATE TABLE births_clean AS
WITH ranked AS (
  SELECT b.*,
         ROW_NUMBER() OVER (
           PARTITION BY child_id
           ORDER BY (gestational_age_weeks IS NOT NULL) +
                    (birth_weight_g IS NOT NULL) +
                    (maternal_age IS NOT NULL) DESC,
                    b.rowid
         ) AS rn
  FROM births b
)
SELECT child_id, mother_id, dob, sex,
       CASE WHEN gestational_age_weeks < 20 OR gestational_age_weeks > 44
            THEN NULL ELSE gestational_age_weeks END          AS gestational_age_weeks,
       CASE WHEN birth_weight_g BETWEEN 0.3 AND 6.5 THEN ROUND(birth_weight_g * 1000)
            WHEN birth_weight_g < 300 OR birth_weight_g > 6500 THEN NULL
            ELSE birth_weight_g END                           AS birth_weight_g,
       delivery_mode, parity,
       CASE WHEN maternal_age < 12 OR maternal_age > 55
            THEN NULL ELSE maternal_age END                   AS maternal_age,
       maternal_asthma, maternal_smoking_pregnancy, seifa_quintile, season_of_birth
FROM ranked
WHERE rn = 1;

-- name: clean_gp_visits
-- Drop duplicate rows, orphans, pre-birth and post-extract encounters.
CREATE TABLE gp_visits_clean AS
SELECT DISTINCT g.gp_record_id, g.child_id, g.visit_date, g.reason_code
FROM gp_visits g
JOIN births_clean b ON b.child_id = g.child_id
WHERE julianday(g.visit_date) >= julianday(b.dob)
  AND julianday(g.visit_date) <= julianday('{extract_end}');

-- name: clean_ed_presentations
CREATE TABLE ed_presentations_clean AS
SELECT DISTINCT e.ed_record_id, e.child_id, e.presentation_date,
       e.diagnosis_code, e.triage_category
FROM ed_presentations e
JOIN births_clean b ON b.child_id = e.child_id
WHERE julianday(e.presentation_date) >= julianday(b.dob)
  AND julianday(e.presentation_date) <= julianday('{extract_end}')
  AND e.triage_category BETWEEN 1 AND 5;
