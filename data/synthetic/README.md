# data/synthetic

This directory holds only files written by `R/00_simulate_cohort.R`:

| File | Contents |
|---|---|
| `births.csv` | synthetic birth registration records |
| `gp_visits.csv` | synthetic general practice encounters |
| `ed_presentations.csv` | synthetic emergency department presentations |
| `injected_faults.csv` | answer key for the injected linkage faults |
| `linked.sqlite` | the same tables in SQLite, plus the clean and feature tables built by later steps |

They are generated from the seed in `config.yml`, are not committed, and
contain **no real data**. Run `Rscript R/00_simulate_cohort.R` (or
`Rscript run_all.R`) to create them.
