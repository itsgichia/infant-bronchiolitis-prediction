# Shared helpers used by every pipeline script.
#
# Scripts are run from anywhere; find_root() walks up from the working
# directory until it finds config.yml, so paths in config.yml are always
# resolved relative to the project root.

suppressPackageStartupMessages({
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(yaml)
})

#' Locate the project root (the directory containing config.yml).
find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  for (i in 1:10) {
    if (file.exists(file.path(d, "config.yml"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not find config.yml above ", start, ". Run scripts from the project root.")
}

#' Read config.yml and resolve path entries against the project root.
load_config <- function(root = find_root()) {
  cfg <- yaml::read_yaml(file.path(root, "config.yml"))
  cfg$root <- root
  cfg$paths <- lapply(cfg$paths, function(p) file.path(root, p))
  for (p in c("synthetic_dir", "output_dir", "figures_dir")) {
    dir.create(cfg$paths[[p]], showWarnings = FALSE, recursive = TRUE)
  }
  cfg
}

#' Open (or create) the SQLite database holding the linked tables.
open_db <- function(cfg) DBI::dbConnect(RSQLite::SQLite(), cfg$paths$sqlite)

#' Timestamped progress message.
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", sprintf(...), "\n")

#' Read a file of named SQL statements.
#'
#' Statements are separated by header lines of the form `-- name: some_name`.
#' Returns a named character vector so the same SQL can be run from R and
#' read on its own in the sql/ folder.
read_sql_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("^--\\s*name:\\s*", lines)
  if (!length(starts)) stop("No '-- name:' headers in ", path)
  ends <- c(starts[-1] - 1L, length(lines))
  out <- character(0)
  for (i in seq_along(starts)) {
    nm <- trimws(sub("^--\\s*name:\\s*", "", lines[starts[i]]))
    body <- lines[(starts[i] + 1L):ends[i]]
    # Drop section banners that follow the statement's closing semicolon.
    last_sql <- max(grep(";", body), 0L)
    if (last_sql > 0L) body <- body[seq_len(last_sql)]
    out[nm] <- paste(body, collapse = "\n")
  }
  out
}

#' Substitute `{placeholder}` tokens in a SQL string.
fill_sql <- function(sql, params) {
  for (nm in names(params)) {
    sql <- gsub(paste0("\\{", nm, "\\}"), as.character(params[[nm]]), sql)
  }
  sql
}

#' Fast rank-based AUC (equivalent to the Mann-Whitney statistic).
fast_auc <- function(y, p) {
  y <- as.integer(y)
  n1 <- sum(y == 1L); n0 <- length(y) - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

logit <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  log(p / (1 - p))
}

#' Brier score (mean squared error of the probability forecast).
brier <- function(y, p) mean((as.numeric(y) - p)^2)

#' Calibration intercept (with slope fixed at 1) and slope.
calibration_stats <- function(y, p) {
  lp <- logit(p)
  y <- as.integer(y)
  slope <- tryCatch(coef(glm(y ~ lp, family = binomial()))[["lp"]], error = function(e) NA_real_)
  intercept <- tryCatch(coef(glm(y ~ 1 + offset(lp), family = binomial()))[[1]], error = function(e) NA_real_)
  c(intercept = intercept, slope = slope)
}

#' Net benefit of acting on `p >= threshold`, per Vickers & Elkin (2006).
#'
#' NB = TP/N - FP/N * t/(1 - t). "Treat all" is the same expression with
#' every person treated; "treat none" is zero by definition.
net_benefit <- function(y, p, thresholds) {
  y <- as.integer(y); n <- length(y)
  vapply(thresholds, function(t) {
    treat <- p >= t
    tp <- sum(treat & y == 1L); fp <- sum(treat & y == 0L)
    tp / n - fp / n * t / (1 - t)
  }, numeric(1))
}

net_benefit_all <- function(y, thresholds) {
  prev <- mean(as.integer(y))
  prev - (1 - prev) * thresholds / (1 - thresholds)
}
