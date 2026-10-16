#!/usr/bin/env Rscript
# Run from the repository root. No private archived files are read at runtime.
suppressMessages(pkgload::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
cfg <- read_config()
suffix <- if ("--two-anchor" %in% args) "-two-anchor" else ""
if ("--initial-rate" %in% args && nzchar(suffix)) {
  stop("Choose either --initial-rate or --two-anchor")
}
if (nzchar(suffix)) {
  cfg$pilot$prior_sets <- paste0("anchored_", cfg$pilot$prior_sets)
}
if ("--initial-rate" %in% args) {
  suffix <- "-initial-rate"
  cfg$pilot$prior_sets <- "initial_rate"
}
pilot <- run_prior_pilot(
  cfg, cache_dir = paste0("cache/prior-pilot", suffix),
  results_dir = paste0("results/prior-pilot", suffix),
  refit = "--refit" %in% args
)
cat("\nCompleted ", length(pilot$fits), " fits; ",
    sum(vapply(pilot$fits, function(x) x$diagnostics$converged, logical(1))),
    " passed diagnostics.\n", sep = "")
cat("Results: results/prior-pilot", suffix, "/pilot.rds\n", sep = "")
