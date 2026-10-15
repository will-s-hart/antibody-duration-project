#!/usr/bin/env Rscript
# Run from the repository root. No private archived files are read at runtime.
suppressMessages(pkgload::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
pilot <- run_prior_pilot(refit = "--refit" %in% args)
cat("\nCompleted ", length(pilot$fits), " fits; ",
    sum(vapply(pilot$fits, function(x) x$diagnostics$converged, logical(1))),
    " passed diagnostics.\n", sep = "")
cat("Results: results/prior-pilot/pilot.rds\n")
