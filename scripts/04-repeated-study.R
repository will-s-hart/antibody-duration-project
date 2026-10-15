#!/usr/bin/env Rscript
# Resumable fixed-truth repeated fits; all posterior scoring is predeclared.
suppressMessages(pkgload::load_all(quiet = TRUE))
cfg <- read_config()
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  x <- args[startsWith(args, paste0("--", name, "="))]
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
cfg$repeated$n_replicates <- as.integer(option("replicates",
                                             cfg$repeated$n_replicates))
workers <- as.integer(option("workers", cfg$repeated$workers))
if (workers < 1L || cfg$repeated$n_replicates < 1L) {
  stop("Positive replicate and worker counts are required")
}
out <- option("results", "results/repeated-default")
cache <- option("cache", "cache/repeated-default")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
yaml::write_yaml(cfg, file.path(out, "configuration.yaml"))
if (!"--collect-only" %in% args) {
  paths <- parallel::mclapply(seq_len(cfg$repeated$n_replicates),
    function(r) run_repeated_replicate(r, cfg, cache, out),
    mc.cores = workers, mc.preschedule = FALSE, mc.set.seed = FALSE)
  if (any(vapply(paths, inherits, logical(1), "try-error"))) {
    stop("A replicate worker failed; completed fit caches are retained")
  }
}
study <- collect_repeated_study(cfg, out)
cat("Completed ", nrow(study$status), " planned fits; ",
      sum(study$status$converged), " passed convergence.\n", sep = "")
