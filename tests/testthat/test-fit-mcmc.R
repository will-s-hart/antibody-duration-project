test_that("implemented engines validate inputs and reject pooled fits", {
  expect_error(fit_mcmc(NULL, "biphasic", NULL), "sim_dataset")
  cfg <- read_config()
  data <- mock_sim_dataset(mock_design_cell(n_participants = 2L))
  expect_error(fit_mcmc(data, cfg = cfg), "one participant")
})

test_that("atomic durations have meaningful convergence diagnostics", {
  skip_if_not_installed("coda")
  skip_if_not_installed("posterior")
  cfg <- read_config()
  seed <- derive_seed("atomic-diagnostic-test", 1L)
  x <- .with_seed(seed, array(rnorm(2000 * 4), c(2000, 4)))
  chains <- do.call(coda::mcmc.list, lapply(1:4, function(i) {
    coda::mcmc(matrix(x[, i], ncol = 1,
                       dimnames = list(NULL, "free_parameter")))
  }))
  d <- rep(Inf, length(x))
  d[seq(1, length(x), by = 10)] <- 100
  dg <- .early_diagnostics(chains, list(severe = d), cfg$fitting$sampler)
  expect_type(dg$converged, "logical")
  expect_false(is.na(dg$converged))
  expect_true(all(is.finite(dg$table$rhat)))
  constant <- .early_diagnostics(chains, list(severe = rep(Inf, length(x))),
                                 cfg$fitting$sampler)
  expect_equal(constant$constant_indicators, "severe")
  expect_true(is.na(constant$table$ess_bulk[
    constant$table$variable == "duration_severe"
  ]))
})

test_that("real JAGS fits emit valid draws without truth leakage", {
  skip_if(!nzchar(Sys.which("jags")), "JAGS command line unavailable")
  skip_if_not_installed("coda")
  skip_if_not_installed("posterior")
  skip_if_not_installed("loo")
  cfg <- read_config()
  sampler <- cfg$fitting$sampler
  sampler$adapt <- 100
  sampler$warmup <- 100
  sampler$draws <- 200
  data <- mock_sim_dataset(mock_design_cell(n_participants = 1L))
  fit <- fit_mcmc(data, "exponential", cfg, sampler = sampler)
  expect_silent(validate_fit_result(fit))
  expect_equal(nrow(fit$tstar_draws), sampler$draws * sampler$chains)
  expect_identical(fit$tstar_draws$crossed,
                   is.finite(fit$tstar_draws$tstar_days))
  expect_true(all(fit$threshold_draws$symptomatic <=
                    fit$threshold_draws$severe))
  data$truth$c0 <- 10000
  data$truth$tstar_days <- 1
  data$meta$population_tstar_days <- 1
  repeated <- fit_mcmc(data, "exponential", cfg, sampler = sampler)
  expect_equal(fit$posterior, repeated$posterior)
})

test_that("initial-rate JAGS curves agree with the independent R transform", {
  skip_if(!nzchar(Sys.which("jags")), "JAGS command line unavailable")
  skip_if_not_installed("coda")
  skip_if_not_installed("posterior")
  cfg <- read_config()
  sampler <- cfg$fitting$sampler
  sampler$adapt <- 100
  sampler$warmup <- 100
  sampler$draws <- 200
  data <- mock_sim_dataset(mock_design_cell(n_participants = 1L))
  for (model in CANDIDATE_MODELS) {
    fit <- fit_mcmc(data, model, cfg, prior_set = "initial_rate",
                     sampler = sampler, compute_loo = FALSE)
    expect_silent(validate_fit_result(fit))
    # Monitoring deterministic nodes leaves the sampler's RNG stream intact.
    command_file <- file.path(fit$work_dir, "run.cmd")
    commands <- readLines(command_file)
    at <- which(commands == paste("update", sampler$draws))
    commands <- append(commands, "monitor mu", after = at - 1L)
    writeLines(commands, command_file)
    log_file <- file.path(fit$work_dir, "check.log")
    status <- system2(Sys.which("jags"), shQuote(command_file),
                       stdout = log_file, stderr = log_file)
    expect_equal(status, 0L)
    chains <- lapply(seq_len(sampler$chains), function(i) {
      as.matrix(coda::read.coda(
        file.path(fit$work_dir, paste0("sampleschain", i, ".txt")),
        file.path(fit$work_dir, "samplesindex.txt"), quiet = TRUE
      ))
    })
    draws <- do.call(rbind, chains)
    free <- early_free_names(model, "initial_rate")
    expect_equal(unname(draws[, free]),
                   unname(as.matrix(fit$posterior[free])))
    mu <- early_log_curve(fit$native_draws, model, data$obs$time_days)
    expect_equal(unname(draws[, paste0("mu[", seq_len(ncol(mu)), "]")]),
                   unname(mu), tolerance = 1e-4)
  }
})
