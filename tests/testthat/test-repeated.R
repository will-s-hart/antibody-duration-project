test_that("empirical CRPS integrates the CDF without tail cancellation", {
  expect_equal(duration_crps(c(0, 1, 2), 1), 2/9)
  expect_equal(duration_crps(rep(10, 5), 8), 2)
  expect_identical(duration_crps(rep(0, 4), 0), 0)
  seed <- derive_seed("crps-integral-test", 1L)
  x <- .with_seed(seed, exp(rnorm(40)))
  y <- 2
  brute <- mean(abs(x-y)) - mean(abs(outer(x, x, "-")))/2
  expect_equal(duration_crps(x, y), brute, tolerance = 1e-12)
  expect_equal(duration_crps(c(1, 1e200), 2), 2.5e199)
  expect_identical(duration_crps(c(1, Inf), 2), Inf)
  expect_error(duration_crps(c(1, NaN), 2))
})

test_that("bounded and restricted scores retain the non-crossing mass", {
  x <- c(0, 100, Inf)
  z <- duration_score_transform(x)
  expect_equal(z, c(0, 100/465.25, 1))
  expect_true(all(is.finite(z)))
  expect_equal(duration_crps(rep(1, 4), z[2]), 1-z[2])
  expect_equal(duration_crps(pmin(rep(Inf, 4), 3652.5), 100), 3552.5)
})

test_that("coverage bounds account for failures and infinite estimates", {
  cfg <- read_config()
  cfg$repeated$n_replicates <- 2L
  cfg$repeated$models <- "plateau"
  cfg$design$followup_days <- 730
  cell <- mock_design_cell(n_participants = 1L, followup_days = 730)
  data <- mock_sim_dataset(cell)
  fit <- mock_fit_result(cell, model = "plateau", n_draws = 20L,
                           participant_ids = 1L)
  fit$tstar_draws <- fit$tstar_draws[
    fit$tstar_draws$estimand == "individual", ]
  fit$prior_set <- "default"
  fit$diagnostics$converged <- TRUE
  fit$diagnostics$ess_tail_min <- NA_real_
  draws <- c(100, 100, rep(Inf, 18))
  fit$threshold_draws <- list(symptomatic = draws, severe = draws)
  fit$log_lik <- matrix(0, 20, nrow(data$obs))
  truth <- c(symptomatic = 138.7, severe = 1513.9)
  scores <- repeated_score_fit(fit, data, cfg, truth)
  expect_silent(validate_score_rows(scores))
  expect_true(all(scores$covered))
  expect_true(all(is.infinite(scores$bias)))
  expect_true(all(is.finite(scores$crps_bounded)))
  status <- data.frame(model = "plateau", followup_days = 730,
                         completed = c(TRUE, FALSE),
                         converged = c(TRUE, FALSE))
  summary <- summarise_repeated_study(scores, status, cfg)
  expect_equal(summary$n_ok, c(1L, 1L))
  expect_equal(summary$coverage_all_lower, c(.5, .5))
  expect_equal(summary$coverage_all_upper, c(1, 1))
  expect_true(all(is.infinite(summary$bias_days)))
  expect_true(all(is.infinite(summary$median_width_days)))
  expect_true(all(is.finite(summary$crps_restricted_days)))
  expect_true(all(is.na(summary$bias_days_mcse)))
  expect_gt(repeated_wilson(0, 100)[2], 0)
  expect_lt(repeated_wilson(100, 100)[1], 1)
})

test_that("skipping LOO leaves posterior draws and inference unchanged", {
  skip_if(!nzchar(Sys.which("jags")))
  skip_if_not_installed("coda")
  skip_if_not_installed("posterior")
  skip_if_not_installed("loo")
  cfg <- read_config()
  sampler <- cfg$fitting$sampler
  sampler$adapt <- 100
  sampler$warmup <- 100
  sampler$draws <- 200
  data <- mock_sim_dataset(mock_design_cell(n_participants = 1L))
  full <- fit_mcmc(data, "exponential", cfg, sampler = sampler)
  no_loo <- fit_mcmc(data, "exponential", cfg, sampler = sampler,
                      compute_loo = FALSE)
  expect_identical(full$posterior, no_loo$posterior)
  expect_identical(full$threshold_draws, no_loo$threshold_draws)
  expect_true(is.na(no_loo$ic$loo))
  expect_true(is.na(no_loo$loo_reliable))
  expect_silent(validate_fit_result(no_loo))
})
