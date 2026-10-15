test_that("flexible families have identical joint priors at all anchors", {
  cfg <- read_config()
  for (set in c("anchored_default", "anchored_broader")) {
    seed <- derive_seed(paste("two-anchor-test", set), 1L)
    anchors <- list()
    for (model in c("biphasic", "powerlaw", "plateau")) {
      free <- draw_early_prior(model, cfg, set, 500, seed)
      p <- early_transform(free, model)
      mu <- early_log_curve(p, model, c(0, 90, 365))
      expect_equal(mu[, 1], free$log_c0, tolerance = 1e-10)
      expect_equal(mu[, 2], free$log_c0 - exp(free$log_D),
                     tolerance = 1e-10)
      expected <- free$log_c0 - exp(free$log_D) *
        (1 + (365 / 90 - 1) * plogis(free$z_rho))
      expect_equal(mu[, 3], expected, tolerance = 1e-9)
      anchors[[model]] <- mu
      expect_true(all(early_tstar(p, model, .202) <=
                        early_tstar(p, model, .03)))
      for (j in 1:5) {
        native <- curve_eval(model, c(0, 90, 365),
                              as.list(p[j, curve_params(model)]))
        expect_equal(log(native), as.numeric(mu[j, ]), tolerance = 1e-8)
      }
      if (model == "biphasic") {
        expect_true(all(p$k1 > p$k2 & p$k2 > 0))
        expect_true(all(p$f > 0 & p$f < 1))
        expect_true(any(p$f > .5))
      }
    }
    expect_equal(anchors$biphasic, anchors$powerlaw, tolerance = 1e-9)
    expect_equal(anchors$biphasic, anchors$plateau, tolerance = 1e-9)
    a <- draw_early_prior("exponential", cfg, set, 50, seed)
    b <- draw_early_prior("exponential", cfg,
                          sub("anchored_", "", set), 50, seed)
    expect_identical(a, b)
  }
})

test_that("two-anchor mixture recovers the calibrated generating curve", {
  cfg <- read_config()
  b <- cfg$calibration$biphasic
  k1 <- log(2) / b$h1_days
  k2 <- log(2) / b$h2_days
  f <- plogis(-(k1 - k2) * b$ts_days)
  ratio <- function(t) (1-f) * exp(-k1*t) + f * exp(-k2*t)
  D <- -log(ratio(90))
  rho <- (-log(ratio(365)) - D) / ((365/90-1)*D)
  d <- data.frame(log_c0 = log(b$c0), log_D = log(D),
                   log_sigma = log(.2), z_rho = qlogis(rho),
                   u_s = k2*90/(rho*D))
  p <- early_transform(d, "biphasic")
  expect_equal(p$f, f, tolerance = 1e-10)
  expect_equal(p$h1, b$h1_days, tolerance = 1e-8)
  expect_equal(p$h2, b$h2_days, tolerance = 1e-8)
  expect_equal(p$ts, b$ts_days, tolerance = 1e-8)
})

test_that("extreme finite power-law times stay distinct from non-crossing", {
  d <- data.frame(log_c0 = log(.92), log_D = log(1.16),
                   log_sigma = log(.2), z_rho = qlogis(.001))
  p <- early_transform(d, "powerlaw")
  log_time <- early_log_tstar(p, "powerlaw", .03)
  expect_true(is.finite(log_time))
  expect_gt(log_time, log(.Machine$double.xmax))
  expect_error(early_tstar(p, "powerlaw", .03), "numeric range")
  expect_identical(early_log_tstar(p, "powerlaw", 1), -Inf)
})

test_that("JAGS implements the two-anchor transformations consistently", {
  skip_if(!nzchar(Sys.which("jags")))
  skip_if_not_installed("posterior")
  skip_if_not_installed("coda")
  skip_if_not_installed("loo")
  cfg <- read_config()
  sampler <- cfg$fitting$sampler
  sampler$adapt <- 100
  sampler$warmup <- 100
  sampler$draws <- 200
  data <- mock_sim_dataset(mock_design_cell(n_participants = 1L))
  for (model in c("biphasic", "powerlaw", "plateau")) {
    fit <- fit_mcmc(data, model, cfg, prior_set = "anchored_default",
                     sampler = sampler)
    expect_silent(validate_fit_result(fit))
    expect_equal(ncol(fit$posterior),
                   length(early_free_names(model, "anchored_default")) + 1L)
    # Monitor the JAGS means on a reproducible repeat, then independently
    # reconstruct them from the free draws in R. CODA writes six digits.
    cmd <- file.path(fit$work_dir, "run.cmd")
    commands <- readLines(cmd)
    commands <- sub("monitor log_c0", "monitor mu\nmonitor log_c0",
                      commands, fixed = TRUE)
    writeLines(commands, cmd)
    status <- system2(Sys.which("jags"), shQuote(cmd),
                        stdout = FALSE, stderr = FALSE)
    expect_equal(status, 0L)
    stem <- file.path(fit$work_dir, "samples")
    chains <- lapply(seq_len(sampler$chains), function(i) {
      as.matrix(coda::read.coda(paste0(stem, "chain", i, ".txt"),
                                 paste0(stem, "index.txt"), quiet = TRUE))
    })
    samples <- as.data.frame(do.call(rbind, chains))
    p <- early_transform(samples, model)
    expected <- early_log_curve(p, model, data$obs$time_days)
    observed <- as.matrix(samples[paste0("mu[", seq_len(nrow(data$obs)),
                                          "]")])
    expect_equal(unname(observed), unname(expected), tolerance = 2e-4)
  }
})
