test_that("both prior sets share the early anchors in every family", {
  cfg <- read_config()
  for (set in c("default", "broader")) {
    for (model in CANDIDATE_MODELS) {
      seed <- derive_seed(paste("prior-test", set, model), 1L)
      free <- draw_early_prior(model, cfg, set, 200, seed)
      p <- early_transform(free, model)
      mu <- early_log_curve(p, model, c(0, 30, 90, 365))
      expect_equal(mu[, 1], free$log_c0, tolerance = 1e-10)
      expect_equal(mu[, 3], free$log_c0 - exp(free$log_D),
                   tolerance = 1e-10)
      expect_true(all(mu[, 2] <= mu[, 1]))
      expect_true(all(mu[, 4] <= mu[, 3]))
      for (j in seq_len(10)) {
        ordinary <- log(curve_eval(model, c(0, 30, 90, 365),
                                    as.list(p[j, curve_params(model)])))
        expect_equal(as.numeric(mu[j, ]), ordinary, tolerance = 1e-8)
      }
    }
  }
})

test_that("durations preserve ordering, infinite mass and long finite tails", {
  cfg <- read_config()
  for (model in CANDIDATE_MODELS) {
    seed <- derive_seed(paste("duration-test", model), 1L)
    p <- early_transform(draw_early_prior(model, cfg, "broader", 100,
                                          seed), model)
    higher <- early_tstar(p, model, .202)
    lower <- early_tstar(p, model, .03)
    expect_true(all(higher <= lower))
    finite <- which(is.finite(lower) & lower > 0)
    for (j in utils::head(finite, 10)) {
      expect_equal(as.numeric(exp(early_log_curve(p[j, ], model, lower[j]))),
                   .03, tolerance = 1e-8)
    }
    if (model != "plateau") expect_true(all(is.finite(lower)))
  }
  p <- early_transform(data.frame(log_c0 = log(.92), log_D = log(1e-5),
                                  log_sigma = log(.2)), "exponential")
  expect_gt(early_tstar(p, "exponential", .03), 1e6)
  p <- early_transform(data.frame(log_c0 = log(.92), log_D = log(1.16),
                                  log_sigma = log(.2), u_inf = .5),
                        "plateau")
  expect_identical(early_tstar(p, "plateau", .03), Inf)
  expect_true(is.finite(early_tstar(p, "plateau", .202)))
  expect_identical(early_tstar(p, "plateau", 1), 0)
})

test_that("prior RNG use is reproducible and restores the caller's state", {
  cfg <- read_config()
  seed <- derive_seed("prior-rng-test", 1L)
  a <- draw_early_prior("biphasic", cfg, "broader", 20, seed)
  b <- draw_early_prior("biphasic", cfg, "broader", 20, seed)
  expect_identical(a, b)
  expect_equal(early_prior(cfg, "broader")$slow_fraction_beta, c(1, 1.5))
})
