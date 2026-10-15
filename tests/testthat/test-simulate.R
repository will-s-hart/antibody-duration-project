test_that("reviews contain the same measurements and no future visits", {
  cfg <- read_config()
  cell <- mock_design_cell(n_participants = 1L, followup_days = 730)
  full <- simulate_full(cell, cfg, horizon_days = 730)
  shorter <- cell
  shorter$followup_days <- 91
  early <- reduce_to_review(full, shorter)
  late <- reduce_to_review(full, cell)
  expect_silent(validate_sim_dataset(early))
  expect_equal(early$obs, late$obs[late$obs$time_days <= 91, ])
  expect_true(all(early$obs$time_days <= 91))
  expect_equal(early$truth, late$truth)
  expect_false(isTRUE(early$meta$mock))
  expect_equal(full$obs, simulate_full(cell, cfg, horizon_days = 730)$obs)
})

test_that("censoring records missing values rather than assay limits", {
  cfg <- read_config()
  cell <- mock_design_cell(n_participants = 1L, followup_days = 365)
  full <- simulate_full(cell, cfg)
  cell$lloq <- .2
  cell$uloq <- .5
  review <- reduce_to_review(full, cell)
  expect_true(any(review$obs$censor == "left"))
  expect_true(any(review$obs$censor == "right"))
  expect_true(all(is.na(review$obs$y_obs[review$obs$censor != "none"])))
  expect_silent(validate_sim_dataset(review))
})
