# A one-person illustration of handout Section 6, steps 1--3.
# It deliberately does not estimate repeated-simulation coverage or adequacy.

.pilot_atomic_save <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path))
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not save ", path)
}

.pilot_endpoint_fit <- function(fit, endpoint, cfg, cell) {
  out <- fit
  out$cell_id <- cell$cell_id
  tt <- fit$threshold_draws[[endpoint]]
  out$tstar_draws$tstar_days <- tt
  out$tstar_draws$crossed <- is.finite(tt)
  validate_fit_result(out)
}

pilot_score_fit <- function(fit, cell, truth_days, cfg, endpoint) {
  s <- summarise_tstar(fit, "individual", cfg$fitting$interval_level)
  width <- if (is.infinite(s$tstar_lo) && is.infinite(s$tstar_hi)) {
    NA_real_
  } else s$tstar_hi - s$tstar_lo
  bias <- s$tstar_med - truth_days
  row <- data.frame(
    cell_id = cell$cell_id, replicate = cell$replicate, method = fit$method,
    model = fit$model, estimand = "individual", truth_model = cell$truth_model,
    followup_days = cell$followup_days,
    visits_per_year = cell$visits_per_year,
    n_participants = cell$n_participants, sigma_log = cell$sigma_log,
    heterogeneity = cell$heterogeneity, c_thr = cell$c_thr,
    tstar_true_days = truth_days,
    tstar_med = s$tstar_med, tstar_lo = s$tstar_lo, tstar_hi = s$tstar_hi,
    interval_level = s$interval_level, bias = bias, abs_err = abs(bias),
    sq_err = bias^2, covered = s$tstar_lo <= truth_days &&
      s$tstar_hi >= truth_days, width = width, rel_width = width / truth_days,
    prob_no_cross = s$prob_no_cross, adequate = NA,
    prior_set = fit$prior_set, endpoint = endpoint,
    converged = fit$diagnostics$converged,
    rhat_max = fit$diagnostics$rhat_max,
    ess_bulk_min = fit$diagnostics$ess_bulk_min,
    ess_tail_min = fit$diagnostics$ess_tail_min,
    n_obs = ncol(fit$log_lik), runtime_s = fit$diagnostics$runtime_s
  )
  score_rows(row)
}

.pilot_future_score <- function(fit, review, full, next_day) {
  obs <- full$obs[full$obs$time_days > review$cell$followup_days &
                    full$obs$time_days <= next_day, , drop = FALSE]
  if (!nrow(obs)) return(data.frame())
  mu <- early_log_curve(fit$native_draws, fit$model, obs$time_days)
  ll <- .early_loglik(mu, fit$native_draws$sigma, obs, review$cell)
  lpd <- apply(ll, 2, function(z) {
    m <- max(z)
    m + log(mean(exp(z - m)))
  })
  data.frame(model = fit$model, prior_set = fit$prior_set,
             review_days = review$cell$followup_days,
             time_days = obs$time_days, log_predictive_density = lpd,
             converged = fit$diagnostics$converged)
}

run_prior_pilot <- function(cfg = read_config(),
                            cache_dir = "cache/prior-pilot",
                            results_dir = "results/prior-pilot",
                            refit = FALSE) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  pilot <- cfg$pilot
  review_days <- cfg$design$followup_days
  base_cell <- design_cell(
    followup_days = max(review_days),
    visits_per_year = pilot$visits_per_year,
    n_participants = pilot$n_participants,
    sigma_log = cfg$calibration$observation$sigma_log,
    c_thr = cfg$calibration$c_thr, replicate = pilot$replicate,
    base_seed = cfg$run$seed
  )
  validate_design_cell(base_cell)
  full <- simulate_full(base_cell, cfg, horizon_days = max(review_days))
  complete <- reduce_to_review(full, base_cell)
  .pilot_atomic_save(complete, file.path(results_dir, "complete-history.rds"))
  utils::write.csv(complete$obs, file.path(results_dir, "observations.csv"),
                   row.names = FALSE)
  b <- cfg$calibration$biphasic
  truth <- vapply(cfg$calibration$thresholds, function(thr) {
    tstar("biphasic", list(c0 = b$c0, h1 = b$h1_days, h2 = b$h2_days,
                           ts = b$ts_days), thr)
  }, numeric(1))
  fits <- list()
  scores <- list()
  future <- list()
  comparisons <- list()
  attempts <- list()
  failures <- list()
  implementation <- tools::md5sum(c("R/priors.R", "R/fit-mcmc.R"))
  for (day in review_days) {
    cell <- design_cell(
      followup_days = day, visits_per_year = pilot$visits_per_year,
      n_participants = pilot$n_participants, sigma_log = base_cell$sigma_log,
      c_thr = base_cell$c_thr, replicate = base_cell$replicate,
      base_seed = cfg$run$seed
    )
    data <- reduce_to_review(full, cell)
    for (prior_set in pilot$prior_sets) {
      for (model in pilot$models) {
        key <- sprintf("%s-%s-%dd", model, prior_set, day)
        path <- file.path(cache_dir, paste0(key, ".rds"))
        signature <- list(cfg = cfg, obs = data$obs,
                           implementation = implementation)
        cached <- if (!refit && file.exists(path)) readRDS(path) else NULL
        if (!is.null(cached) && identical(cached$signature, signature)) {
          fit <- cached$fit
          cat("Reused ", key, "\n", sep = "")
        } else {
          fit <- NULL
          sampler <- cfg$fitting$sampler
          for (attempt in 1:3) {
            cat("Fitting ", key, " attempt ", attempt, ": ",
                sampler$draws, " draws/chain\n", sep = "")
            result <- tryCatch(fit_mcmc(
              data, model, cfg, prior_set = prior_set, sampler = sampler,
              work_dir = file.path(cache_dir, key, paste0("attempt", attempt))
            ), error = function(e) e)
            if (inherits(result, "error")) {
              failures[[length(failures) + 1L]] <- data.frame(
                key = key, attempt = attempt, message = conditionMessage(result)
              )
              cat("Failed: ", conditionMessage(result), "\n", sep = "")
              break
            }
            fit <- result
            attempts[[length(attempts) + 1L]] <- data.frame(
              key = key, attempt = attempt, draws_per_chain = sampler$draws,
              converged = fit$diagnostics$converged,
              rhat_max = fit$diagnostics$rhat_max,
              ess_bulk_min = fit$diagnostics$ess_bulk_min,
              ess_tail_min = fit$diagnostics$ess_tail_min
            )
            .pilot_atomic_save(list(signature = signature, fit = fit), path)
            if (fit$diagnostics$converged) break
            sampler$draws <- sampler$draws * 4L
            sampler$warmup <- sampler$warmup * 2L
          }
        }
        if (is.null(fit)) next
        fits[[key]] <- fit
        comparisons[[key]] <- cbind(
          followup_days = day, compare_models(list(fit), cfg)
        )
        for (endpoint in names(cfg$calibration$thresholds)) {
          endpoint_cell <- design_cell(
            followup_days = day, visits_per_year = pilot$visits_per_year,
            n_participants = pilot$n_participants,
            sigma_log = base_cell$sigma_log,
            c_thr = cfg$calibration$thresholds[[endpoint]],
            replicate = base_cell$replicate, base_seed = cfg$run$seed
          )
          endpoint_fit <- .pilot_endpoint_fit(fit, endpoint, cfg, endpoint_cell)
          # Include the prior scenario in result keys, allowing a combined
          # score table without changing the shared schema's uniqueness key.
          endpoint_cell$cell_id <- paste0(endpoint_cell$cell_id,
                                          "-prior-", prior_set)
          endpoint_fit$cell_id <- endpoint_cell$cell_id
          validate_design_cell(endpoint_cell)
          validate_fit_result(endpoint_fit)
          scores[[paste(key, endpoint)]] <- pilot_score_fit(
            endpoint_fit, endpoint_cell, truth[[endpoint]], cfg, endpoint
          )
        }
        next_day <- review_days[which(review_days == day) + 1L]
        if (!is.na(next_day)) {
          future[[key]] <- .pilot_future_score(fit, data, complete, next_day)
        }
        # Persist each completed fit's summaries, even if a later fit fails.
        utils::write.csv(score_rows(do.call(rbind, scores)),
                         file.path(results_dir, "duration-summary.csv"),
                         row.names = FALSE)
        utils::write.csv(do.call(rbind, comparisons),
                         file.path(results_dir, "model-comparison.csv"),
                         row.names = FALSE)
        if (length(attempts)) utils::write.csv(do.call(rbind, attempts),
          file.path(results_dir, "sampling-attempts.csv"), row.names = FALSE)
        if (length(failures)) utils::write.csv(do.call(rbind, failures),
          file.path(results_dir, "failures.csv"), row.names = FALSE)
      }
    }
  }
  out <- list(cfg = cfg, complete = complete, truth_days = truth,
               fits = fits, scores = do.call(rbind, scores),
               comparison = do.call(rbind, comparisons),
               future = do.call(rbind, future),
               attempts = do.call(rbind, attempts),
               failures = if (length(failures)) do.call(rbind, failures) else
                 data.frame(key = character(), attempt = integer(),
                            message = character()))
  utils::write.csv(out$failures, file.path(results_dir, "failures.csv"),
                   row.names = FALSE)
  if (length(future)) utils::write.csv(out$future,
    file.path(results_dir, "future-predictions.csv"), row.names = FALSE)
  .pilot_atomic_save(out, file.path(results_dir, "pilot.rds"))
  out
}
