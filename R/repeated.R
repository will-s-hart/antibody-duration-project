# Repeated individual fits and scoring from handout Section 6.
# The general workstream driver/scoring interfaces remain separate.

duration_crps <- function(draws, truth) {
  stopifnot(length(draws) > 0, !anyNA(draws), all(draws >= 0),
              length(truth) == 1L, is.finite(truth), truth >= 0)
  if (any(is.infinite(draws))) return(Inf)
  x <- sort(draws)
  knots <- sort(unique(c(0, x, truth)))
  if (length(knots) == 1L) return(0)
  left <- utils::head(knots, -1)
  # Non-negative CDF integrals remain stable even for very long finite tails.
  discrepancy <- findInterval(left, x) / length(x) - (left >= truth)
  sum(diff(knots) * discrepancy^2)
}

duration_score_transform <- function(x, scale_days = 365.25) {
  stopifnot(!anyNA(x), all(x >= 0), is.finite(scale_days), scale_days > 0)
  1 / (1 + scale_days / x)
}

repeated_score_fit <- function(fit, data, cfg, truth_days) {
  validate_fit_result(fit)
  rows <- lapply(names(cfg$calibration$thresholds), function(endpoint) {
    cell <- data$cell
    cell$c_thr <- cfg$calibration$thresholds[[endpoint]]
    cell$cell_id <- paste0(cell$cell_id, "-endpoint-", endpoint)
    validate_design_cell(cell)
    endpoint_fit <- .pilot_endpoint_fit(fit, endpoint, cfg, cell)
    s <- pilot_score_fit(endpoint_fit, cell, truth_days[[endpoint]],
                           cfg, endpoint)
    d <- fit$threshold_draws[[endpoint]]
    r <- cfg$repeated
    g <- function(x) duration_score_transform(x, r$transform_scale_days)
    s$bias_bounded <- g(s$tstar_med) - g(truth_days[[endpoint]])
    s$width_bounded <- g(s$tstar_hi) - g(s$tstar_lo)
    s$unbounded_interval <- is.infinite(s$tstar_lo) |
      is.infinite(s$tstar_hi)
    s$infinite_median <- is.infinite(s$tstar_med)
    s$crps_bounded <- duration_crps(g(d), g(truth_days[[endpoint]]))
    s$crps_unrestricted_empirical_days <- duration_crps(d,
                                                       truth_days[[endpoint]])
    horizons <- c(r$crps_horizon_days,
                    unlist(r$sensitivity_horizons_days))
    for (j in seq_along(horizons)) {
      name <- c("crps_restricted_days", "crps_short_days",
                   "crps_long_days")[j]
      s[[name]] <- duration_crps(pmin(d, horizons[j]),
                                  min(truth_days[[endpoint]], horizons[j]))
    }
    score_rows(s)
  })
  out <- score_rows(dplyr::bind_rows(rows))
  numeric_cols <- vapply(out, is.numeric, logical(1))
  if (any(vapply(out[numeric_cols], function(x) any(is.nan(x)), logical(1)))) {
    stop("Undefined repeated-study score")
  }
  out
}

.repeated_cell <- function(cfg, day, replicate) {
  design_cell(followup_days = day,
    visits_per_year = cfg$pilot$visits_per_year, n_participants = 1L,
    sigma_log = cfg$calibration$observation$sigma_log,
    c_thr = cfg$calibration$c_thr, replicate = replicate,
    base_seed = cfg$run$seed)
}

run_repeated_replicate <- function(replicate, cfg = read_config(),
                                   cache_dir = "cache/repeated-default",
                                   results_dir = "results/repeated-default") {
  stopifnot(cfg$pilot$n_participants == 1L)
  days <- cfg$design$followup_days
  cell <- .repeated_cell(cfg, max(days), replicate)
  full <- simulate_full(cell, cfg, horizon_days = max(days))
  complete <- reduce_to_review(full, cell)
  validate_sim_dataset(complete)
  truth <- vapply(cfg$calibration$thresholds, function(thr) {
    b <- cfg$calibration$biphasic
    tstar("biphasic", list(c0 = b$c0, h1 = b$h1_days,
                            h2 = b$h2_days, ts = b$ts_days), thr)
  }, numeric(1))
  code <- tools::md5sum(c("R/priors.R", "R/fit-mcmc.R", "R/repeated.R",
                           "R/pilot.R"))
  scores <- list()
  statuses <- list()
  all_attempts <- list()
  for (day in days) {
    data <- reduce_to_review(full, .repeated_cell(cfg, day, replicate))
    for (model in unlist(cfg$repeated$models)) {
      key <- sprintf("rep-%03d-%s-%dd", replicate, model, day)
      path <- file.path(cache_dir, paste0(key, ".rds"))
      signature <- list(cell = data$cell, obs = data$obs,
        calibration = cfg$calibration, fitting = cfg$fitting,
        scoring = cfg$repeated[c("prior_set", "transform_scale_days",
          "crps_horizon_days", "sensitivity_horizons_days")], code = code)
      cached <- if (file.exists(path)) readRDS(path) else NULL
      if (!is.null(cached) && identical(cached$signature, signature)) {
        result <- cached$result
      } else {
        fit <- NULL
        attempts <- list()
        sampler <- cfg$fitting$sampler
        error_message <- ""
        for (attempt in 1:3) {
          cat(key, " attempt ", attempt, "\n", sep = "")
          candidate <- tryCatch(fit_mcmc(data, model, cfg,
            prior_set = cfg$repeated$prior_set, sampler = sampler,
            compute_loo = FALSE,
            work_dir = file.path(cache_dir, key, paste0("attempt", attempt))),
            error = function(e) e)
          error_message <- if (inherits(candidate, "error")) {
            conditionMessage(candidate)
          } else ""
          attempts[[attempt]] <- data.frame(replicate = replicate,
            model = model, followup_days = day, attempt = attempt,
            draws_per_chain = sampler$draws,
            status = if (nzchar(error_message)) "error" else
              if (candidate$diagnostics$converged) "converged" else
                "unconverged",
            message = error_message)
          if (nzchar(error_message)) break
          fit <- candidate
          if (fit$diagnostics$converged) break
          sampler$draws <- sampler$draws * 4L
          sampler$warmup <- sampler$warmup * 2L
        }
        status <- data.frame(replicate = replicate, model = model,
          followup_days = day, prior_set = cfg$repeated$prior_set,
          completed = !is.null(fit),
          converged = !is.null(fit) && fit$diagnostics$converged,
          status = if (is.null(fit)) "error" else
            if (fit$diagnostics$converged) "converged" else "unconverged",
          rhat_max = if (is.null(fit)) NA_real_ else fit$diagnostics$rhat_max,
          ess_bulk_min = if (is.null(fit)) NA_real_ else
            fit$diagnostics$ess_bulk_min,
          message = error_message)
        s <- if (is.null(fit)) NULL else {
          repeated_score_fit(fit, data, cfg, truth)
        }
        if (!is.null(fit)) {
          # Preserve validated native parameter/duration draws and raw CODA,
          # without storing duplicate wide draws or likelihood matrices.
          fit[c("posterior", "native_draws", "log_lik", "loo")] <- NULL
          validate_fit_result(fit)
        }
        result <- list(status = status, scores = s,
          attempts = dplyr::bind_rows(attempts), fit = fit)
        .pilot_atomic_save(list(signature = signature, result = result), path)
      }
      statuses[[key]] <- result$status
      if (!is.null(result$scores)) scores[[key]] <- result$scores
      all_attempts[[key]] <- result$attempts
    }
  }
  out <- list(replicate = replicate, cfg = cfg, observations = complete$obs,
    truth_days = truth, status = dplyr::bind_rows(statuses),
    scores = score_rows(dplyr::bind_rows(scores)),
    attempts = dplyr::bind_rows(all_attempts))
  path <- file.path(results_dir, "replicates",
                      sprintf("rep-%03d.rds", replicate))
  .pilot_atomic_save(out, path)
  cat("Completed replicate ", replicate, "\n", sep = "")
  path
}

repeated_wilson <- function(success, n, level = .95) {
  if (n == 0L) return(c(lo = NA_real_, hi = NA_real_))
  z <- stats::qnorm(1 - (1-level)/2)
  p <- success/n
  centre <- (p+z^2/(2*n))/(1+z^2/n)
  half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n)
  c(lo = max(0, centre-half), hi = min(1, centre+half))
}

.repeated_mean <- function(x) {
  if (!length(x)) return(c(mean = NA_real_, mcse = NA_real_))
  c(mean = mean(x), mcse = if (length(x) > 1L && all(is.finite(x))) {
    stats::sd(x)/sqrt(length(x))
  } else NA_real_)
}

summarise_repeated_study <- function(scores, status, cfg) {
  validate_score_rows(scores)
  output <- list()
  for (model in unlist(cfg$repeated$models)) {
    for (day in cfg$design$followup_days) {
      st <- status[status$model == model & status$followup_days == day, ]
      n <- sum(st$converged)
      planned <- cfg$repeated$n_replicates
      for (endpoint in names(cfg$calibration$thresholds)) {
        s <- scores[scores$model == model & scores$followup_days == day &
                      scores$endpoint == endpoint & scores$converged, ]
        stopifnot(nrow(s) == n, nrow(st) == planned)
        covered <- sum(s$covered)
        wi <- repeated_wilson(covered, n)
        row <- data.frame(model = model, followup_days = day,
          endpoint = endpoint, prior_set = cfg$repeated$prior_set,
          n_planned = planned, n_completed = sum(st$completed), n_ok = n,
          n_failed = planned-n, n_covered = covered,
          coverage = if (n) covered/n else NA_real_,
          coverage_mcse = if (n) sqrt(covered/n*(1-covered/n)/n) else NA_real_,
          coverage_lo = wi[1], coverage_hi = wi[2],
          coverage_all_lower = covered/planned,
          coverage_all_upper = (covered+planned-n)/planned)
        metrics <- list(bias_days = s$bias,
          bias_bounded = s$bias_bounded, width_bounded = s$width_bounded,
          prob_unbounded = as.numeric(s$unbounded_interval),
          prob_infinite_median = as.numeric(s$infinite_median),
          prob_no_cross = s$prob_no_cross,
          crps_bounded = s$crps_bounded,
          crps_restricted_days = s$crps_restricted_days,
          crps_short_days = s$crps_short_days,
          crps_long_days = s$crps_long_days,
          crps_unrestricted_empirical_days =
            s$crps_unrestricted_empirical_days)
        for (nm in names(metrics)) {
          ms <- .repeated_mean(metrics[[nm]])
          row[[nm]] <- ms[1]
          row[[paste0(nm, "_mcse")]] <- ms[2]
        }
        widths <- ifelse(s$unbounded_interval, Inf, s$width)
        row$median_width_days <- if (n) {
          as.numeric(quantile(widths, .5, type = 1))
        } else NA_real_
        row$unbounded_lo <- repeated_wilson(sum(s$unbounded_interval), n)[1]
        row$unbounded_hi <- repeated_wilson(sum(s$unbounded_interval), n)[2]
        output[[length(output)+1L]] <- row
      }
    }
  }
  out <- dplyr::bind_rows(output)
  if (any(vapply(out[vapply(out, is.numeric, logical(1))],
                    function(x) any(is.nan(x)), logical(1)))) {
    stop("Undefined repeated-study summary")
  }
  out
}

collect_repeated_study <- function(cfg = read_config(),
                                    results_dir = "results/repeated-default") {
  paths <- file.path(results_dir, "replicates", sprintf("rep-%03d.rds",
    seq_len(cfg$repeated$n_replicates)))
  if (!all(file.exists(paths))) stop("Replicates are not all complete")
  reps <- lapply(paths, readRDS)
  status <- dplyr::bind_rows(lapply(reps, `[[`, "status"))
  scores <- score_rows(dplyr::bind_rows(lapply(reps, `[[`, "scores")))
  attempts <- dplyr::bind_rows(lapply(reps, `[[`, "attempts"))
  summary <- summarise_repeated_study(scores, status, cfg)
  out <- list(cfg = cfg, scores = scores, status = status, attempts = attempts,
    summary = summary, truth_days = reps[[1]]$truth_days)
  for (name in c("scores", "status", "attempts", "summary")) {
    utils::write.csv(out[[name]], file.path(results_dir, paste0(name, ".csv")),
                       row.names = FALSE)
  }
  .pilot_atomic_save(out, file.path(results_dir, "study.rds"))
  capture.output(sessionInfo(),
                   file = file.path(results_dir, "session-info.txt"))
  out
}
