# ---------------------------------------------------------------------------
# Simulation: truth models, observation process, and reduction to a review.
#
# OWNER: workstream (b), simulation and inference. Everyone downstream consumes
# what this produces, so changes to the returned shape go through group review.
# Contract in INTERFACES.md.
#
# The central idea (handout Section 3): simulate one complete antibody history
# per participant, then recreate what a trial review would have seen by hiding
# everything after the review date. Keeping `simulate_full()` independent of
# the follow-up length is what makes the follow-up comparison clean -- the same
# underlying history is truncated at 3, 6, 12 months, so differences between
# reviews are information, not noise.
#
# The biphasic generator lives here; the mechanistic generators live in
# R/calibrate.R and belong to workstream (c). Both return a `sim_dataset()`
# over the full horizon, so `reduce_to_review()` treats them alike.
# ---------------------------------------------------------------------------

#' Simulate complete antibody histories for one design cell.
#'
#' Deterministic given `cell$seed`. Must not depend on `cell$followup_days` or
#' `cell$visits_per_year`: those belong to [reduce_to_review()].
#'
#' @param cell a [design_cell()].
#' @param cfg configuration from [read_config()].
#' @param horizon_days length of the simulated history, long enough that every
#'   participant's true duration is resolved.
#' @param dense_visits_per_year sampling density of the stored history.
#' @return a [sim_dataset()] over the full horizon, passing
#'   [validate_sim_dataset()]. `$truth$tstar_days` comes from [tstar()] and is
#'   `Inf` where a participant never crosses; `$meta$population_tstar_days`
#'   holds the population target.
simulate_full <- function(cell, cfg, horizon_days = 1200,
                          dense_visits_per_year = 52) {
  validate_design_cell(cell)
  if (cell$truth_model != "biphasic" || cell$heterogeneity) {
    stop("simulate_full(): this initial implementation supports ",
         "homogeneous biphasic truth only")
  }
  if (horizon_days < cell$followup_days) {
    stop("horizon_days must include the last review")
  }
  # Store daily measurements, plus the requested dense grid. Reduction can
  # then select rounded-day visits without moving their dates or redrawing.
  times <- sort(unique(c(seq(0, horizon_days, by = 1), horizon_days,
                         round(seq(0, horizon_days,
                                   by = 365.25 / dense_visits_per_year)))))
  times <- times[times <= horizon_days]
  b <- cfg$calibration$biphasic
  latent <- curve_biphasic(times, b$c0, b$h1_days, b$h2_days, b$ts_days)
  # Use the complete-history design seed so changing a review never changes
  # the history. The driver supplies this same full object to every review.
  obs <- .with_seed(cell$seed, do.call(rbind, lapply(
    seq_len(cell$n_participants), function(i) {
      o <- observe(latent, cell$sigma_log)
      data.frame(participant_id = as.integer(i), time_days = times,
                 y_obs = o$y_obs, censor = o$censor)
    }
  )))
  truth <- data.frame(
    participant_id = seq_len(cell$n_participants), c0 = b$c0,
    h1 = b$h1_days, h2 = b$h2_days, ts = b$ts_days,
    tstar_days = tstar("biphasic", list(c0 = b$c0, h1 = b$h1_days,
                                      h2 = b$h2_days, ts = b$ts_days),
                        cell$c_thr)
  )
  full_cell <- cell
  full_cell$followup_days <- as.numeric(horizon_days)
  full_cell$visits_per_year <- as.numeric(dense_visits_per_year)
  validate_sim_dataset(sim_dataset(
    full_cell, obs, truth, truth$tstar_days[1],
    extra_meta = list(mock = FALSE, complete_history = TRUE)
  ))
}

#' Reduce a complete history to the data available at a trial review.
#'
#' Truncate at `cell$followup_days`, thin to `cell$visits_per_year`, and apply
#' the assay limits in `cell$lloq` / `cell$uloq` as censoring.
#'
#' @param full the output of [simulate_full()].
#' @param cell the [design_cell()] describing the review.
#' @return a [sim_dataset()] carrying `cell`, passing [validate_sim_dataset()].
#'   `$truth` is unchanged: the truth does not depend on what was observed.
reduce_to_review <- function(full, cell) {
  validate_sim_dataset(full)
  validate_design_cell(cell)
  if (cell$n_participants != full$cell$n_participants ||
      cell$followup_days > full$cell$followup_days) {
    stop("Review is not contained in the complete history")
  }
  times <- unique(round(seq(0, cell$followup_days + 0.5,
                            by = 365.25 / cell$visits_per_year)))
  times <- times[times <= cell$followup_days]
  if (!all(times %in% full$obs$time_days)) {
    stop("Complete history does not contain the planned visits")
  }
  obs <- full$obs[full$obs$time_days %in% times, , drop = FALSE]
  if (!is.na(cell$lloq)) {
    left <- obs$y_obs < cell$lloq
    obs$censor[left] <- "left"
    obs$y_obs[left] <- NA_real_
  }
  if (!is.na(cell$uloq)) {
    right <- !is.na(obs$y_obs) & obs$y_obs > cell$uloq
    obs$censor[right] <- "right"
    obs$y_obs[right] <- NA_real_
  }
  validate_sim_dataset(sim_dataset(cell, obs, full$truth,
                                   full$meta$population_tstar_days,
                                   extra_meta = list(mock = FALSE)))
}

#' Simulate and reduce in one step.
#'
#' The convenience path most callers want.
make_dataset <- function(cell, cfg, ...) {
  reduce_to_review(simulate_full(cell, cfg, ...), cell)
}

#' Apply the log-normal observation model of handout Equation (5).
#'
#' `log y = log c(t) + e`, `e ~ N(0, sigma_log^2)`, then censoring at the
#' quantification limits. Censored rows carry `NA` in `y_obs`; see
#' [sim_dataset()] for why the limit is not substituted for the value.
observe <- function(c_true, sigma_log, lloq = NA_real_, uloq = NA_real_) {
  if (any(!is.finite(c_true) | c_true <= 0) || sigma_log <= 0) {
    stop("observe(): positive titres and log-error SD are required")
  }
  y <- exp(log(c_true) + rnorm(length(c_true), 0, sigma_log))
  censor <- rep("none", length(y))
  if (!is.na(lloq)) censor[y < lloq] <- "left"
  if (!is.na(uloq)) censor[y > uloq] <- "right"
  y[censor != "none"] <- NA_real_
  data.frame(y_obs = y, censor = censor)
}

#' Draw correlated per-participant random effects on `c0` and `h2`.
#'
#' Handout Equation (6): log-normal effects with log-SDs `omega_0`, `omega_h`
#' and correlation `rho_0h`, all read from `cfg$calibration$heterogeneity`.
#'
#' @return a data frame of `participant_id`, `eta_0`, `eta_h`.
draw_random_effects <- function(n, cfg) {
  .not_implemented("draw_random_effects", "(b) simulation and inference")
}
