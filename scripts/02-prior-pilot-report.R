#!/usr/bin/env Rscript
# Rebuild figures and a self-contained HTML report from completed pilot fits.
suppressMessages(pkgload::load_all(quiet = TRUE))
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
anchored <- "--two-anchor" %in% args
out_dir <- if (anchored) "results/prior-pilot-two-anchor" else
  "results/prior-pilot"
pilot <- readRDS(file.path(out_dir, "pilot.rds"))
cfg <- pilot$cfg
prior_levels <- unlist(cfg$pilot$prior_sets)
colours <- c(exponential = "#C44335", biphasic = "#176B9B",
             powerlaw = "#7B5196", plateau = "#2D8263")
theme_set(theme_minimal(base_size = 12))
times <- sort(unique(c(seq(0, 2200, by = 15), 90, 365)))

curve_summary <- function(params, model, prior_set, followup_days = NA_real_) {
  log_mu <- early_log_curve(params, model, times)
  q <- apply(log_mu, 2, quantile, c(.05, .5, .95))
  data.frame(time_days = times, lo = pmax(exp(q[1, ]), 1e-300),
             med = pmax(exp(q[2, ]), 1e-300),
             hi = exp(q[3, ]), model = model, prior_set = prior_set,
             followup_days = followup_days)
}

prior_curves <- list()
prior_durations <- list()
anchor_checks <- list()
for (set in cfg$pilot$prior_sets) {
  # Identical RNG streams give identical common draws across all families.
  seed <- derive_seed(paste("prior-predictive", set), cfg$pilot$replicate,
                      cfg$run$seed)
  for (model in cfg$pilot$models) {
    free <- draw_early_prior(model, cfg, set, 4000, seed)
    p <- early_prior(cfg, set)
    pars <- early_transform(free, model, p$reference_day,
                              if (anchored) p$later_day else 365)
    key <- paste(set, model)
    prior_curves[[key]] <- curve_summary(pars, model, set)
    if (anchored && model != "exponential") {
      mu <- early_log_curve(pars, model, c(0, 90, 365))
      anchor_checks[[key]] <- data.frame(
        prior_set = set, model = model,
        max_day90_log_error = max(abs(mu[, 2] - free$log_c0 + pars$D)),
        max_day365_log_error = max(abs(mu[, 3] - free$log_c0 + pars$D_b)),
        prob_f_above_half = if (model == "biphasic") mean(pars$f > .5)
          else NA_real_
      )
    }
    for (endpoint in names(cfg$calibration$thresholds)) {
      log_tt <- early_log_tstar(pars, model,
                                 cfg$calibration$thresholds[[endpoint]])
      # Type 1 preserves atoms at zero/Inf and avoids Inf interpolation.
      q <- exp(quantile(log_tt, c(.05, .5, .95), names = FALSE, type = 1))
      prior_durations[[paste(key, endpoint)]] <- data.frame(
        prior_set = set, model = model, endpoint = endpoint,
        lo_days = q[1], median_days = q[2], hi_days = q[3],
        prob_no_cross = mean(log_tt == Inf),
        prob_finite_beyond_numeric_range = mean(is.finite(log_tt) &
          log_tt > log(.Machine$double.xmax))
      )
    }
  }
}
prior_curves <- bind_rows(prior_curves)
prior_curves$prior_set <- factor(prior_curves$prior_set,
                                 levels = prior_levels)
prior_curves$model <- factor(prior_curves$model, levels = names(colours))
prior_durations <- bind_rows(prior_durations)
write.csv(prior_durations, file.path(out_dir, "prior-duration-summary.csv"),
           row.names = FALSE)
if (length(anchor_checks)) {
  write.csv(bind_rows(anchor_checks), file.path(out_dir, "anchor-checks.csv"),
             row.names = FALSE)
}
thresholds <- data.frame(
  titre = unlist(cfg$calibration$thresholds),
  endpoint = names(cfg$calibration$thresholds)
)
prior_plot <- ggplot(prior_curves, aes(time_days, med)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = .2) +
  geom_line(aes(colour = model), linewidth = .7) +
  geom_hline(data = thresholds, aes(yintercept = titre),
             linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = if (anchored) c(90, 365) else 90,
               linetype = "dotted", colour = "grey55") +
  facet_grid(prior_set ~ model) + scale_y_log10() +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  coord_cartesian(ylim = c(.001, 5)) +
  labs(title = "Common early priors, different tail assumptions",
        subtitle = "Latent titre: median and 90% prior interval",
        x = "Days since initial measurement", y = "Convalescent-relative titre",
        caption = paste0(
          if (anchored) "Dotted lines: days 90 and 365. " else
            "Dotted line: day 90. ", "Dashed lines: 0.202 and 0.030 ",
          "thresholds."
        )) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "prior-predictive.png"), prior_plot,
        width = 13, height = 7, dpi = 160)

scores <- pilot$scores
scores$prior_set <- factor(scores$prior_set, levels = prior_levels)
scores$endpoint <- factor(scores$endpoint, levels = c("symptomatic", "severe"))
scores$review_months <- cfg$design$followup_months[
  match(scores$followup_days, cfg$design$followup_days)
]
offset <- setNames(c(-.6, -.2, .2, .6), names(colours))
scores$x <- scores$review_months + offset[scores$model]
cap <- 7305
scores$plot_lo <- pmax(1, pmin(scores$tstar_lo, cap))
scores$plot_hi <- pmin(scores$tstar_hi, cap)
scores$plot_med <- pmin(scores$tstar_med, cap)
scores$upper_label <- ifelse(is.infinite(scores$tstar_hi), "Inf",
                              ifelse(scores$tstar_hi > cap, ">", ""))
truth <- data.frame(endpoint = names(pilot$truth_days),
                     truth_days = as.numeric(pilot$truth_days))
truth$endpoint <- factor(truth$endpoint, levels = c("symptomatic", "severe"))
duration_plot <- ggplot(scores, aes(x, plot_med, colour = model)) +
  geom_hline(data = truth, aes(yintercept = truth_days),
             linetype = "dashed", colour = "grey35") +
  geom_linerange(aes(ymin = plot_lo, ymax = plot_hi), linewidth = .6) +
  geom_point(data = scores[is.finite(scores$tstar_med), ],
               aes(shape = converged), size = 2) +
  geom_text(aes(y = plot_hi, label = upper_label), vjust = -.5, size = 2.6,
              show.legend = FALSE) +
  facet_grid(prior_set ~ endpoint) + scale_y_log10() +
  scale_colour_manual(values = colours) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4)) +
  guides(shape = "none") +
  scale_x_continuous(breaks = cfg$design$followup_months) +
  coord_cartesian(ylim = c(70, cap * 1.4), clip = "off") +
  labs(title = "Individual duration estimates as follow-up accumulates",
        subtitle = "Posterior medians and 90% intervals; one simulated person",
        x = "Available follow-up (months)",
        y = "Threshold-crossing time (days; logarithmic axis)",
        caption = paste("Dashed lines: true crossings. Inf denotes an infinite",
                        "upper limit; > denotes a finite limit beyond",
                        "20 years.")) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "duration-estimates.png"), duration_plot,
        width = 12, height = 8, dpi = 160)

chosen_days <- cfg$design$followup_days[c(1, 2, 4, 6)]
posterior_curves <- bind_rows(lapply(names(pilot$fits), function(key) {
  fit <- pilot$fits[[key]]
  day <- as.numeric(sub(".*-([0-9]+)d$", "\\1", key))
  if (!day %in% chosen_days) return(NULL)
  curve_summary(fit$native_draws, fit$model, fit$prior_set, day)
}))
b <- cfg$calibration$biphasic
truth_curve <- data.frame(time_days = times,
                          titre = curve_biphasic(times, b$c0, b$h1_days,
                                                  b$h2_days, b$ts_days))
observations <- bind_rows(lapply(chosen_days, function(day) {
  obs <- pilot$complete$obs[pilot$complete$obs$time_days <= day, ]
  obs$followup_days <- day
  obs
}))
posterior_curves$review <- factor(posterior_curves$followup_days,
  levels = chosen_days,
  labels = paste(cfg$design$followup_months[c(1, 2, 4, 6)], "months"))
posterior_curves$prior_set <- factor(posterior_curves$prior_set,
                                     levels = prior_levels)
observations$review <- factor(observations$followup_days, levels = chosen_days,
  labels = levels(posterior_curves$review))
trajectory_plot <- ggplot(posterior_curves,
                          aes(time_days, med, colour = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = .09,
               colour = NA) + geom_line(linewidth = .65) +
  geom_line(data = truth_curve, aes(time_days, titre), inherit.aes = FALSE,
             colour = "black", linetype = "dashed", linewidth = .7) +
  geom_point(data = observations, aes(time_days, y_obs), inherit.aes = FALSE,
              colour = "black", size = 1.1) +
  geom_hline(data = thresholds, aes(yintercept = titre),
              colour = "grey50", linetype = "dotted") +
  facet_grid(prior_set ~ review) + scale_y_log10() +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  coord_cartesian(ylim = c(.005, 2)) +
  labs(title = "Similar early fits can imply very different persistence",
        subtitle = "Latent titre: posterior median and 90% interval",
        x = "Days since initial measurement", y = "Convalescent-relative titre",
        caption = paste0(
          "Black dashed: truth. Black points: data available at that ",
          "review."
        )) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "posterior-trajectories.png"), trajectory_plot,
        width = 14, height = 7, dpi = 160)

comparison <- pilot$comparison
comparison$review_months <- cfg$design$followup_months[
  match(comparison$followup_days, cfg$design$followup_days)
]
future <- pilot$future |>
  group_by(prior_set, model, review_days) |>
  summarise(n_future = n(), mean_log_score = mean(log_predictive_density),
             .groups = "drop")
write.csv(future, file.path(out_dir, "future-summary.csv"), row.names = FALSE)
diagnostics <- bind_rows(lapply(names(pilot$fits), function(key) {
  f <- pilot$fits[[key]]
  data.frame(key = key, converged = f$diagnostics$converged,
               rhat_max = f$diagnostics$rhat_max,
               ess_bulk_min = f$diagnostics$ess_bulk_min,
               ess_tail_min = f$diagnostics$ess_tail_min,
               draws_per_chain = f$sampler$draws,
               pareto_k_max = f$pareto_k_max, loo_reliable = f$loo_reliable)
}))
write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
yaml::write_yaml(cfg, file.path(out_dir, "configuration.yaml"))

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
html_table <- function(x) {
  fmt <- function(v) {
    if (is.numeric(v)) v <- format(signif(v, 4), trim = TRUE)
    escape_html(as.character(v))
  }
  header <- paste0("<th>", escape_html(names(x)), "</th>", collapse = "")
  rows <- apply(as.data.frame(lapply(x, fmt)), 1, function(row) {
    paste0("<tr><td>", paste(row, collapse = "</td><td>"), "</td></tr>")
  })
  paste0("<div class='table'><table><thead><tr>", header,
          "</tr></thead><tbody>", paste(rows, collapse = ""),
          "</tbody></table></div>")
}
image_data <- function(path) {
  bytes <- as.integer(readBin(path, "raw", file.info(path)$size))
  pad <- (3 - length(bytes) %% 3) %% 3
  m <- matrix(c(bytes, rep(0L, pad)), nrow = 3)
  indices <- as.vector(rbind(m[1, ] %/% 4,
    (m[1, ] %% 4) * 16 + m[2, ] %/% 16,
    (m[2, ] %% 16) * 4 + m[3, ] %/% 64, m[3, ] %% 64))
  alphabet <- strsplit(paste0(LETTERS, collapse = ""), "")[[1]]
  alphabet <- c(alphabet, letters, as.character(0:9), "+", "/")
  encoded <- alphabet[indices + 1L]
  if (pad) encoded[seq.int(length(encoded) - pad + 1L, length(encoded))] <- "="
  paste0("<img src='data:image/png;base64,", paste0(encoded, collapse = ""),
          "' alt='", basename(path), "'>")
}
selected <- scores[scores$followup_days %in% chosen_days,
                    c("review_months", "prior_set", "model", "endpoint",
                       "tstar_med", "tstar_lo", "tstar_hi", "prob_no_cross",
                       "covered", "converged")]
early <- scores[scores$followup_days == min(cfg$design$followup_days) &
                  scores$model == "exponential" & scores$endpoint == "severe",
                  ]
late <- scores[scores$followup_days == max(cfg$design$followup_days) &
                 scores$model == "biphasic" & scores$endpoint == "severe", ]
findings <- c(
  "<h2>What this realisation shows</h2><ul>",
  sprintf(paste0(
    "<li>At three months, the exponential severe-duration ",
    "estimate is %.0f [%.0f, %.0f] days under the default ",
    "prior and %.0f [%.0f, %.0f] under the broader prior. Both ",
    "90%% intervals exclude the true crossing of %.0f ",
    "days.</li>"
  ),
    early$tstar_med[1], early$tstar_lo[1], early$tstar_hi[1],
    early$tstar_med[2], early$tstar_lo[2], early$tstar_hi[2],
    pilot$truth_days["severe"]),
  sprintf(paste0(
    "<li>At 24 months, biphasic severe-duration estimates ",
    "remain uncertain: %.0f [%.0f, %.0f] days with the default ",
    "prior and %.0f [%.0f, %.0f] with the broader prior. Extra ",
    "follow-up does not guarantee monotonic interval narrowing ",
    "in a single noisy history.</li>"
  ),
    late$tstar_med[1], late$tstar_lo[1], late$tstar_hi[1],
    late$tstar_med[2], late$tstar_lo[2], late$tstar_hi[2]),
  paste0(
    "<li>The flexible families agree much more closely on the ",
    "higher crossing after longer follow-up, while differing ",
    "substantially about the lower crossing. The plateau ",
    "retains high severe non-crossing probability, and the ",
    "power-law severe estimate can be very long or ",
    "prior-sensitive at early reviews.</li>"
  ),
  paste0(
    "<li>Broadening the priors preserves the shared early ",
    "information across families. It mainly exposes additional ",
    "uncertainty about the unobserved tail; it does not make ",
    "the exponential's constant decline assumption ",
    "flexible.</li></ul>"
  )
)
if (anchored) {
  findings <- c(
    "<h2>What this realisation shows</h2><ul>",
    sprintf(paste0(
      "<li>The exponential retains its original prior and chains. ",
      "At three months its severe-duration estimate is %.0f [%.0f, %.0f] ",
      "days under the default and %.0f [%.0f, %.0f] under the broader ",
      "prior; the true crossing is %.0f days.</li>"
    ), early$tstar_med[1], early$tstar_lo[1], early$tstar_hi[1],
       early$tstar_med[2], early$tstar_lo[2], early$tstar_hi[2],
       pilot$truth_days["severe"]),
    sprintf(paste0(
      "<li>At 24 months the biphasic severe-duration estimates are ",
      "%.0f [%.0f, %.0f] days under the default and ",
      "%.0f [%.0f, %.0f] under the broader prior. These are individual ",
      "posterior intervals from one noisy history.</li>"
    ), late$tstar_med[1], late$tstar_lo[1], late$tstar_hi[1],
       late$tstar_med[2], late$tstar_lo[2], late$tstar_hi[2]),
    paste0(
      "<li>The three flexible families now share the joint prior on ",
      "initial, day-90 and day-365 latent titres. They need not agree ",
      "between these anchors or about later threshold crossings. ",
      "The day-365 prior supplies information beyond early reviews.</li>",
      "<li>The extra biphasic variable is the slow rate as a fraction ",
      "of the mean log-decay rate between days 90 and 365. The initial ",
      "slow-component weight f is derived and may exceed one half.</li></ul>"
    )
  )
  power_prior <- prior_durations[prior_durations$model == "powerlaw" &
                                   prior_durations$endpoint == "severe", ]
  findings <- c(findings, sprintf(paste0(
    "<p class='note'>The anchor construction is coherent, but these ",
    "exploratory priors are not yet calibrated for plausible long-term ",
    "behaviour. The power-law severe-duration prior has a 95th percentile ",
    "of %.2g days under the default and %.2g under the broader option. ",
    "Matching first-year titres does not equalise information about the ",
    "tail. Keep this as a sensitivity formulation before choosing priors ",
    "for the main simulation study.</p>"
  ), power_prior$hi_days[1], power_prior$hi_days[2]))
}

baseline_html <- character()
baseline_path <- "results/prior-pilot/pilot.rds"
if (anchored && file.exists(baseline_path)) {
  baseline <- readRDS(baseline_path)
  if (!identical(baseline$complete$obs, pilot$complete$obs)) {
    stop("The baseline comparison must use identical observations")
  }
  old <- baseline$scores
  old$prior_strength <- old$prior_set
  old$formulation <- "Early only"
  new <- pilot$scores
  new$prior_strength <- sub("anchored_", "", new$prior_set)
  new$formulation <- "Days 90 and 365"
  both <- bind_rows(old, new)
  both$review_months <- cfg$design$followup_months[
    match(both$followup_days, cfg$design$followup_days)
  ]
  both$endpoint <- factor(both$endpoint, levels = c("symptomatic", "severe"))
  both$prior_strength <- factor(both$prior_strength,
                                 levels = c("default", "broader"))
  both$x <- both$review_months + offset[both$model] +
    ifelse(both$formulation == "Early only", -.09, .09)
  both$plot_med <- pmin(both$tstar_med, cap)
  both$plot_lo <- pmax(1, pmin(both$tstar_lo, cap))
  both$plot_hi <- pmin(both$tstar_hi, cap)
  both$upper_label <- ifelse(is.infinite(both$tstar_hi), "Inf",
                              ifelse(both$tstar_hi > cap, ">", ""))
  label_key <- with(both, paste(prior_strength, endpoint, model,
                                 review_months, upper_label))
  both$upper_label[duplicated(label_key)] <- ""
  comparison_plot <- ggplot(both, aes(x, plot_med, colour = model,
                                      shape = formulation)) +
    geom_hline(data = truth, aes(yintercept = truth_days),
                 inherit.aes = FALSE, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = plot_lo, ymax = plot_hi,
                         linetype = formulation), linewidth = .45) +
    geom_point(data = both[is.finite(both$tstar_med), ], size = 1.8) +
    geom_text(aes(y = plot_hi, label = upper_label), size = 2,
                 vjust = -.5, show.legend = FALSE) +
    facet_grid(prior_strength ~ endpoint) + scale_y_log10() +
    scale_colour_manual(values = colours) +
    scale_x_continuous(breaks = cfg$design$followup_months) +
    coord_cartesian(ylim = c(70, cap*1.4), clip = "off") +
    labs(title = "Effect of matching the second anchor",
           subtitle = paste("Identical observations; individual medians",
                              "and 90% intervals"),
           x = "Available follow-up (months)", y = "Crossing time (days)",
           caption = paste("Inf: infinite upper limit; >: finite upper",
                             "limit beyond 20 years.")) +
    theme(legend.position = "bottom")
  ggsave(file.path(out_dir, "prior-formulation-comparison.png"),
           comparison_plot, width = 14, height = 9, dpi = 160)
  cols <- c("formulation", "prior_strength", "review_months", "model",
              "endpoint", "tstar_med", "tstar_lo", "tstar_hi",
              "prob_no_cross", "converged")
  selected_both <- both[both$review_months %in% c(3, 12, 24), cols]
  write.csv(selected_both,
              file.path(out_dir, "prior-formulation-comparison.csv"),
              row.names = FALSE)
  baseline_html <- c("<h2>Comparison with the early-only formulation</h2>",
    paste0("<p>Both analyses use exactly the same simulated observations. ",
             "Inf is retained; flagged fits require caution.</p>"),
    image_data(file.path(out_dir, "prior-formulation-comparison.png")),
    html_table(selected_both))
}
if (is.null(pilot$failures)) {
  pilot$failures <- data.frame(key = character(), attempt = integer(),
                               message = character())
}
write.csv(pilot$failures, file.path(out_dir, "failures.csv"), row.names = FALSE)
html <- c("<!doctype html><html lang='en'><meta charset='utf-8'>",
  "<title>One-person antibody-duration pilot</title><style>",
  paste0(
    "body{font:16px/1.55 system-ui;max-width:1300px;",
    "margin:40px auto;padding:0 24px;color:#20303e}"
  ),
  paste0(
    "h1,h2{color:#173f5d}img{width:100%;height:auto}",
    "table{border-collapse:collapse;font-size:13px}"
  ),
  paste0(
    "th,td{padding:7px 10px;border-bottom:1px solid #dde3e8;",
    "text-align:left}th{background:#eef3f7}"
  ),
  paste0(
    ".table{overflow:auto;margin:20px 0}",
    ".note{background:#eef3f7;padding:18px;border-left:4px ",
    "solid #176B9B}"
  ),
  "</style><h1>One-person antibody-duration pilot</h1>",
  paste0(
    "<p>First three steps of handout Section 6: recover a ",
    "biphasic truth, reveal longer follow-up, and compare four ",
    "candidate families under default and broader priors.</p>"
  ),
  sprintf(paste0(
    "<p>One simulated individual; %d monthly measurements ",
    "including day zero; log-error SD 0.20. True individual ",
    "crossings: %.1f days (symptomatic anchor) and %.1f days ",
    "(severe anchor).</p>"
  ),
    nrow(pilot$complete$obs), pilot$truth_days["symptomatic"],
    pilot$truth_days["severe"]),
  paste0(
    "<p class='note'>This is one noisy realisation, not a ",
    "coverage study. An interval containing the truth here ",
    "does not establish calibration or trial adequacy. Severe ",
    "protection is represented only by the specified antibody ",
    "threshold.</p>"
  ),
  sprintf(paste0(
    "<p>%d/%d completed fits passed the declared convergence checks ",
    "out of %d planned fits. %d fits failed to produce a result ",
    "(see failures below). %d ",
    "pointwise LOO comparisons have Pareto-k above 0.7 or ",
    "unresolved convergence and are flagged as unreliable.</p>"
  ),
    sum(diagnostics$converged), nrow(diagnostics),
    length(cfg$pilot$models) * length(cfg$pilot$prior_sets) *
      length(cfg$design$followup_days), length(unique(pilot$failures$key)),
    sum(!diagnostics$loo_reliable)),
  findings,
  baseline_html,
  "<h2>Duration estimates</h2>",
  image_data(file.path(out_dir, "duration-estimates.png")),
  "<h2>Posterior trajectories</h2>",
  image_data(file.path(out_dir, "posterior-trajectories.png")),
  "<h2>Prior predictions</h2>",
  paste0(
    "<p>Median and 90% intervals of latent trajectories. Both ",
    "columns preserve the common initial/day-90 prior across ",
    "models; ", if (anchored) paste0(
      "the flexible families additionally share the day-365 anchor. "
    ) else "", "Broadening changes tails and floor probabilities ",
    "as well as dispersion. No future measurements enter these ",
    "priors.</p>"
  ),
  image_data(file.path(out_dir, "prior-predictive.png")),
  html_table(prior_durations),
  if (anchored) paste0(
    "<p>Finite power-law durations can exceed floating-point range ",
    "under the broader prior. Quantiles are computed on log time; the ",
    "separate numeric-range column records that finite tail mass. It ",
    "is not included in the probability of never crossing.</p>"
  ),
  if (anchored) html_table(bind_rows(anchor_checks)),
  "<h2>Selected individual crossing estimates</h2>",
  paste0(
    "<p>All intervals are 90%. Inf is retained as an outcome. ",
    "Covered records whether this single interval contains the ",
    "generating truth.</p>"
  ),
  html_table(selected),
  "<h2>Sampler diagnostics</h2>", html_table(diagnostics),
  if (nrow(pilot$failures)) c("<h2>Failures</h2>",
    paste0("<p>These fits have no duration estimate in the figures or ",
             "tables. Raw JAGS output is retained in the cache; no draws ",
             "were dropped or recoded as never crossing.</p>"),
    html_table(pilot$failures)),
  paste0(
    "<p>Four independent JAGS chains; initially 2000 ",
    "adaptation, 2000 warm-up and 4000 retained draws per ",
    "chain. Failed diagnostics trigger longer runs, with each ",
    "attempt retained. Divergences do not apply to this ",
    "sampler. Rank R-hat must be below 1.01; bulk and ",
    "applicable tail ESS must be at least 400. Constant ",
    "outcome/quantile indicators are marked not applicable.</p>"
  ),
  "<h2>Observed-visit model comparison</h2>",
  paste0(
    "<p>Higher expected log predictive density is better. ",
    "Pointwise PSIS-LOO concerns interpolation among observed ",
    "visits, not extrapolation to the severe threshold. ",
    "Unreliable comparisons must not be used to declare a ",
    "winner.</p>"
  ),
  html_table(comparison),
  "<h2>Withheld subsequent visits</h2>",
  paste0(
    "<p>Mean log predictive density for the visits between a ",
    "review and the next review. Higher is better. This uses ",
    "future observations only for evaluation after fitting. ",
    "There are few observations and one individual.</p>"
  ),
  html_table(future),
  paste0(
    "<h2>Reproduce</h2><pre>Rscript ",
    "scripts/01-prior-pilot.R", if (anchored) " --two-anchor" else "",
    "\nRscript scripts/02-prior-pilot-report.R",
    if (anchored) " --two-anchor" else "", "</pre>"
  ),
  paste0(
    "<p>All source is R. JAGS must be installed separately; R ",
    "dependencies are recorded in renv.lock. Configuration, ",
    "session information, observations, posterior fits and ",
    "sampling attempts are saved alongside this report or in ",
    if (anchored) "cache/prior-pilot-two-anchor. " else "cache/prior-pilot. ",
    "No private archived files are read by ",
    "these scripts.</p></html>"
  ))
writeLines(html, file.path(out_dir, "report.html"))
cat("Wrote ", out_dir, "/report.html and figures.\n", sep = "")
