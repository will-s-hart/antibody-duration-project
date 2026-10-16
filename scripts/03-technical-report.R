#!/usr/bin/env Rscript
# Assemble report figures and small tables from the three completed pilots.
# Run the fitting scripts first; this script never fits or simulates data.
suppressMessages(pkgload::load_all(quiet = TRUE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

paths <- c("results/prior-pilot/pilot.rds",
            "results/prior-pilot-two-anchor/pilot.rds",
            "results/prior-pilot-initial-rate/pilot.rds")
if (!all(file.exists(paths))) stop("Run all three prior pilots first")
pilots <- lapply(paths, readRDS)
stopifnot(all(vapply(pilots, function(p) {
  identical(pilots[[1]]$complete$obs, p$complete$obs) &&
    identical(pilots[[1]]$truth_days, p$truth_days)
}, logical(1))))
for (p in pilots) {
  validate_sim_dataset(p$complete)
  validate_score_rows(p$scores)
  invisible(lapply(p$fits, validate_fit_result))
}
cfg <- pilots[[3]]$cfg
out <- "results/technical-report"
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
sets <- c("default", "broader", "anchored_default", "initial_rate")
labels <- c("Day 90\ndefault", "Day 90\nbroader",
             "Two anchors\ndefault", "Initial rate\ndefault")
# Keep four rows legible; the broader two-anchor pilot remains saved.
pilots <- lapply(pilots, function(p) {
  p$fits <- Filter(function(f) f$prior_set %in% sets, p$fits)
  p$scores <- p$scores[p$scores$prior_set %in% sets, ]
  p$future <- p$future[p$future$prior_set %in% sets, ]
  keep_key <- function(key) {
    vapply(key, function(k) any(vapply(sets, function(set) {
      grepl(paste0("-", set, "-"), k, fixed = TRUE)
    }, logical(1))), logical(1))
  }
  p$attempts <- p$attempts[keep_key(p$attempts$key), ]
  p$failures <- p$failures[keep_key(p$failures$key), ]
  p
})
models <- c("exponential", "biphasic", "powerlaw", "plateau")
model_labels <- c("Exponential", "Biphasic", "Power law", "Plateau")
colours <- setNames(c("#C44335", "#176B9B", "#7B5196", "#2D8263"), models)
decorate <- function(x) {
  x$prior <- factor(x$prior_set, levels = sets, labels = labels)
  if ("model" %in% names(x)) x$model <- factor(x$model, levels = models)
  x
}
theme_set(theme_minimal(base_size = 9) + theme(
  panel.grid.minor = element_blank(),
  strip.text = element_text(colour = "#20303E", size = 9),
  strip.text.y.right = element_text(angle = 0),
  legend.position = "bottom", legend.title = element_blank(),
  plot.margin = margin(5, 8, 5, 5)
))
model_scales <- list(
  scale_colour_manual(values = colours, breaks = models,
                        labels = model_labels, drop = FALSE),
  scale_fill_manual(values = colours, breaks = models,
                      labels = model_labels, drop = FALSE)
)
save_plot <- function(plot, name, height) {
  ggsave(file.path(out, "figures", paste0(name, ".pdf")), plot,
           width = 7.2, height = height, device = grDevices::pdf,
           useDingbats = FALSE)
  ggsave(file.path(out, "figures", paste0(name, ".png")), plot,
           width = 7.2, height = height, dpi = 160)
}
times <- sort(unique(c(seq(0, 2200, 15), 90, 365)))
summarise_curve <- function(pars, model, set, day = NA_real_) {
  mu <- early_log_curve(pars, model, times)
  q <- apply(mu, 2, quantile, probs = c(.05, .5, .95))
  data.frame(time_days = times, lo = pmax(exp(q[1, ]), 1e-300),
               med = pmax(exp(q[2, ]), 1e-300), hi = exp(q[3, ]),
               model = model, prior_set = set, followup_days = day)
}

prior_curves <- list()
prior_times <- list()
for (set in sets) {
  seed <- derive_seed(paste("prior-predictive", set), cfg$pilot$replicate,
                        cfg$run$seed)
  for (model in models) {
    free <- draw_early_prior(model, cfg, set, 4000, seed)
    pars <- early_transform(free, model)
    key <- paste(set, model)
    prior_curves[[key]] <- summarise_curve(pars, model, set)
    for (endpoint in names(cfg$calibration$thresholds)) {
      lt <- early_log_tstar(pars, model,
                             cfg$calibration$thresholds[[endpoint]])
      prior_times[[paste(key, endpoint)]] <- data.frame(
        prior_set = set, model = model, endpoint = endpoint,
        prob_no_cross = mean(lt == Inf),
        lo_days = exp(quantile(lt, .05, type = 1)),
        median_days = exp(quantile(lt, .5, type = 1)),
        hi_days = exp(quantile(lt, .95, type = 1)),
        prob_numeric_overflow = mean(is.finite(lt) &
                                      lt > log(.Machine$double.xmax))
      )
    }
  }
}
prior_curves <- decorate(bind_rows(prior_curves))
prior_times <- bind_rows(prior_times)
thresholds <- unlist(cfg$calibration$thresholds)
prior_plot <- ggplot(prior_curves, aes(time_days, med, colour = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = .18,
                 colour = NA) + geom_line(linewidth = .5) +
  geom_hline(yintercept = thresholds, linetype = "dashed",
                 colour = "grey45", linewidth = .3) +
  geom_vline(xintercept = c(90, 365), linetype = "dotted",
                 colour = "grey60", linewidth = .3) +
  facet_grid(prior ~ model,
    labeller = labeller(model = setNames(model_labels, models))) +
  scale_y_log10() + model_scales + coord_cartesian(ylim = c(.001, 5)) +
  scale_x_continuous(breaks = c(0, 1000, 2000)) +
  labs(x = "Days since initial measurement",
         y = "Convalescent-relative latent titre") +
  theme(legend.position = "none")
save_plot(prior_plot, "prior-trajectories", 7.0)

scores <- decorate(bind_rows(lapply(pilots, `[[`, "scores")))
scores$endpoint <- factor(scores$endpoint,
                            levels = c("symptomatic", "severe"),
                            labels = c("Symptomatic threshold",
                                         "Severe threshold"))
scores$review_months <- cfg$design$followup_months[
  match(scores$followup_days, cfg$design$followup_days)
]
offset <- setNames(c(-.6, -.2, .2, .6), models)
scores$x <- scores$review_months + offset[as.character(scores$model)]
cap <- 7305
scores$lo <- pmax(1, pmin(scores$tstar_lo, cap))
scores$hi <- pmin(scores$tstar_hi, cap)
scores$med <- pmin(scores$tstar_med, cap)
scores$upper <- ifelse(is.infinite(scores$tstar_hi), "Inf",
                         ifelse(scores$tstar_hi > cap, ">", ""))
truth <- data.frame(endpoint = levels(scores$endpoint),
                       value = as.numeric(pilots[[1]]$truth_days))
truth$endpoint <- factor(truth$endpoint, levels = levels(scores$endpoint))
planned <- expand.grid(prior_set = sets, model = models,
                        followup_days = cfg$design$followup_days,
                        stringsAsFactors = FALSE)
failed <- anti_join(planned, scores,
                     by = c("prior_set", "model", "followup_days"))
failed <- decorate(merge(failed,
                          data.frame(endpoint = levels(scores$endpoint))))
failed$endpoint <- factor(failed$endpoint, levels = levels(scores$endpoint))
failed$x <- cfg$design$followup_months[
  match(failed$followup_days, cfg$design$followup_days)
] + offset[as.character(failed$model)]
duration_plot <- ggplot(scores, aes(x, med, colour = model)) +
  geom_hline(data = truth, aes(yintercept = value), inherit.aes = FALSE,
                 linetype = "dashed", colour = "grey40", linewidth = .35) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .45) +
  geom_point(data = scores[is.finite(scores$tstar_med), ], size = 1.4) +
  geom_point(data = scores[is.infinite(scores$tstar_med), ],
                 shape = 2, size = 1.7) +
  geom_point(data = failed, aes(x, 80), inherit.aes = FALSE,
                 shape = 4, size = 1.8, colour = "grey30") +
  geom_text(aes(y = cap, label = upper), vjust = -.7, size = 2.1,
                show.legend = FALSE) +
  facet_grid(prior ~ endpoint) + model_scales[1] + scale_y_log10() +
  scale_x_continuous(breaks = cfg$design$followup_months) +
  coord_cartesian(ylim = c(70, cap*1.45), clip = "off") +
  labs(x = "Available follow-up (months)",
         y = "Individual crossing time (days; logarithmic scale)")
save_plot(duration_plot, "posterior-durations", 8.0)

reviews <- cfg$design$followup_days[c(1, 2, 4, 6)]
fits <- unlist(lapply(pilots, `[[`, "fits"), recursive = FALSE)
posterior <- bind_rows(lapply(names(fits), function(key) {
  fit <- fits[[key]]
  day <- as.numeric(sub(".*-([0-9]+)d$", "\\1", key))
  if (!day %in% reviews) return(NULL)
  summarise_curve(fit$native_draws, fit$model, fit$prior_set, day)
})) |> decorate()
review_labels <- paste(cfg$design$followup_months[c(1, 2, 4, 6)], "months")
posterior$review <- factor(posterior$followup_days, levels = reviews,
                              labels = review_labels)
b <- cfg$calibration$biphasic
truth_curve <- data.frame(time_days = times,
  value = curve_biphasic(times, b$c0, b$h1_days, b$h2_days, b$ts_days))
obs <- bind_rows(lapply(seq_along(reviews), function(j) {
  x <- pilots[[1]]$complete$obs
  x <- x[x$time_days <= reviews[j], ]
  x$review <- factor(review_labels[j], levels = review_labels)
  x
}))
posterior_plot <- ggplot(posterior, aes(time_days, med, colour = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = .12,
                 colour = NA) + geom_line(linewidth = .5) +
  geom_line(data = truth_curve, aes(time_days, value), inherit.aes = FALSE,
                 colour = "black", linetype = "dashed", linewidth = .5) +
  geom_point(data = obs, aes(time_days, y_obs), inherit.aes = FALSE,
                 colour = "black", size = .7) +
  geom_hline(yintercept = thresholds, linetype = "dotted",
                 colour = "grey45", linewidth = .3) +
  facet_grid(prior ~ review) + model_scales + scale_y_log10() +
  scale_x_continuous(breaks = c(0, 1000, 2000)) +
  coord_cartesian(ylim = c(.005, 2)) +
  labs(x = "Days since initial measurement",
         y = "Convalescent-relative latent titre")
save_plot(posterior_plot, "posterior-trajectories", 8.0)

non_cross <- scores
non_cross$stage <- "Posterior"
prior_cross <- decorate(prior_times)
prior_cross$endpoint <- factor(prior_cross$endpoint,
  levels = c("symptomatic", "severe"), labels = levels(scores$endpoint))
prior_cross$review_months <- 0
prior_cross$stage <- "Prior"
non_cross <- bind_rows(non_cross, prior_cross)
non_cross$x <- non_cross$review_months +
  offset[as.character(non_cross$model)] * .55
cross_plot <- ggplot(non_cross,
                       aes(x, prob_no_cross, colour = model)) +
  geom_line(data = non_cross[non_cross$stage == "Posterior", ],
                linewidth = .4) +
  geom_point(aes(shape = stage), size = 1.4) +
  scale_shape_manual(values = c(Prior = 5, Posterior = 16)) +
  facet_grid(prior ~ endpoint) + model_scales[1] +
  scale_y_continuous(breaks = c(0, .5, 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, cfg$design$followup_months),
                        labels = c("Prior", cfg$design$followup_months)) +
  labs(x = "Available follow-up (months)", y = "Probability of never crossing")
save_plot(cross_plot, "non-crossing", 7.0)

future <- bind_rows(lapply(pilots, `[[`, "future")) |>
  group_by(prior_set, model, review_days) |>
  summarise(mean_log_score = mean(log_predictive_density),
               n_visits = n(), .groups = "drop") |> decorate()
future$review_months <- cfg$design$followup_months[
  match(future$review_days, cfg$design$followup_days)
]
future_plot <- ggplot(future,
  aes(review_months, mean_log_score, colour = model)) +
  geom_line(linewidth = .5) + geom_point(size = 1.6) +
  facet_wrap(~ prior, ncol = 2) + model_scales[1] +
  scale_x_continuous(breaks = cfg$design$followup_months[-6]) +
  labs(x = "Review at which predictions were made (months)",
         y = "Mean log predictive density at subsequent visits")
save_plot(future_plot, "future-prediction", 4.0)

# Numerical tables describe inputs only; results are presented as figures.
fmt <- function(x) format(x, trim = TRUE, scientific = FALSE)
normal <- function(centre, sd) {
  sprintf("$N(%s,\\,%s^2)$", centre, fmt(sd))
}
beta <- function(x) sprintf("$\\mathrm{Beta}(%s,%s)$", fmt(x[1]), fmt(x[2]))
write_rows <- function(rows, name) {
  calibration <- name == "calibration"
  cols <- if (calibration) "Xrl" else "lXX"
  header <- if (calibration) "Quantity & Value & Source / basis" else
    "Free variable & Default & Broader"
  writeLines(c(paste0("\\begin{tabularx}{\\textwidth}{@{}", cols, "@{}}"),
    paste0("\\toprule ", header, " \\\\ \\midrule"),
    paste0(apply(rows, 1, paste, collapse = " & "), " \\\\"),
    "\\bottomrule", "\\end{tabularx}"),
    file.path(out, "tables", paste0(name, ".tex")))
}
ep <- cfg$fitting$priors
ap <- cfg$fitting$priors$two_anchor
calibration <- rbind(
  c("Initial titre $c_0$", fmt(b$c0), "Hogan et al. (2023)"),
  c("Fast half-life $h_1$ (days)", fmt(b$h1_days), "Hogan et al. (2023)"),
  c("Slow half-life $h_2$ (days)", fmt(b$h2_days), "Hogan et al. (2023)"),
  c("Equal-contribution time $t_s$ (days)", fmt(b$ts_days),
       "Hogan et al. (2023)"),
  c("Log-error SD $\\sigma_{\\log}$",
       fmt(cfg$calibration$observation$sigma_log),
       "Workshop choice"),
  c("Symptomatic threshold $c_{\\mathrm{sym}}$", fmt(thresholds[1]),
       "Khoury et al. (2021)"),
  c("Severe threshold $c_{\\mathrm{sev}}$", fmt(thresholds[2]),
       "Khoury et al. (2021)"),
  c("Anchor days $t_a,t_b$", paste(ep$reference_day, ap$later_day, sep = ", "),
       "Workshop choice")
)
write_rows(calibration, "calibration")
common <- rbind(
  c("$\\log c_0$", normal(paste0("\\log ", ep$c0_median),
      ep$default$c0_log_sd), normal(paste0("\\log ", ep$c0_median),
      ep$broader$c0_log_sd)),
  c("$\\log D_a$", normal(paste0("\\log ", ep$decline_median),
      ep$default$decline_log_sd), normal(paste0("\\log ", ep$decline_median),
      ep$broader$decline_log_sd)),
  c("$\\log\\sigma_{\\log}$", normal(paste0("\\log ", ep$sigma_median),
      ep$default$sigma_log_sd), normal(paste0("\\log ", ep$sigma_median),
      ep$broader$sigma_log_sd))
)
write_rows(common, "common-priors")
shapes <- rbind(
  c("Biphasic: $u_f$", beta(unlist(ep$default$slow_fraction_beta)),
       beta(unlist(ep$broader$slow_fraction_beta))),
  c("Biphasic: $z_q$", normal(paste0("\\log ", ep$rate_ratio_median-1),
      ep$default$rate_ratio_log_sd), normal(paste0("\\log ",
      ep$rate_ratio_median-1), ep$broader$rate_ratio_log_sd)),
  c("Power law: $\\log\\beta$", normal(paste0("\\log ", ep$beta_median),
      ep$default$beta_log_sd), normal(paste0("\\log ", ep$beta_median),
      ep$broader$beta_log_sd)),
  c("Plateau: $u_\\infty$", beta(unlist(ep$default$floor_fraction_beta)),
       beta(unlist(ep$broader$floor_fraction_beta)))
)
write_rows(shapes, "early-shape-priors")
anchors <- rbind(
  c("Flexible families: $\\operatorname{logit}\\rho$",
       normal("\\operatorname{logit}0.25", ap$default$rho_logit_sd),
       normal("\\operatorname{logit}0.25", ap$broader$rho_logit_sd)),
  c("Biphasic: $u_s$", beta(unlist(ap$default$slow_rate_fraction_beta)),
       beta(unlist(ap$broader$slow_rate_fraction_beta)))
)
write_rows(anchors, "two-anchor-priors")
rp <- ep$initial_rate
rate_rows <- rbind(
  c("$\\log c_0$", normal(paste0("\\log ", ep$c0_median),
                            ep$default$c0_log_sd)),
  c("$\\log r$ (days$^{-1}$)",
       normal(paste0("\\log ", rp$rate_median_per_day), rp$rate_log_sd)),
  c("$\\log\\sigma_{\\log}$", normal(paste0("\\log ", ep$sigma_median),
                                        ep$default$sigma_log_sd)),
  shapes[1:3, 1:2],
  c("Plateau: $p$", beta(unlist(rp$floor_fraction_beta)))
)
writeLines(c("\\begin{tabularx}{\\textwidth}{@{}lX@{}}",
  "\\toprule Free variable & Prior \\\\ \\midrule",
  paste0(apply(rate_rows, 1, paste, collapse = " & "), " \\\\"),
  "\\bottomrule", "\\end{tabularx}"),
  file.path(out, "tables", "initial-rate-priors.tex"))

# A small macro file keeps reported counts and selected findings tied to
# the saved fits, without duplicating large numerical tables in the report.
macro <- function(name, value) {
  sprintf("\\newcommand{\\%s}{%s}", name, value)
}
macros <- c(macro("TruthSym", sprintf("%.1f", pilots[[1]]$truth_days[1])),
  macro("TruthSev", sprintf("%.1f", pilots[[1]]$truth_days[2])),
  macro("FitCount", length(fits)),
  macro("PlannedCount", nrow(planned)),
  macro("FailureCount", nrow(failed) / length(levels(scores$endpoint))),
  macro("ConvergedCount", sum(vapply(fits, function(f) {
    f$diagnostics$converged
  }, logical(1)))),
  macro("AttemptsCount", sum(vapply(pilots, function(p) {
    nrow(p$attempts) + if (is.null(p$failures)) 0L else nrow(p$failures)
  }, 1L))),
  macro("ReliableLooCount", sum(vapply(fits, function(f) f$loo_reliable,
                                        logical(1)))),
  macro("ObservationCount", nrow(pilots[[1]]$complete$obs)))
for (set in sets) {
  tag <- c(default = "EarlyDefault", broader = "EarlyBroader",
             anchored_default = "AnchorDefault",
             initial_rate = "RateDefault")[[set]]
  for (model in c("exponential", "biphasic", "powerlaw")) {
    z <- scores[scores$prior_set == set &
                  scores$model == model & scores$followup_days == 730 &
                  scores$endpoint == "Severe threshold", ]
    value <- sprintf("%.0f [%.0f, %.0f]", z$tstar_med, z$tstar_lo, z$tstar_hi)
    macros <- c(macros, macro(paste0(tag, switch(model,
      exponential = "Exp", biphasic = "Bip", powerlaw = "Power")), value))
  }
}
rate_plateau <- prior_times[prior_times$prior_set == "initial_rate" &
                             prior_times$model == "plateau", ]
for (endpoint in c("symptomatic", "severe")) {
  macros <- c(macros, macro(paste0("RatePriorNoCross",
    if (endpoint == "symptomatic") "Sym" else "Sev"), sprintf("%.2f",
      rate_plateau$prob_no_cross[rate_plateau$endpoint == endpoint])))
}
rate_fits <- Filter(function(f) f$prior_set == "initial_rate", fits)
macros <- c(macros, macro("RateFitCount", length(rate_fits)),
  macro("RateConvergedCount", sum(vapply(rate_fits, function(f) {
    f$diagnostics$converged
  }, logical(1)))))
writeLines(macros, file.path(out, "tables", "results-macros.tex"))
write.csv(prior_times, file.path(out, "prior-duration-checks.csv"),
              row.names = FALSE)
write.csv(scores, file.path(out, "duration-summary.csv"), row.names = FALSE)
write.csv(future, file.path(out, "future-summary.csv"), row.names = FALSE)
write.csv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))),
              file.path(out, "input-manifest.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
cat("Wrote technical-report figures, input tables and provenance.\n")
