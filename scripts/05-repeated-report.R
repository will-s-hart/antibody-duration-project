#!/usr/bin/env Rscript
# Make figures and numerical macros after the repeated study has completed.
suppressMessages(pkgload::load_all(quiet = TRUE))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  x <- args[startsWith(args, paste0("--", name, "="))]
  if (length(x)) sub(paste0("^--", name, "="), "", x[1]) else default
}
study_path <- option("input", "results/repeated-default/study.rds")
study <- readRDS(study_path)
cfg <- study$cfg
expected <- cfg$repeated$n_replicates * length(cfg$repeated$models) *
  length(cfg$design$followup_days)
stopifnot(nrow(study$status) == expected)
validate_score_rows(study$scores)
out <- option("output", "results/technical-report")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
s <- study$summary
models <- c("exponential", "biphasic", "powerlaw", "plateau")
model_labels <- c("Exponential", "Biphasic", "Power law", "Plateau")
colours <- setNames(c("#C44335", "#176B9B", "#7B5196", "#2D8263"), models)
s$model <- factor(s$model, levels = models)
s$endpoint <- factor(s$endpoint, levels = c("symptomatic", "severe"),
                        labels = c("Symptomatic threshold", "Severe threshold"))
s$months <- cfg$design$followup_months[
  match(s$followup_days, cfg$design$followup_days)
]
offset <- setNames(c(-.36, -.12, .12, .36), models)
s$x <- s$months + offset[as.character(s$model)]
theme_set(theme_minimal(base_size = 10) + theme(
  panel.grid.minor = element_blank(), legend.position = "bottom",
  legend.title = element_blank(), plot.margin = margin(7, 8, 7, 5)
))
common <- list(
  scale_colour_manual(values = colours, breaks = models,
                        labels = model_labels, drop = FALSE),
  scale_x_continuous(breaks = cfg$design$followup_months),
  facet_wrap(~ endpoint, nrow = 1, scales = "free_y"),
  labs(x = "Available follow-up (months)")
)
save_plot <- function(plot, name, height = 3.3) {
  ggsave(file.path(out, "figures", paste0(name, ".pdf")), plot,
           width = 7.2, height = height, device = grDevices::pdf,
           useDingbats = FALSE)
  ggsave(file.path(out, "figures", paste0(name, ".png")), plot,
           width = 7.2, height = height, dpi = 160)
}
coverage <- ggplot(s, aes(x, coverage, colour = model)) +
  geom_hline(yintercept = .9, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(ymin = coverage_all_lower, ymax = coverage_all_upper),
                   linewidth = 3, alpha = .16) +
  geom_errorbar(aes(ymin = coverage_lo, ymax = coverage_hi), width = .14,
                  linewidth = .35) +
  geom_line(linewidth = .5) + geom_point(size = 1.5) + common +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, .9, 1)) +
  labs(y = "Coverage of central 90% interval")
save_plot(coverage, "repeated-coverage", 3.7)

width <- s
cap <- max(unlist(cfg$repeated$sensitivity_horizons_days))
width$value <- pmin(width$median_width_days, cap)
width$label <- ifelse(is.infinite(width$median_width_days), "Inf",
                         ifelse(width$median_width_days > cap, ">", ""))
width_plot <- ggplot(width, aes(x, value, colour = model)) +
  geom_line(linewidth = .5) + geom_point(size = 1.5) +
  geom_text(aes(label = label), vjust = -.8, size = 2.5,
                show.legend = FALSE) + common +
  scale_y_log10() + coord_cartesian(ylim = c(10, cap*1.5), clip = "off") +
  labs(y = "Median 90% interval width (days)")
save_plot(width_plot, "repeated-width")
unbounded <- ggplot(s, aes(x, prob_unbounded, colour = model)) +
  geom_errorbar(aes(ymin = unbounded_lo, ymax = unbounded_hi), width = .14,
                  linewidth = .35) +
  geom_line(linewidth = .5) + geom_point(size = 1.5) + common +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
  labs(y = "Fraction of unbounded 90% intervals")
save_plot(unbounded, "repeated-unbounded")

mean_plot <- function(column, ylabel, cap = NULL, lower_zero = FALSE) {
  d <- s
  d$value <- d[[column]]
  se <- d[[paste0(column, "_mcse")]]
  d$lo <- ifelse(is.finite(se), d$value-1.96*se, d$value)
  d$hi <- ifelse(is.finite(se), d$value+1.96*se, d$value)
  if (lower_zero) d$lo <- pmax(d$lo, 0)
  d$label <- ""
  if (!is.null(cap)) {
    d$label <- ifelse(is.infinite(d$value), "Inf",
                         ifelse(d$value > cap, ">", ""))
    d$value <- pmax(-cap, pmin(d$value, cap))
    d$lo <- pmax(-cap, pmin(d$lo, cap))
    d$hi <- pmax(-cap, pmin(d$hi, cap))
  }
  p <- ggplot(d, aes(x, value, colour = model)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .14,
                    linewidth = .35) +
    geom_line(linewidth = .5) + geom_point(size = 1.5) + common +
    labs(y = ylabel)
  if (!lower_zero) p <- p + geom_hline(yintercept = 0, linetype = "dashed",
                                       colour = "grey45", linewidth = .3)
  if (!is.null(cap)) {
    p <- p + geom_text(aes(label = label), vjust = -.8, size = 2.5,
                         show.legend = FALSE) +
      coord_cartesian(ylim = c(-cap, cap*1.3), clip = "off")
  }
  p
}
save_plot(mean_plot("bias_days", "Mean posterior-median error (days)",
                      cap = cfg$repeated$crps_horizon_days), "repeated-bias")
save_plot(mean_plot("bias_bounded", "Mean posterior-median error on g scale"),
             "repeated-bias-bounded")
save_plot(mean_plot("crps_bounded", "Mean CRPS on g scale", lower_zero = TRUE),
             "repeated-crps")
save_plot(mean_plot("crps_restricted_days",
                      "Mean 10-year restricted CRPS (days)", lower_zero = TRUE),
             "repeated-crps-restricted")

horizons <- c(cfg$repeated$crps_horizon_days,
                unlist(cfg$repeated$sensitivity_horizons_days))
sensitivity <- bind_rows(lapply(seq_along(horizons), function(j) {
  d <- s
  name <- c("crps_restricted_days", "crps_short_days", "crps_long_days")[j]
  d$score <- d[[name]]
  d$horizon <- paste(format(horizons[j]/365.25, trim = TRUE), "years")
  d
}))
sensitivity$horizon <- factor(sensitivity$horizon,
                                 levels = c("5 years", "10 years", "20 years"))
hp <- ggplot(sensitivity, aes(x, score, colour = model)) +
  geom_line(linewidth = .5) + geom_point(size = 1.4) +
  scale_colour_manual(values = colours, breaks = models,
                         labels = model_labels) +
  scale_x_continuous(breaks = cfg$design$followup_months) +
  facet_wrap(~ horizon + endpoint, ncol = 2, scales = "free_y") +
  labs(x = "Available follow-up (months)",
         y = "Mean restricted CRPS (days)")
save_plot(hp, "repeated-horizon", 6.8)

macro <- function(name, value) sprintf("\\newcommand{\\%s}{%s}", name, value)
pct <- function(x) paste0(sprintf("%.0f", 100*x), "\\%")
selected <- function(model, day, endpoint = "Severe threshold") {
  s[s$model == model & s$followup_days == day & s$endpoint == endpoint, ]
}
macros <- c(macro("RepeatN", cfg$repeated$n_replicates),
  macro("RepeatPlanned", expected),
  macro("RepeatCompleted", sum(study$status$completed)),
  macro("RepeatConverged", sum(study$status$converged)),
  macro("RepeatFailed", sum(!study$status$converged)),
  macro("RepeatErrors", sum(!study$status$completed)),
  macro("RepeatUnconverged", sum(study$status$completed &
                                   !study$status$converged)),
  macro("RepeatAttempts", nrow(study$attempts)),
  macro("RepeatCoverageMaxSe", sprintf("%.3f",
    1/(2*sqrt(cfg$repeated$n_replicates)))),
  macro("RepeatCoverageNominalSe", sprintf("%.3f",
    sqrt(.9*.1/cfg$repeated$n_replicates))),
  macro("RepeatMinDenominator", min(s$n_ok)),
  macro("RepeatMaxDenominator", max(s$n_ok)))
for (model in c("exponential", "biphasic", "powerlaw", "plateau")) {
  tag <- c(exponential = "Exp", biphasic = "Bip", powerlaw = "Power",
             plateau = "Plateau")[[model]]
  for (day in c(91, 730)) {
    d <- selected(model, day)
    suffix <- if (day == 91) "Early" else "Late"
    macros <- c(macros,
      macro(paste0("Repeat", tag, suffix, "Coverage"), pct(d$coverage)),
      macro(paste0("Repeat", tag, suffix, "Bias"),
              if (is.infinite(d$bias_days)) "$\\infty$" else
                sprintf("%.0f", d$bias_days)),
      macro(paste0("Repeat", tag, suffix, "Crps"),
              sprintf("%.3f", d$crps_bounded)))
  }
}
writeLines(macros, file.path(out, "tables", "repeated-macros.tex"))
write.csv(data.frame(path = study_path,
  md5 = unname(tools::md5sum(study_path))),
  file.path(out, "repeated-input-manifest.csv"), row.names = FALSE)
cat("Wrote repeated-study figures and report macros.\n")
