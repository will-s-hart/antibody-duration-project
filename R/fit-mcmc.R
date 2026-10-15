# ---------------------------------------------------------------------------
# Alternate models, model selection, and Bayesian inference.
#
# OWNER: workstream (b), simulation and inference. Contract in INTERFACES.md.
#
# Two jobs. First, fit the full set of candidate families from handout
# Section 4 -- exponential, biphasic, power law, plateau -- without assuming
# which is correct, and quantify how much of the uncertainty in the duration
# comes from the choice of family rather than from within any one of them.
# Second, provide the Bayesian fits that give honest intervals where the
# curvature-based ones cannot.
#
# THE ENGINE IS THIS WORKSTREAM'S CHOICE. Stan through cmdstanr, JAGS through
# rjags, brms, or a hand-written sampler are all acceptable; all of them are
# in Suggests rather than Imports so that nobody else has to install a
# toolchain. What is fixed is the output: posterior draws in
# `fit_result$tstar_draws`, and whatever diagnostics the engine provides
# recorded in `$diagnostics`. Keep priors fixed across design cells, so that
# differences between cells reflect the data and not per-cell tuning.
#
# The plateau family needs a decision before it is included in comparisons:
# its duration is undefined whenever the set-point sits above the threshold,
# so a naive model average over families is polluted by infinities. Record
# what you decide and why.
# ---------------------------------------------------------------------------

#' Fit one candidate family by MCMC.
#'
#' @param data a [sim_dataset()].
#' @param model one of [CANDIDATE_MODELS].
#' @param cfg configuration from [read_config()].
#' @param engine identifier recorded in `$provenance$engine`.
#' @return a [fit_result()] with `method = "mcmc"`. Populate `$diagnostics`
#'   with `rhat_max`, `ess_bulk_min` and `divergences`, and set `converged`
#'   from an explicit rule agreed with the group -- filtering after the fact
#'   on an undeclared rule changes what the aggregate means.
fit_mcmc <- function(data, model = "biphasic", cfg, engine = "jags",
                     prior_set = "default", sampler = cfg$fitting$sampler,
                     work_dir = tempfile("antibody-jags-"), ...) {
  validate_sim_dataset(data)
  model <- match.arg(model, CANDIDATE_MODELS)
  if (engine != "jags") stop("The implemented engine is 'jags'")
  if (length(unique(data$obs$participant_id)) != 1L) {
    stop("fit_mcmc(): this pilot fits one participant at a time")
  }
  binary <- Sys.which("jags")
  if (!nzchar(binary)) stop("Install the JAGS command-line sampler")
  for (pkg in c("coda", "posterior", "loo")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg)
  }
  p <- early_prior(cfg, prior_set)
  obs <- data$obs
  if (any(obs$y_obs[obs$censor == "none"] <= 0)) {
    stop("Uncensored log-normal observations must be positive")
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir <- normalizePath(work_dir)
  start <- proc.time()[[3]]
  jdata <- list(
    N = nrow(obs), t = obs$time_days, ta = p$reference_day,
    c0_mean = log(p$c0_median), c0_prec = 1 / p$c0_log_sd^2,
    D_mean = log(p$decline_median), D_prec = 1 / p$decline_log_sd^2,
    sigma_mean = log(p$sigma_median), sigma_prec = 1 / p$sigma_log_sd^2
  )
  if (model == "biphasic") {
    jdata <- c(jdata, list(f_a = p$slow_fraction_beta[1],
                           f_b = p$slow_fraction_beta[2],
                           q_mean = log(p$rate_ratio_median - 1),
                           q_prec = 1 / p$rate_ratio_log_sd^2))
  }
  if (model == "powerlaw") {
    jdata <- c(jdata, list(beta_mean = log(p$beta_median),
                           beta_prec = 1 / p$beta_log_sd^2))
  }
  if (model == "plateau") {
    jdata <- c(jdata, list(u_a = p$floor_fraction_beta[1],
                           u_b = p$floor_fraction_beta[2]))
  }
  blocks <- character()
  for (kind in c("none", "left", "right")) {
    idx <- which(obs$censor == kind)
    if (!length(idx)) next
    label <- if (kind == "none") "unc" else kind
    jdata[[paste0("n_", label)]] <- length(idx)
    jdata[[paste0(label, "_idx")]] <- idx
    if (kind == "none") {
      jdata$log_y <- log(obs$y_obs[idx])
      blocks <- c(blocks, "for (j in 1:n_unc) {
        log_y[j] ~ dnorm(mu[unc_idx[j]], tau)
      }")
    } else {
      limit <- if (kind == "left") data$cell$lloq else data$cell$uloq
      if (is.na(limit)) stop("Missing assay limit for censored observation")
      jdata[[paste0(label, "_one")]] <- rep(1L, length(idx))
      jdata[[paste0("log_", label, "_limit")]] <- log(limit)
      z <- if (kind == "left") {
        "(log_left_limit - mu[left_idx[j]]) / sigma"
      } else {
        "(mu[right_idx[j]] - log_right_limit) / sigma"
      }
      blocks <- c(blocks, sprintf(
        "for (j in 1:n_%s) { %s_one[j] ~ dbern(phi(%s)) }",
        label, label, z
      ))
    }
  }
  model_file <- file.path(work_dir, "model.bug")
  writeLines(.early_jags_model(model, blocks), model_file)
  data_file <- file.path(work_dir, "data.R")
  dump(names(jdata), file = data_file, envir = list2env(jdata))
  commands <- c(sprintf("model in \"%s\"", model_file),
                sprintf("data in \"%s\"", data_file),
                sprintf("compile, nchains(%d)", sampler$chains))
  chain_seeds <- integer(sampler$chains)
  for (chain in seq_len(sampler$chains)) {
    seed <- derive_seed(paste(data$cell$cell_id, model, prior_set, chain),
                        data$cell$replicate, data$cell$seed)
    chain_seeds[chain] <- seed
    init <- as.list(draw_early_prior(model, cfg, prior_set, 1, seed)[1, ])
    init$.RNG.name <- "base::Wichmann-Hill"
    init$.RNG.seed <- max(1L, seed)
    init_file <- file.path(work_dir, sprintf("init%d.R", chain))
    dump(names(init), file = init_file, envir = list2env(init))
    commands <- c(commands, sprintf("parameters in \"%s\", chain(%d)",
                                     init_file, chain))
  }
  stem <- file.path(work_dir, "samples")
  commands <- c(commands, "initialize", sprintf("adapt %d", sampler$adapt),
                sprintf("update %d", sampler$warmup),
                paste("monitor", early_free_names(model)),
                sprintf("update %d", sampler$draws),
                sprintf("coda *, stem(\"%s\")", stem), "exit")
  command_file <- file.path(work_dir, "run.cmd")
  writeLines(commands, command_file)
  log_file <- file.path(work_dir, "jags.log")
  status <- system2(binary, shQuote(command_file), stdout = log_file,
                    stderr = log_file)
  if (status != 0L) {
    stop("JAGS failed: ", paste(utils::tail(readLines(log_file), 12),
                                  collapse = "\n"))
  }
  chains <- lapply(seq_len(sampler$chains), function(i) {
    coda::read.coda(paste0(stem, "chain", i, ".txt"),
                    paste0(stem, "index.txt"), quiet = TRUE)
  })
  chains <- do.call(coda::mcmc.list, chains)
  draws <- as.data.frame(as.matrix(chains))
  draws$chain_id <- rep(seq_len(sampler$chains), each = sampler$draws)
  pars <- early_transform(draws, model, p$reference_day)
  durations <- lapply(cfg$calibration$thresholds, function(thr) {
    early_tstar(pars, model, thr)
  })
  diagnostics <- .early_diagnostics(chains, durations, sampler)
  log_mu <- early_log_curve(pars, model, obs$time_days)
  loglik <- .early_loglik(log_mu, pars$sigma, obs, data$cell)
  # Account for serial dependence in the importance-sampling diagnostics.
  shifted <- sweep(loglik, 2, apply(loglik, 2, max), "-")
  r_eff <- loo::relative_eff(exp(shifted), chain_id = draws$chain_id)
  loo_fit <- suppressWarnings(loo::loo(loglik, r_eff = r_eff, cores = 1))
  pareto_k <- loo::pareto_k_values(loo_fit)
  names_native <- curve_params(model)
  long_params <- data.frame(
    draw_id = rep(seq_len(nrow(pars)), each = length(names_native)),
    param = rep(names_native, times = nrow(pars)),
    value = as.numeric(t(as.matrix(pars[names_native])))
  )
  selected <- names(cfg$calibration$thresholds)[
    vapply(cfg$calibration$thresholds, function(z) {
      isTRUE(all.equal(z, data$cell$c_thr))
    }, logical(1))
  ]
  tt <- if (length(selected)) durations[[selected[1]]] else {
    early_tstar(pars, model, data$cell$c_thr)
  }
  fit <- fit_result(
    data$cell$cell_id, data$cell$replicate, "mcmc", model, long_params,
    data.frame(draw_id = seq_len(nrow(pars)), estimand = "individual",
               participant_id = unique(obs$participant_id),
               tstar_days = tt, crossed = is.finite(tt)),
    diagnostics = default_diagnostics(
      converged = diagnostics$converged,
      rhat_max = max(diagnostics$table$rhat, na.rm = TRUE),
      ess_bulk_min = min(diagnostics$table$ess_bulk, na.rm = TRUE),
      ess_tail_min = min(diagnostics$table$ess_tail, na.rm = TRUE),
      divergences = NA_integer_, n_draws = nrow(pars),
      runtime_s = proc.time()[[3]] - start,
      message = if (diagnostics$converged) "Diagnostic thresholds passed" else
        "Failed diagnostic thresholds; inspect diagnostics_table"
    ),
    ic = default_ic(n_par = length(early_free_names(model)),
                    loo = loo_fit$estimates["looic", "Estimate"]),
    provenance = default_provenance(
      engine = "JAGS CLI", engine_version = readLines(log_file)[1],
      seed = data$cell$seed
    )
  )
  fit$posterior <- draws
  fit$native_draws <- pars
  fit$threshold_draws <- durations
  fit$diagnostics_table <- diagnostics$table
  fit$constant_indicators <- diagnostics$constant_indicators
  fit$log_lik <- loglik
  fit$loo <- loo_fit
  fit$pareto_k_max <- max(pareto_k)
  fit$loo_reliable <- diagnostics$converged && all(pareto_k <= 0.7)
  fit$prior_set <- prior_set
  fit$sampler <- sampler
  fit$chain_seeds <- chain_seeds
  fit$work_dir <- work_dir
  validate_fit_result(fit)
}

.early_jags_model <- function(model, likelihood) {
  shared <- "log_c0 ~ dnorm(c0_mean, c0_prec)
    log_D ~ dnorm(D_mean, D_prec)
    log_sigma ~ dnorm(sigma_mean, sigma_prec)
    D <- exp(log_D)
    sigma <- exp(log_sigma)
    tau <- exp(-2 * log_sigma)"
  shape <- switch(model,
    exponential = "k <- D / ta",
    biphasic = "u_f ~ dbeta(f_a, f_b)
      z_q ~ dnorm(q_mean, q_prec)
      f <- u_f / 2
      q <- 1 + exp(z_q)
      lower[1] <- D / q
      upper[1] <- D
      for (b in 1:50) {
        midpoint[b] <- (lower[b] + upper[b]) / 2
        a[b] <- log(1-f) - q * midpoint[b]
        c[b] <- log(f) - midpoint[b]
        m[b] <- max(a[b], c[b])
        value[b] <- m[b] + log(exp(a[b]-m[b])+exp(c[b]-m[b]))
        above[b] <- step(value[b] + D)
        lower[b+1] <- above[b]*midpoint[b]+(1-above[b])*lower[b]
        upper[b+1] <- above[b]*upper[b]+(1-above[b])*midpoint[b]
      }
      s <- (lower[51] + upper[51]) / (2 * ta)",
    powerlaw = "log_beta ~ dnorm(beta_mean, beta_prec)
      beta <- exp(log_beta)
      v <- D / beta
      log_alpha <- v + log(1-exp(-v)) - log(ta)",
    plateau = "u_inf ~ dbeta(u_a, u_b)
      log_r <- log(u_inf) - D
      k <- (D-log(1-u_inf)+log(1-exp(log_r))) / ta"
  )
  curve <- switch(model,
    exponential = "mu[i] <- log_c0 - k*t[i]",
    biphasic = "fast[i] <- log(1-f)-q*s*t[i]
      slow[i] <- log(f)-s*t[i]
      mix[i] <- max(fast[i],slow[i])
      mu[i] <- log_c0+mix[i]+
        log(exp(fast[i]-mix[i])+exp(slow[i]-mix[i]))",
    powerlaw = "z[i] <- log_alpha + log(max(t[i],1.0E-100))
      mu[i] <- log_c0-beta*(max(z[i],0)+log(1+exp(-abs(z[i]))))",
    plateau = "fall[i] <- log(1-exp(log_r))-k*t[i]
      mix[i] <- max(log_r,fall[i])
      mu[i] <- log_c0+mix[i]+
        log(exp(log_r-mix[i])+exp(fall[i]-mix[i]))"
  )
  paste("model {", shared, shape, "for (i in 1:N) {", curve, "}",
        paste(likelihood, collapse = "\n"), "}", sep = "\n")
}

.early_loglik <- function(log_mu, sigma, obs, cell) {
  ll <- matrix(NA_real_, nrow(log_mu), ncol(log_mu))
  for (j in seq_len(nrow(obs))) {
    ll[, j] <- switch(obs$censor[j],
      none = stats::dnorm(log(obs$y_obs[j]), log_mu[, j], sigma, log = TRUE),
      left = stats::pnorm(log(cell$lloq), log_mu[, j], sigma, log.p = TRUE),
      right = stats::pnorm(log(cell$uloq), log_mu[, j], sigma,
                     lower.tail = FALSE, log.p = TRUE)
    )
  }
  ll
}

.early_diagnostics <- function(chains, durations, sampler) {
  a <- as.array(chains)
  if (length(dim(a)) == 2L) {
    a <- array(a, c(nrow(a), 1L, ncol(a)),
                dimnames = list(NULL, colnames(chains[[1]]), NULL))
  }
  a <- aperm(a, c(1, 3, 2))
  variables <- lapply(seq_len(dim(a)[3]), function(j) a[, , j])
  names(variables) <- dimnames(a)[[3]]
  constant <- character()
  for (nm in names(durations)) {
    d <- durations[[nm]]
    variables[[paste0("duration_", nm)]] <- matrix(
      early_duration_transform(d), nrow = dim(a)[1]
    )
    indicator <- as.numeric(is.infinite(d))
    if (length(unique(indicator)) > 1L) {
      variables[[paste0("no_cross_", nm)]] <- matrix(
        indicator, nrow = dim(a)[1]
      )
    } else {
      constant <- c(constant, nm)
    }
  }
  table <- do.call(rbind, lapply(names(variables), function(nm) {
    v <- variables[[nm]]
    if (length(unique(as.numeric(v))) == 1L) {
      # Only a constant scientific outcome can bypass a convergence test;
      # a frozen free parameter is a sampler failure.
      is_outcome <- startsWith(nm, "duration_")
      return(data.frame(variable = nm, rhat = if (is_outcome) 1 else Inf,
                         ess_bulk = if (is_outcome) NA_real_ else 0,
                         ess_tail = if (is_outcome) NA_real_ else 0,
                         constant = TRUE))
    }
    tail_ess <- vapply(c(0.05, 0.95), function(prob) {
      at_quantile <- v <= quantile(v, prob, names = FALSE)
      if (length(unique(as.numeric(at_quantile))) == 1L) return(NA_real_)
      posterior::ess_quantile(v, probs = prob)
    }, numeric(1))
    data.frame(variable = nm, rhat = posterior::rhat(v),
               ess_bulk = posterior::ess_bulk(v),
               ess_tail = if (all(is.na(tail_ess))) NA_real_ else
                 min(tail_ess, na.rm = TRUE), constant = FALSE)
  }))
  tested <- table[!table$constant, ]
  passed <- all(is.finite(tested$rhat)) &&
    all(tested$rhat < sampler$rhat_max) &&
    all(tested$ess_bulk >= sampler$ess_min) &&
    all(is.na(tested$ess_tail) | tested$ess_tail >= sampler$ess_min) &&
    !any(is.infinite(table$rhat))
  list(converged = passed, table = table, constant_indicators = constant)
}

#' Fit every candidate family to the same dataset.
#'
#' @return a named list of [fit_result()], one per family.
fit_all_models <- function(data, cfg, models = CANDIDATE_MODELS, ...) {
  out <- lapply(models, function(model) fit_mcmc(data, model, cfg, ...))
  stats::setNames(out, models)
}

#' Compare candidate families on one dataset.
#'
#' Return the information criteria per family and the selected family. Whether
#' selection uses AIC, BIC, WAIC or LOO is this workstream's call; record it in
#' the returned table so a figure can state it.
compare_models <- function(fits, cfg) {
  rows <- lapply(fits, function(fit) data.frame(
    model = fit$model, prior_set = fit$prior_set,
    elpd_loo = fit$loo$estimates["elpd_loo", "Estimate"],
    elpd_loo_se = fit$loo$estimates["elpd_loo", "SE"],
    p_loo = fit$loo$estimates["p_loo", "Estimate"],
    pareto_k_max = fit$pareto_k_max,
    converged = fit$diagnostics$converged,
    loo_reliable = fit$loo_reliable
  ))
  do.call(rbind, rows)
}

#' Pool duration draws across families by their weights.
#'
#' Model averaging is not automatically an improvement: it can be worse than
#' the best single family, and it is sensitive to how non-crossing draws are
#' treated. Report it beside the individual families rather than in place of
#' them.
#'
#' @return a [fit_result()] whose `$tstar_draws` are the pooled draws.
model_average_tstar <- function(fits, cfg, weights = "aic") {
  .not_implemented("model_average_tstar",
                   "(b) simulation, inference and MCMC")
}
