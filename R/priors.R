# Priors and deterministic transformations from handout Section 5.
# Free variables are independent; derived trajectory parameters are not.

early_prior <- function(cfg, prior_set = "default") {
  prior_set <- match.arg(prior_set, c("default", "broader",
    "anchored_default", "anchored_broader", "initial_rate"))
  initial_rate <- prior_set == "initial_rate"
  anchored <- startsWith(prior_set, "anchored_")
  prior_set <- sub("^anchored_", "", prior_set)
  if (initial_rate) prior_set <- "default"
  p <- cfg$fitting$priors
  if (is.null(p[[prior_set]])) stop("Missing fitting prior configuration")
  out <- c(p[c("reference_day", "c0_median", "decline_median", "sigma_median",
                "rate_ratio_median", "beta_median")], p[[prior_set]])
  out$slow_fraction_beta <- unlist(out$slow_fraction_beta)
  out$floor_fraction_beta <- unlist(out$floor_fraction_beta)
  out$two_anchor <- anchored
  out$initial_rate <- initial_rate
  if (initial_rate) {
    out[names(p$initial_rate)] <- p$initial_rate
    out$floor_fraction_beta <- unlist(out$floor_fraction_beta)
  }
  if (anchored) {
    a <- p$two_anchor
    out <- c(out, a[c("later_day", "rho_median")], a[[prior_set]])
    out$slow_rate_fraction_beta <- unlist(out$slow_rate_fraction_beta)
  }
  out
}

early_free_names <- function(model, prior_set = "default") {
  if (prior_set == "initial_rate") {
    return(c("log_c0", "log_rate", "log_sigma", switch(
      match.arg(model, CANDIDATE_MODELS),
      exponential = character(), biphasic = c("u_f", "z_q"),
      powerlaw = "log_beta", plateau = "u_inf"
    )))
  }
  if (startsWith(prior_set, "anchored_") && model != "exponential") {
    return(c("log_c0", "log_D", "log_sigma", "z_rho",
             if (model == "biphasic") "u_s"))
  }
  c("log_c0", "log_D", "log_sigma", switch(
    match.arg(model, CANDIDATE_MODELS),
    exponential = character(), biphasic = c("u_f", "z_q"),
    powerlaw = "log_beta", plateau = "u_inf"
  ))
}

draw_early_prior <- function(model, cfg, prior_set, n, seed) {
  p <- early_prior(cfg, prior_set)
  .with_seed(seed, {
    x <- data.frame(
      log_c0 = rnorm(n, log(p$c0_median), p$c0_log_sd),
      log_D = rnorm(n, log(p$decline_median), p$decline_log_sd),
      log_sigma = rnorm(n, log(p$sigma_median), p$sigma_log_sd)
    )
    if (p$initial_rate) {
      # Use the second common RNG stream for the shared rate, not a decline.
      x$log_rate <- log(p$rate_median_per_day) +
        (x$log_D - log(p$decline_median)) *
        p$rate_log_sd / p$decline_log_sd
      x$log_D <- NULL
    }
    if (p$two_anchor && model != "exponential") {
      x$z_rho <- rnorm(n, qlogis(p$rho_median), p$rho_logit_sd)
      if (model == "biphasic") {
        x$u_s <- stats::rbeta(n, p$slow_rate_fraction_beta[1],
                              p$slow_rate_fraction_beta[2])
      }
      return(x)
    }
    if (model == "biphasic") {
      x$u_f <- stats::rbeta(n, p$slow_fraction_beta[1],
                            p$slow_fraction_beta[2])
      x$z_q <- rnorm(n, log(p$rate_ratio_median - 1),
                     p$rate_ratio_log_sd)
    }
    if (model == "powerlaw") {
      x$log_beta <- rnorm(n, log(p$beta_median), p$beta_log_sd)
    }
    if (model == "plateau") {
      x$u_inf <- stats::rbeta(n, p$floor_fraction_beta[1],
                       p$floor_fraction_beta[2])
    }
    x
  })
}

.early_logadd <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

.early_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# Solve in dimensionless time. Vector bisection avoids a horizon cutoff and
# remains accurate when the slow rate is far below the fast one.
early_transform <- function(draws, model, reference_day = 90,
                            later_day = 365) {
  d <- as.data.frame(draws)
  if ("log_rate" %in% names(d)) {
    return(.initial_rate_transform(d, model, reference_day))
  }
  if ("z_rho" %in% names(d) && model != "exponential") {
    return(.two_anchor_transform(d, model, reference_day, later_day))
  }
  missing <- setdiff(early_free_names(model), names(d))
  if (length(missing)) stop("Missing free parameters: ", paste(missing,
                                                             collapse = ", "))
  x <- data.frame(c0 = exp(d$log_c0), D = exp(d$log_D),
                  sigma = exp(d$log_sigma))
  if (model == "exponential") x$k <- x$D / reference_day
  if (model == "biphasic") {
    x$f <- d$u_f / 2
    x$q <- 1 + exp(d$z_q)
    lo <- x$D / x$q
    hi <- x$D
    for (i in seq_len(60)) {
      mid <- (lo + hi) / 2
      above <- .early_logadd(log1p(-x$f) - x$q * mid,
                             log(x$f) - mid) > -x$D
      lo <- ifelse(above, mid, lo)
      hi <- ifelse(above, hi, mid)
    }
    x$k2 <- (lo + hi) / (2 * reference_day)
    x$k1 <- x$q * x$k2
    x$h1 <- log(2) / x$k1
    x$h2 <- log(2) / x$k2
    x$ts <- (log1p(-x$f) - log(x$f)) / (x$k1 - x$k2)
  }
  if (model == "powerlaw") {
    x$beta <- exp(d$log_beta)
    v <- x$D / x$beta
    x$log_alpha <- v + log(-expm1(-v)) - log(reference_day)
    x$alpha <- exp(x$log_alpha)
  }
  if (model == "plateau") {
    x$log_r <- log(d$u_inf) - x$D
    x$A_inf <- exp(d$log_c0 + x$log_r)
    x$k <- (x$D - log1p(-d$u_inf) + log1p(-exp(x$log_r))) /
      reference_day
  }
  if (anyNA(x) || any(vapply(x, function(z) any(!is.finite(z)), logical(1)))) {
    stop("Non-finite derived parameter; do not silently discard this draw")
  }
  x
}

.initial_rate_transform <- function(d, model, reference_day) {
  missing <- setdiff(early_free_names(model, "initial_rate"), names(d))
  if (length(missing)) stop("Missing initial-rate free parameters")
  x <- data.frame(c0 = exp(d$log_c0), rate = exp(d$log_rate),
                  sigma = exp(d$log_sigma))
  if (model == "exponential") x$k <- x$rate
  if (model == "biphasic") {
    x$f <- d$u_f / 2
    x$q <- 1 + exp(d$z_q)
    x$k2 <- x$rate / ((1 - x$f) * x$q + x$f)
    x$k1 <- x$q * x$k2
    x$h1 <- log(2) / x$k1
    x$h2 <- log(2) / x$k2
    x$ts <- (log1p(-x$f) - log(x$f)) / (x$k1 - x$k2)
  }
  if (model == "powerlaw") {
    x$beta <- exp(d$log_beta)
    x$log_alpha <- d$log_rate - d$log_beta
    x$alpha <- exp(x$log_alpha)
  }
  if (model == "plateau") {
    x$log_r <- log(d$u_inf)
    x$A_inf <- x$c0 * d$u_inf
    x$k <- x$rate / (1 - d$u_inf)
  }
  x$D <- as.numeric(log(x$c0) -
                       early_log_curve(x, model, reference_day))
  if (anyNA(x) || any(vapply(x, function(z) any(!is.finite(z)), logical(1)))) {
    stop("Non-finite initial-rate parameter; no draws may be discarded")
  }
  x
}

# Matching the anchors first avoids conditioning their distribution on a
# separately chosen, potentially incompatible mixture fraction.
.two_anchor_transform <- function(d, model, ta, tb) {
  stopifnot(tb > ta, ta > 0)
  missing <- setdiff(early_free_names(model, "anchored_default"), names(d))
  if (length(missing)) stop("Missing two-anchor free parameters")
  x <- data.frame(c0 = exp(d$log_c0), D = exp(d$log_D),
                  sigma = exp(d$log_sigma), rho = plogis(d$z_rho))
  m <- tb / ta
  x$D_b <- x$D * (1 + (m - 1) * x$rho)
  if (model == "powerlaw") {
    lo <- rep(0, nrow(x))
    hi <- log(m) / ((m - 1) * x$rho)
    for (i in seq_len(60)) {
      v <- (lo + hi) / 2
      above <- log1p((m - 1) * (-expm1(-v))) >
        (m - 1) * x$rho * v
      lo <- ifelse(above, v, lo)
      hi <- ifelse(above, hi, v)
    }
    v <- (lo + hi) / 2
    x$beta <- x$D / v
    x$log_alpha <- v + log(-expm1(-v)) - log(ta)
    x$alpha <- exp(x$log_alpha)
  } else {
    v <- if (model == "biphasic") d$u_s * x$rho * x$D else 0
    lo <- rep(0, nrow(x))
    hi <- rep(1, nrow(x))
    for (i in seq_len(60)) {
      w <- (lo + hi) / 2
      log_den <- log1p(-exp(-x$D + v) * w)
      log_f <- -x$D + v + log1p(-w) - log_den
      log_fast <- log(-expm1(-x$D + v)) - log_den
      log_rb <- .early_logadd(log_fast + m * (-x$D + log(w)),
                               log_f - m * v)
      above <- log_rb > -x$D_b
      lo <- ifelse(above, w, lo)
      hi <- ifelse(above, hi, w)
    }
    w <- (lo + hi) / 2
    log_f <- -x$D + v + log1p(-w) -
      log1p(-exp(-x$D + v) * w)
    k <- (x$D - log(w)) / ta
    if (model == "biphasic") {
      x$f <- exp(log_f)
      x$k1 <- k
      x$k2 <- v / ta
      x$q <- x$k1 / x$k2
      x$h1 <- log(2) / x$k1
      x$h2 <- log(2) / x$k2
      x$ts <- (log1p(-x$f) - log_f) / (x$k1 - x$k2)
    } else {
      x$log_r <- log_f
      x$A_inf <- exp(d$log_c0 + log_f)
      x$k <- k
    }
  }
  if (anyNA(x) || any(vapply(x, function(z) any(!is.finite(z)), logical(1)))) {
    stop("Non-finite two-anchor parameter; no draws may be discarded")
  }
  x
}

early_log_curve <- function(params, model, time_days) {
  p <- params
  n <- nrow(p)
  base <- matrix(log(p$c0), n, length(time_days))
  switch(model,
    exponential = base - outer(p$k, time_days),
    biphasic = base + .early_logadd(
      matrix(log1p(-p$f), n, length(time_days)) - outer(p$k1, time_days),
      matrix(log(p$f), n, length(time_days)) - outer(p$k2, time_days)
    ),
    powerlaw = {
      z <- outer(p$log_alpha, log(pmax(time_days, .Machine$double.xmin)),
                 "+")
      base - p$beta * .early_softplus(z)
    },
    plateau = base + .early_logadd(
      matrix(p$log_r, n, length(time_days)),
      matrix(log1p(-exp(p$log_r)), n, length(time_days)) -
        outer(p$k, time_days)
    )
  )
}

# All positive-rate biphasic draws eventually cross. Solve on a scaled
# interval with an analytic upper bound instead of declaring a long tail Inf.
early_log_tstar <- function(params, model, threshold) {
  if (model != "powerlaw") return(log(early_tstar(params, model, threshold)))
  L <- log(params$c0 / threshold)
  out <- rep(-Inf, nrow(params))
  idx <- L > 0
  z <- L[idx] / params$beta[idx]
  out[idx] <- z + log(-expm1(-z)) - params$log_alpha[idx]
  out
}

early_tstar <- function(params, model, threshold) {
  stopifnot(length(threshold) == 1L, is.finite(threshold), threshold > 0)
  p <- params
  L <- log(p$c0 / threshold)
  if (model == "exponential") out <- pmax(L, 0) / p$k
  if (model == "powerlaw") {
    log_t <- early_log_tstar(p, model, threshold)
    if (any(log_t > log(.Machine$double.xmax))) {
      stop("Finite power-law duration exceeds numeric range")
    }
    out <- exp(log_t)
  }
  if (model == "plateau") {
    out <- rep(Inf, nrow(p))
    out[L <= 0] <- 0
    idx <- L > 0 & p$A_inf < threshold
    out[idx] <- (log(p$c0[idx] - p$A_inf[idx]) -
                   log(threshold - p$A_inf[idx])) / p$k[idx]
  }
  if (model == "biphasic") {
    L <- pmax(L, 0)
    lo <- rep(0, length(L))
    hi <- rep(1, length(L))
    for (i in seq_len(60)) {
      mid <- (lo + hi) / 2
      log_ratio <- .early_logadd(log1p(-p$f) - p$q * L * mid,
                                 log(p$f) - L * mid)
      above <- log_ratio > -L
      lo <- ifelse(above, mid, lo)
      hi <- ifelse(above, hi, mid)
    }
    out <- L / p$k2 * (lo + hi) / 2
  }
  if (anyNA(out)) stop("Undefined duration; no draws may be discarded")
  as.numeric(out)
}

# The transform is finite even when the scientific estimand is infinite.
early_duration_transform <- function(x) {
  out <- x / (365.25 + x)
  out[is.infinite(x)] <- 1
  out
}
