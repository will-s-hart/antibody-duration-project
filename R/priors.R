# Priors and deterministic transformations from handout Section 5.
# Free variables are independent; derived trajectory parameters are not.

early_prior <- function(cfg, prior_set = c("default", "broader")) {
  prior_set <- match.arg(prior_set)
  p <- cfg$fitting$priors
  if (is.null(p[[prior_set]])) stop("Missing fitting prior configuration")
  out <- c(p[c("reference_day", "c0_median", "decline_median", "sigma_median",
                "rate_ratio_median", "beta_median")], p[[prior_set]])
  out$slow_fraction_beta <- unlist(out$slow_fraction_beta)
  out$floor_fraction_beta <- unlist(out$floor_fraction_beta)
  out
}

early_free_names <- function(model) {
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
early_transform <- function(draws, model, reference_day = 90) {
  d <- as.data.frame(draws)
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
early_tstar <- function(params, model, threshold) {
  stopifnot(length(threshold) == 1L, is.finite(threshold), threshold > 0)
  p <- params
  L <- log(p$c0 / threshold)
  if (model == "exponential") out <- pmax(L, 0) / p$k
  if (model == "powerlaw") {
    out <- rep(0, nrow(p))
    idx <- L > 0
    z <- L[idx] / p$beta[idx]
    log_t <- z + log(-expm1(-z)) - p$log_alpha[idx]
    if (any(log_t > log(.Machine$double.xmax))) {
      stop("Finite power-law duration exceeds numeric range")
    }
    out[idx] <- exp(log_t)
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
