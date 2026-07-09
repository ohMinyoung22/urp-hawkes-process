rm(list = ls())
gc()

# ============================================================
# 0. True parameters
# ============================================================

set.seed(123)

T_max <- 1000
K <- 3

mean_mu_true <- c(0.2, 0.1, 0.1)
sigma_mu_true <- c(0.30, 0.30, 0.30)
ell_mu_true <- 30

eta_true <- c(5, 3, 2)
phi_true <- c(2, 5, 3)

alpha_true <- matrix(
  c(
    0.57, 0.00, 0.00,
    0.00, 0.55, 0.26,
    0.00, 0.00, 0.73
  ),
  nrow = 3,
  byrow = TRUE
)

gamma_true <- matrix(
  c(
    0.00, 0.00, 0.26,
    0.00, 0.00, 0.00,
    0.14, 0.00, 0.00
  ),
  nrow = 3,
  byrow = TRUE
)

trpar_true <- list(
  alpha = alpha_true,
  gamma = gamma_true,
  eta = eta_true,
  phi = phi_true
)

cseg <- seq(0, T_max, by = 1)


# ============================================================
# 1. Simulate independent log-GP baseline
# ============================================================

simulate_log_mu_gp <- function(cseg, mean_mu, sigma_mu, ell_mu) {
  
  K <- length(mean_mu)
  G <- length(cseg)
  
  if (length(sigma_mu) == 1) {
    sigma_mu <- rep(sigma_mu, K)
  }
  
  lmu_mean <- log(mean_mu) - 0.5 * sigma_mu^2
  
  lmu <- matrix(NA_real_, nrow = G, ncol = K)
  lmu[1, ] <- lmu_mean + sigma_mu * rnorm(K)
  
  for (r in 2:G) {
    delta_s <- cseg[r] - cseg[r - 1]
    rho <- exp(-delta_s / ell_mu)
    innovation_scale <- sqrt(-expm1(-2 * delta_s / ell_mu))
    
    lmu[r, ] <-
      lmu_mean +
      rho * (lmu[r - 1, ] - lmu_mean) +
      sigma_mu * innovation_scale * rnorm(K)
  }
  
  colnames(lmu) <- paste0("process_", 1:K)
  
  list(
    log_mu = lmu,
    mu = exp(lmu),
    lmu_mean = lmu_mean,
    sigma_mu = sigma_mu,
    ell_mu = ell_mu,
    cseg = cseg
  )
}


# ============================================================
# 2. Simulate MHPWI events
# ============================================================

simulate_MHPWI_GP <- function(log_mu_grid, cseg, trpar, maxT) {
  
  alpha <- trpar$alpha
  gamma <- trpar$gamma
  eta <- trpar$eta
  phi <- trpar$phi
  
  K <- nrow(alpha)
  
  alpha_scaled <- sweep(alpha, 1, eta, "/")
  
  baseline_at <- function(tt) {
    left <- findInterval(tt, cseg, all.inside = TRUE)
    
    if (left >= length(cseg)) {
      return(exp(log_mu_grid[nrow(log_mu_grid), ]))
    }
    
    weight <- (tt - cseg[left]) / (cseg[left + 1] - cseg[left])
    lmu_tt <- (1 - weight) * log_mu_grid[left, ] +
      weight * log_mu_grid[left + 1, ]
    
    exp(lmu_tt)
  }
  
  maxmu <- exp(apply(log_mu_grid, 2, max))
  alpha_row_total <- rowSums(alpha_scaled)
  
  t_now <- 0
  events <- vector("list", K)
  for (k in 1:K) events[[k]] <- numeric(0)
  
  all_time <- numeric(0)
  all_mark <- integer(0)
  
  while (t_now < maxT) {
    
    Rg_now <- numeric(K)
    
    for (l in 1:K) {
      if (length(events[[l]]) > 0) {
        ages <- t_now - events[[l]]
        Rg_now[l] <- sum(exp(-ages / eta[l]))
      }
    }
    
    maxintensity <- sum(maxmu) + sum(alpha_row_total * Rg_now)
    
    if (!is.finite(maxintensity) || maxintensity <= 0) break
    
    t_candidate <- t_now + rexp(1, rate = maxintensity)
    
    if (t_candidate > maxT) break
    
    Rg <- numeric(K)
    Rh <- numeric(K)
    
    for (l in 1:K) {
      if (length(events[[l]]) > 0) {
        ages <- t_candidate - events[[l]]
        Rg[l] <- sum(exp(-ages / eta[l]))
        Rh[l] <- sum(exp(-ages / phi[l]))
      }
    }
    
    mu_tt <- baseline_at(t_candidate)
    
    excitation <- mu_tt + as.vector(crossprod(alpha_scaled, Rg))
    inhibition <- as.vector(crossprod(gamma, Rh))
    
    lambda <- excitation * exp(-inhibition)
    total_lambda <- sum(lambda)
    
    if (total_lambda > maxintensity * (1 + 1e-8)) {
      stop("Thinning upper bound violated.")
    }
    
    if (total_lambda > 0 && runif(1) < total_lambda / maxintensity) {
      k_star <- sample(1:K, size = 1, prob = lambda / total_lambda)
      
      events[[k_star]] <- c(events[[k_star]], t_candidate)
      all_time <- c(all_time, t_candidate)
      all_mark <- c(all_mark, k_star)
      
      if (length(all_time) %% 100 == 0) {
        message(length(all_time), " points generated...")
      }
    }
    
    t_now <- t_candidate
  }
  
  data.frame(
    time = all_time,
    mark = all_mark
  )
}


# ============================================================
# 3. Validation functions
# ============================================================

validate_event_data <- function(events, T_max, K) {
  
  checks <- c(
    no_missing = !anyNA(events),
    finite_time = all(is.finite(events$time)),
    time_in_range = all(events$time > 0 & events$time <= T_max),
    strictly_increasing = nrow(events) <= 1 || all(diff(events$time) > 0),
    valid_mark = all(events$mark %in% 1:K)
  )
  
  print(checks)
  
  cat("\nNumber of events:", nrow(events), "\n")
  
  cat("\nEvents by mark:\n")
  print(table(factor(events$mark, levels = 1:K)))
  
  cat("\nInter-event time summary:\n")
  print(summary(diff(events$time)))
  
  invisible(checks)
}


check_gp_innovations <- function(gp) {
  
  lmu <- gp$log_mu
  cseg <- gp$cseg
  lmu_mean <- gp$lmu_mean
  sigma_mu <- gp$sigma_mu
  ell_mu <- gp$ell_mu
  
  G <- nrow(lmu)
  K <- ncol(lmu)
  
  z <- matrix(NA_real_, nrow = G - 1, ncol = K)
  
  for (r in 1:(G - 1)) {
    delta_s <- cseg[r + 1] - cseg[r]
    rho <- exp(-delta_s / ell_mu)
    innovation_scale <- sqrt(-expm1(-2 * delta_s / ell_mu))
    
    z[r, ] <-
      (
        lmu[r + 1, ] -
          lmu_mean -
          rho * (lmu[r, ] - lmu_mean)
      ) /
      (sigma_mu * innovation_scale)
  }
  
  colnames(z) <- paste0("process_", 1:K)
  
  cat("Innovation means:\n")
  print(colMeans(z))
  
  cat("\nInnovation standard deviations:\n")
  print(apply(z, 2, sd))
  
  cat("\nCross-process innovation correlations:\n")
  print(cor(z))
  
  invisible(z)
}


make_true_intensity_function <- function(true_gp, trpar) {
  
  cseg <- true_gp$cseg
  log_mu_grid <- true_gp$log_mu
  
  alpha <- trpar$alpha
  gamma <- trpar$gamma
  eta <- trpar$eta
  phi <- trpar$phi
  
  K <- nrow(alpha)
  alpha_scaled <- sweep(alpha, 1, eta, "/")
  
  baseline_at <- function(tt) {
    left <- findInterval(tt, cseg, all.inside = TRUE)
    
    if (left >= length(cseg)) {
      return(exp(log_mu_grid[nrow(log_mu_grid), ]))
    }
    
    weight <- (tt - cseg[left]) / (cseg[left + 1] - cseg[left])
    lmu_tt <- (1 - weight) * log_mu_grid[left, ] +
      weight * log_mu_grid[left + 1, ]
    
    exp(lmu_tt)
  }
  
  function(tt, history_time = numeric(0), history_mark = integer(0)) {
    
    Rg <- numeric(K)
    Rh <- numeric(K)
    
    for (l in 1:K) {
      idx <- which(history_mark == l & history_time < tt)
      
      if (length(idx) > 0) {
        ages <- tt - history_time[idx]
        Rg[l] <- sum(exp(-ages / eta[l]))
        Rh[l] <- sum(exp(-ages / phi[l]))
      }
    }
    
    excitation <- baseline_at(tt) + as.vector(crossprod(alpha_scaled, Rg))
    inhibition <- as.vector(crossprod(gamma, Rh))
    
    lambda <- excitation * exp(-inhibition)
    names(lambda) <- paste0("lambda_", 1:K)
    
    lambda
  }
}


time_rescaling_check_GL8 <- function(events, intensity_function, breakpoints) {
  
  x <- c(
    -0.9602898564975363, -0.7966664774136267,
    -0.5255324099163290, -0.1834346424956498,
    0.1834346424956498,  0.5255324099163290,
    0.7966664774136267,  0.9602898564975363
  )
  
  w <- c(
    0.1012285362903763, 0.2223810344533745,
    0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763
  )
  
  N <- nrow(events)
  z <- numeric(N)
  
  history_time <- numeric(0)
  history_mark <- integer(0)
  previous_time <- 0
  
  for (i in 1:N) {
    
    current_time <- events$time[i]
    
    internal_cuts <- breakpoints[
      breakpoints > previous_time &
        breakpoints < current_time
    ]
    
    cuts <- sort(unique(c(previous_time, internal_cuts, current_time)))
    
    interval_integral <- 0
    
    for (j in 1:(length(cuts) - 1)) {
      
      a <- cuts[j]
      b <- cuts[j + 1]
      
      midpoint <- 0.5 * (a + b)
      half_width <- 0.5 * (b - a)
      
      quad_time <- midpoint + half_width * x
      
      total_lambda <- vapply(
        quad_time,
        function(tt) {
          sum(intensity_function(
            tt = tt,
            history_time = history_time,
            history_mark = history_mark
          ))
        },
        numeric(1)
      )
      
      interval_integral <- interval_integral +
        half_width * sum(w * total_lambda)
    }
    
    z[i] <- interval_integral
    
    history_time <- c(history_time, current_time)
    history_mark <- c(history_mark, events$mark[i])
    previous_time <- current_time
    
    if (i %% 100 == 0 || i == N) {
      message("Time-rescaling: ", i, " / ", N)
    }
  }
  
  list(
    integrated_hazard = z,
    mean = mean(z),
    variance = var(z),
    lag1_cor = acf(z, plot = FALSE)$acf[2],
    max_z = max(z),
    prob_one_larger = 1 - pexp(max(z))^length(z)
  )
}


# ============================================================
# 4. Generate data
# ============================================================

true_gp <- simulate_log_mu_gp(
  cseg = cseg,
  mean_mu = mean_mu_true,
  sigma_mu = sigma_mu_true,
  ell_mu = ell_mu_true
)

simulated_events <- simulate_MHPWI_GP(
  log_mu_grid = true_gp$log_mu,
  cseg = cseg,
  trpar = trpar_true,
  maxT = T_max
)


# ============================================================
# 5. Validate generated data
# ============================================================

cat("\n================ Event data check ================\n")
validate_event_data(
  events = simulated_events,
  T_max = T_max,
  K = K
)

cat("\n================ GP innovation check ================\n")
z_gp <- check_gp_innovations(true_gp)

cat("\n================ Time-rescaling check ================\n")

true_intensity <- make_true_intensity_function(
  true_gp = true_gp,
  trpar = trpar_true
)

rescaling_result <- time_rescaling_check_GL8(
  events = simulated_events,
  intensity_function = true_intensity,
  breakpoints = true_gp$cseg
)

print(rescaling_result[c(
  "mean",
  "variance",
  "lag1_cor",
  "max_z",
  "prob_one_larger"
)])

cat("\nLargest integrated hazards:\n")

z <- rescaling_result$integrated_hazard
dt <- diff(c(0, simulated_events$time))
idx <- order(z, decreasing = TRUE)[1:5]

print(data.frame(
  event = idx,
  time = simulated_events$time[idx],
  mark = simulated_events$mark[idx],
  waiting_time = dt[idx],
  integrated_hazard = z[idx]
))


# ============================================================
# 6. Quick decision rule
# ============================================================

cat("\n================ Decision ================\n")

if (
  rescaling_result$mean > 0.9 &&
  rescaling_result$mean < 1.1 &&
  rescaling_result$variance > 0.8 &&
  rescaling_result$variance < 1.3 &&
  abs(rescaling_result$lag1_cor) < 0.08 &&
  rescaling_result$max_z < 10
) {
  cat("PASS: This dataset is safe to use for recovery.\n")
} else {
  cat("CAUTION: Dataset is mostly okay, but consider regenerating with another seed.\n")
}


# ============================================================
# 7. Optional plots
# ============================================================

par(mfrow = c(1, 2))

matplot(
  true_gp$cseg / 3600,
  true_gp$mu,
  type = "l",
  lty = 1,
  xlab = "Time (hours)",
  ylab = expression(mu[k](t)),
  main = "True GP baseline"
)
getwd()
abline(h = mean_mu_true, lty = 2)

qqplot(
  qexp(ppoints(length(z))),
  sort(z),
  xlab = "Theoretical Exp(1)",
  ylab = "Observed integrated hazards",
  main = "Time-rescaling QQ plot"
)

abline(0, 1, lty = 2)

par(mfrow = c(1, 1))

write.csv(
  simulated_events,
  file = "simulated_events_1000.csv",
  row.names = FALSE
)
