rm(list = ls()); gc()
# =============================================================================-
library(cmdstanr)
library(bayesplot)
library(tidyverse)

# =============================================================================-
# model ----
# =============================================================================-

MHPWI_model_code = "
functions {
  
  // ===========================================================================
  // lambda_total(u) Function
  // ===========================================================================
  //mu_u는 각 K에서의 절대시간에서의 mu
  
  real lambda_total_u(
    real u,
    vector R_g,
    vector R_h,
    matrix alpha_scaled,
    vector eta,
    matrix gamma,
    vector phi,
    vector mu_u,
    int K
  ){
    vector[K] Rg_u;
    vector[K] Rh_u;
  
    vector[K] excitation;
    vector[K] inhibition;
  
    for (l in 1:K) {
      Rg_u[l] = R_g[l] * exp(-u / eta[l]);
      Rh_u[l] = R_h[l] * exp(-u / phi[l]);
    }
  
    excitation = mu_u + alpha_scaled' * Rg_u;
    inhibition = gamma' * Rh_u;
  
    return dot_product(excitation, exp(-inhibition));
  }
  
  // ===========================================================================
  // Log-Likelihood Function
  // ===========================================================================
  real MHPWI_loglike_lpdf(
    vector t,
    array[] int m,
    int K,
    real T_end,
    matrix alpha,
    vector eta,
    matrix gamma,
    vector phi,
    matrix mu_s, // [K, size_s_mu]
    array[] int idx_event,
    array[,] int idx_quad,
    int Q,
    vector x,
    vector w
  ){
  
    int N = num_elements(t);
    real log_intensity = 0;
    real compensator = 0;
    
    matrix[K, K] alpha_scaled;

    for (l in 1:K) {
      for (k in 1:K) {
        alpha_scaled[l, k] = alpha[l, k] / eta[l];
      }
    }
    
    // =========================================================================
    // R vector
    // =========================================================================
    // vector R for excitation
    // R_g[l] at event i:
    // sum_{j: t_j < t_i, m[j] = l} exp(-(t[i] - t[j]) / eta[l])
    
    // vector R for inhibition
    // R_h[l] at event i:
    // sum_{j: t_j < t_i, m[j] = l} exp(-(t[i] - t[j]) / phi[l])
    
    vector[K] R_g = rep_vector(0, K);
    vector[K] R_h = rep_vector(0, K);
    
    //[0 ~ t1]
    {
      real D = t[1];
      real height = 0;
      vector[K] R_zero = rep_vector(0, K);
      
      for(q in 1 : Q){
        real u_q = 0.5 * D * (x[q] + 1);
        vector[K] mu_q = col(mu_s, idx_quad[1, q]);
        
        height += w[q]* lambda_total_u(
          u_q, R_zero, R_zero, alpha_scaled, eta, gamma, phi, mu_q, K
        );
      }
      
      compensator += 0.5 * D * height;
    }
    
    // ========================================================================= 
    // Main loop
    // =========================================================================
    for(i in 1:N){
    
      int k_now = m[i];
      
      
      // update the accumulated term from t[i-1] to t[i]
      if(i > 1){
        real dt = t[i] - t[i - 1];

        for(l in 1:K){
          R_g[l] *= exp(-dt / eta[l]);
          R_h[l] *= exp(-dt / phi[l]);
        }
      }
      
      // =======================================================================
      // Log-Intensity
      // =======================================================================
    
      vector[K] mu_i = col(mu_s, idx_event[i]);
      real excitation = mu_i[k_now];
      real inhibition = 0;
      
      for(l in 1:K){
        excitation += alpha_scaled[l, k_now] * R_g[l];
        inhibition += gamma[l, k_now] * R_h[l];
      }
      
      log_intensity += log(excitation) - inhibition;
      
      
      // =======================================================================
      // update R after observing current event
      // current event group becomes a past source group for future events
      // =======================================================================
    
      int l_future = k_now;
      R_g[l_future] += 1;
      R_h[l_future] += 1;
      
      // =======================================================================
      // Gauss-Legendre Quadrature (Numerical Integration)
      // =======================================================================
      
      real D;
      real height = 0;
      if (i < N) {
        D = t[i + 1] - t[i];
      } else {
        D = T_end - t[i];
      }
      
      for(q in 1:Q){
        
        real u_q = 0.5 * D * (x[q] + 1);
        vector[K] mu_q = col(mu_s, idx_quad[i + 1, q]);
        
        height += w[q] * lambda_total_u(
          u_q, R_g, R_h, alpha_scaled, eta, gamma, phi, mu_q, K
        );
      }
      
      compensator += 0.5 * D * height;
      
    }
    
    // =========================================================================
    // Return
    // =========================================================================
  
    return log_intensity - compensator;
  }
}

data {
  int<lower=1> N;
  int<lower=1> K;

  vector<lower=0>[N] t;
  array[N] int<lower=1, upper=K> m;

  real<lower=0> T_end;

  // Normal prior scale for a_star
  real<lower=0> a_sd;
  
  // mu_k(t) 평가가 필요한 점들의 전체 시간 배열
  int<lower=2> size_s_mu;
  vector<lower=0>[size_s_mu] s_mu;
  
  // observed event t[i]가 s_mu에서 어떤 인덱스?
  array[N] int<lower=1, upper=size_s_mu> idx_event;
  
  // N+1개 구간, 각 구간 당 8개 node
  // r번째 구간 (t[r-1], t[r])에서의 q번째 노드의
  // Absolute value에 대한, s_mu에서의 인덱스?
  array[N + 1, 8] int<lower=1, upper=size_s_mu> idx_quad;
  
  real<lower=0> ell_mu; // Length scale
}

transformed data {
  int Q = 8;
  vector[8] x;
  vector[8] w;
  vector[size_s_mu - 1] rho_mu;
  vector[size_s_mu - 1] innovation_scale_mu;

  // 8-point Gauss-Legendre nodes on [-1, 1]
  x[1] = -0.9602898564975363;
  x[2] = -0.7966664774136267;
  x[3] = -0.5255324099163290;
  x[4] = -0.1834346424956498;
  x[5] =  0.1834346424956498;
  x[6] =  0.5255324099163290;
  x[7] =  0.7966664774136267;
  x[8] =  0.9602898564975363;

  // 8-point Gauss-Legendre weights
  w[1] = 0.1012285362903763;
  w[2] = 0.2223810344533745;
  w[3] = 0.3137066458778873;
  w[4] = 0.3626837833783620;
  w[5] = 0.3626837833783620;
  w[6] = 0.3137066458778873;
  w[7] = 0.2223810344533745;
  w[8] = 0.1012285362903763;
  
  for(r in 1:(size_s_mu - 1)){
    real delta_s = s_mu[r+1] - s_mu[r];
    
    rho_mu[r] = 
      exp(-delta_s / ell_mu);
    
    innovation_scale_mu[r] = 
      sqrt(-expm1(-2.0 * delta_s / ell_mu));
  }
}


parameters {

  vector[K] lmu_mean;
  vector<lower=0>[K] sigma_mu;
  matrix[K, size_s_mu] z_mu;

  // a_star > 0  -> excitation
  // a_star < 0  -> inhibition
  matrix<lower=-1, upper=1>[K, K] a_star;

  vector[K] leta;
  vector[K] lphi;
}


transformed parameters {
  // exp 씌워진 mu값
  matrix[K, size_s_mu] lmu_s;
  matrix[K, size_s_mu] mu_s;
  
  vector<lower=0>[K] eta;
  vector<lower=0>[K] phi;
  
  for(k in 1 : K){
    lmu_s[k, 1] = 
    lmu_mean[k] + sigma_mu[k] * z_mu[k, 1];
    
    for(r in 1 : (size_s_mu - 1)){
      lmu_s[k, r+1] = 
        lmu_mean[k] 
        + rho_mu[r] * (lmu_s[k, r] - lmu_mean[k])
        + sigma_mu[k]
          * innovation_scale_mu[r]
          * z_mu[k, r+1];
        
    }
  }
  

  mu_s = exp(lmu_s);
  eta = exp(leta);
  phi = exp(lphi);
}

model {

  matrix[K, K] alpha;
  matrix[K, K] gamma;

  // ===========================================================================
  // alpha/gamma construction 
  // alpha[l,k] = max(a_star[l,k], 0)
  // gamma[l,k] = max(-a_star[l,k], 0)
  // ===========================================================================
  for (l in 1:K) {
    for (k in 1:K) {
      alpha[l, k] = a_star[l, k] > 0 ? a_star[l, k] : 0;
      gamma[l, k] = a_star[l, k] < 0 ? -a_star[l, k] : 0;
    }
  }


  // ===========================================================================
  // Priors
  // ===========================================================================
  
  lmu_mean ~ normal(-6, 2);
  sigma_mu ~ normal(0, 0.5);
  to_vector(z_mu) ~ std_normal();
  
  to_vector(a_star) ~ normal(0, a_sd);

  leta ~ normal(6, 1.5); // set the parameters manually
  lphi ~ normal(6, 1.5); // set the parameters manually

  // ===========================================================================
  // Likelihood
  // ===========================================================================
  
  target += MHPWI_loglike_lpdf(
    t | m, K, T_end, alpha, eta, gamma, phi, mu_s, idx_event, idx_quad, Q, x, w
  );
  
}

generated quantities {
  matrix<lower=0, upper=1>[K, K] alpha;
  matrix<lower=0, upper=1>[K, K] gamma;
  matrix[K, K] signed_effect;

  matrix<lower=0, upper=1>[K, K] is_excitation;
  matrix<lower=0, upper=1>[K, K] is_inhibition;

  for (l in 1:K) {
    for (k in 1:K) {
      alpha[l, k] = a_star[l, k] > 0 ? a_star[l, k] : 0;
      gamma[l, k] = a_star[l, k] < 0 ? -a_star[l, k] : 0;

      signed_effect[l, k] = alpha[l, k] - gamma[l, k];

      is_excitation[l, k] = a_star[l, k] > 0 ? 1 : 0;
      is_inhibition[l, k] = a_star[l, k] < 0 ? 1 : 0;
    }
  }
}
"


#################
# =============================================================================
# Time grid for time-varying mu_k(t)
# =============================================================================
T_end = 86400
t_num <- as.numeric(t_num)
T_end <- as.numeric(T_end)

N <- length(t_num)
Q <- 8L

stopifnot(
  N >= 1L,
  all(is.finite(t_num)),
  is.finite(T_end),
  all(t_num >= 0),
  all(diff(t_num) > 0),
  T_end > t_num[N]
)

# 8-point Gauss-Legendre nodes
x_gl <- c(
  -0.9602898564975363,
  -0.7966664774136267,
  -0.5255324099163290,
  -0.1834346424956498,
  0.1834346424956498,
  0.5255324099163290,
  0.7966664774136267,
  0.9602898564975363
)

# 8-point Gauss-Legendre weights
w_gl <- c(
  0.1012285362903763,
  0.2223810344533745,
  0.3137066458778873,
  0.3626837833783620,
  0.3626837833783620,
  0.3137066458778873,
  0.2223810344533745,
  0.1012285362903763
)

interval_left  <- c(0, t_num)
interval_right <- c(t_num, T_end)

interval_length <- interval_right - interval_left

stopifnot(
  length(interval_left) == N + 1L,
  length(interval_right) == N + 1L,
  all(interval_length > 0)
)

quad_times <- matrix(
  NA_real_,
  nrow = N + 1L,
  ncol = Q
)

for (r in seq_len(N + 1L)) {
  quad_times[r, ] <-
    interval_left[r] +
    0.5 * interval_length[r] * (x_gl + 1)
}

raw_s_mu <- c(
  0,
  t_num,
  as.vector(t(quad_times))
)

s_mu <- sort(unique(raw_s_mu))

size_s_mu <- length(s_mu)

idx_event <- match(t_num, s_mu)

idx_event <- as.integer(idx_event)

idx_quad <- matrix(
  match(
    as.vector(t(quad_times)),
    s_mu
  ),
  nrow = N + 1L,
  ncol = Q,
  byrow = TRUE
)

storage.mode(idx_quad) <- "integer"

##################


file_MHPWI = write_stan_file(MHPWI_model_code)
MHPWI_model = cmdstan_model(file_MHPWI)

sim_dat1 = read.csv("simulated_events.csv")
t_num = sim_dat1$time
m = sim_dat1$mark
K = 3

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

eta_true <- c(1342.8, 251.4, 149.4)
phi_true <- c(158.4, 1509.0, 164.4)
ell_mu_true <- 10800
sigma_mu_true <- c(0.30, 0.30, 0.30)

stopifnot(
  ell_mu_true > max(eta_true),
  ell_mu_true > max(phi_true),
  T_end >= 3 * ell_mu_true
)

## construct Stan data ----
stan_data <- list(
  N = as.integer(N),
  K = as.integer(K),
  
  t = t_num,
  m = as.integer(m),
  
  T_end = T_end,
  a_sd = 0.5,
  
  size_s_mu = as.integer(size_s_mu),
  s_mu = s_mu,
  
  idx_event = idx_event,
  idx_quad = idx_quad,
  
  ell_mu = ell_mu_true
)
range(diff(stan_data$s_mu))
min(diff(stan_data$s_mu))
summary(diff(stan_data$s_mu))
stan_data$ell_mu
is.unsorted(s_mu)


dir.create("stan_outputs", showWarnings = FALSE)

saveRDS(fit_400, "fit_400.RDS")

## fit the model ---- 
fit_6000 <- MHPWI_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  adapt_delta = 0.9,
  seed = 123,
  refresh = 10,
  output_dir = "stan_outputs"
)


diag <- fit_2000$diagnostic_summary()
summary_2000 <- fit_2000$summary()
summary_6000 <- fit_6000$summary()
fit_6000$diagnostic_summary()

summary_6000 %>%
  filter(str_detect(variable, "eta"))

draws <- posterior::as_draws_df(
  fit_6000$draws(
    variables = c(
      "alpha[1,1]",
      "a_star[1,1]",
      "leta[1]",
      "eta[1]",
      "sigma_mu[1]"
    )
  )
)

cor(
  draws[["alpha[1,1]"]],
  draws[["leta[1]"]]
)

cor(
  draws[["alpha[1,1]"]],
  draws[["sigma_mu[1]"]]
)

idx <- draws[["a_star[1,1]"]] > 0

cor(
  draws[["a_star[1,1]"]][idx],
  draws[["leta[1]"]][idx]
)
a11 <- draws[["a_star[1,1]"]]
idx <- a11 > 0

c(
  prob_excitation = mean(idx),
  mean_alpha = mean(pmax(a11, 0)),
  mean_given_excitation = mean(a11[idx]),
  median_given_excitation = median(a11[idx])
)
jump11 <- draws[["alpha[1,1]"]] / draws[["eta[1]"]]

quantile(
  jump11,
  c(0.05, 0.5, 0.95)
)
idx <- draws[["a_star[1,1]"]] > 0

jump11_positive <-
  draws[["a_star[1,1]"]][idx] /
  draws[["eta[1]"]][idx]

c(
  mean = mean(jump11_positive),
  median = median(jump11_positive),
  quantile(
    jump11_positive,
    c(0.05, 0.5, 0.95)
  )
)
