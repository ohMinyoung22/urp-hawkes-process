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
  // real<lower=0> a_sd;
  
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
  
  real<lower=0> tau0;
  real<lower=0> slab_scale_effect; // s
  real<lower=0> slab_df_effect; // nu
  
  real<lower=1> hs_df_local;
  real<lower=1> hs_df_global;
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

  vector[K] leta;
  vector[K] lphi;
  
  matrix[K, K] horseshoe_z_effect;
  matrix<lower=0>[K, K] horseshoe_lambda;
  real<lower=0> horseshoe_tau;
  real<lower=0> horseshoe_c2;
}


transformed parameters {
  matrix[K, K] a_star;
  // exp 씌워진 mu값
  matrix[K, size_s_mu] lmu_s;
  matrix[K, size_s_mu] mu_s;
  vector<lower=0>[K] eta;
  vector<lower=0>[K] phi;
  
  for(l in 1:K){
    for(k in 1:K){
      real lambda2;
      real tau2;
      real lambda_tilde;
      
      lambda2 = square(horseshoe_lambda[l, k]);
      tau2 = square(horseshoe_tau);
      
      lambda_tilde = sqrt(
        horseshoe_c2 * lambda2 / (horseshoe_c2 + tau2 * lambda2)
      );
      
      a_star[l,k] = horseshoe_tau * lambda_tilde * horseshoe_z_effect[l, k];
    }
  }
  
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
  to_vector(horseshoe_z_effect) ~ std_normal();
  to_vector(horseshoe_lambda) ~ student_t(hs_df_local, 0, 1);
  horseshoe_tau ~ student_t(hs_df_global, 0, tau0);
  
  horseshoe_c2 ~ inv_gamma(
    0.5 * slab_df_effect,
    0.5 * slab_df_effect * square(slab_scale_effect)
  );

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
  
  // lmu_mean misspecification
  lmu_mean ~ normal(-3, 0.5);
  sigma_mu ~ normal(0, 0.3);
  
  leta ~ normal(log(3.5), 0.35);
  lphi ~ normal(log(3.5), 0.35);
  
  to_vector(z_mu) ~ std_normal();
  //to_vector(a_star) ~ normal(0, a_sd);

  // ===========================================================================
  // Likelihood
  // ===========================================================================
  
  target += MHPWI_loglike_lpdf(
    t | m, K, T_end, alpha, eta, gamma, phi, mu_s, idx_event, idx_quad, Q, x, w
  );
  
}

generated quantities {
  matrix<lower=0>[K, K] alpha;
  matrix<lower=0>[K, K] gamma;
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
setwd("C:/Users/nexen/Desktop/Hawkes_Process/urp-hawkes-process/data")
sim_dat1 = read.csv("simulated_events_1000.csv")
t_num = sim_dat1$time
m = sim_dat1$mark
K = 3

a_star_true <- matrix(
  c(
    0.57, 0.00, -0.26,
    0.00, 0.55, 0.26,
    -0.14, 0.00, 0.73
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

eta_true <- c(5, 3, 2)
phi_true <- c(2, 5, 3)
ell_mu_true <- 30
sigma_mu_true <- c(0.30, 0.30, 0.30)
mean_mu_true <- c(0.2, 0.1, 0.1)
lmu_mean_true <- c(-1.65, -2.35, -2.35)

T_end = 1000
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



## construct Stan data ----
stan_data <- list(
  N = as.integer(N),
  K = as.integer(K),
  
  t = t_num,
  m = as.integer(m),
  
  T_end = T_end,
  
  size_s_mu = as.integer(size_s_mu),
  s_mu = s_mu,
  
  idx_event = idx_event,
  idx_quad = idx_quad,
  
  ell_mu = ell_mu_true,
  
  
  tau0 = 0.05,
  slab_scale_effect = 1,
  slab_df_effect = 4,
  
  hs_df_local = 4,
  hs_df_global = 4
)
range(diff(stan_data$s_mu))
min(diff(stan_data$s_mu))
summary(diff(stan_data$s_mu))
stan_data$ell_mu
is.unsorted(s_mu)

setwd("C:/Users/nexen/Desktop/Hawkes_Process/urp-hawkes-process")
## fit the model ---- 
fit_horse_mis <- MHPWI_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  seed = 123,
  refresh = 10,
  output_dir = "artifact"
)

### Sampler diagnostic
fit_horse_mis$diagnostic_summary()

saveRDS(fit, "fit_normalGP.RDS")
fit$cmdstan_diagnose()

summary <- fit$summary()

### trace plot
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

# 1 <= i, j <= 3인 a_star[i,j] 이름 만들기
vars <- as.vector(outer(
  1:3, 1:3,
  Vectorize(function(i, j) sprintf("a_star[%d,%d]", i, j))
))

# draws 추출
draws <- fit_horse_mis$draws(variables = vars)

# alpha_true가 3x3 matrix라고 가정
true_vals <- as.vector(a_star_true[1:3, 1:3])
names(true_vals) <- vars

# trace plot + true value 수평 점선
p <- mcmc_trace(draws, pars = vars) +
  geom_hline(
    data = data.frame(parameter = vars, a_star_true = true_vals),
    aes(yintercept = a_star_true),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  )

p

# eta[1], eta[2], eta[3] 이름 만들기
vars <- sprintf("lmu_mean[%d]", 1:3)

# draws 추출
draws <- fit_horse_mis$draws(variables = vars)

# eta_true가 길이 3 벡터라고 가정
true_vals <- lmu_mean_true[1:3]

names(true_vals) <- vars

# trace plot + true value 수평 점선
p <- mcmc_trace(draws, pars = vars) +
  geom_hline(
    data = data.frame(parameter = vars, lmu_mean_true = true_vals),
    aes(yintercept = lmu_mean_true),
    linetype = "dashed",
    linewidth = 0.4,
    inherit.aes = FALSE
  )

p


### coverage plot
library(posterior)
library(dplyr)
library(tibble)

lmu_mean_true <-
  log(mean_mu_true) - sigma_mu_true^2 / 2

true_values <- c(
  setNames(
    as.vector(t(a_star_true)),
    paste0(
      "a_star[",
      rep(1:3, each = 3),
      ",",
      rep(1:3, 3),
      "]"
    )
  ),
  setNames(
    sigma_mu_true,
    paste0("sigma_mu[", 1:3, "]")
  ),
  
  setNames(
    eta_true,
    paste0("eta[", 1:3, "]")
  ),
  
  setNames(
    phi_true,
    paste0("phi[", 1:3, "]")
  ),
  
  setNames(
    lmu_mean_true,
    paste0("lmu_mean[", 1:3, "]")
  )
)

library(posterior)
library(dplyr)
library(tidyr)

draws_df <- fit_horse_mis$draws(
  variables = names(true_values),
  format = "draws_df"
)

result <- draws_df %>%
  as.data.frame() %>%
  select(all_of(names(true_values))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "draw"
  ) %>%
  group_by(Parameter) %>%
  summarise(
    `Mean estimate` = mean(draw),
    ci_lower = quantile(draw, 0.025),
    ci_upper = quantile(draw, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    `True value` = true_values[Parameter],
    `credible interval` = paste0(
      "[",
      round(ci_lower, 4),
      ", ",
      round(ci_upper, 4),
      "]"
    ),
    coverage = (
      `True value` >= ci_lower &
        `True value` <= ci_upper
    )
  ) %>%
  select(
    Parameter,
    `True value`,
    `Mean estimate`,
    `credible interval`,
    coverage
  )

print(result, n = Inf)
