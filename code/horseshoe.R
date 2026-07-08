Shoe_model_code = "
functions {
  
  // ===========================================================================
  // lambda_total(u) Function
  // ===========================================================================
  real lambda_total_u(
    real u,
    vector R_g,
    vector R_h,
    matrix alpha_scaled,
    vector eta,
    matrix gamma,
    vector phi,
    vector mu,
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
  
    excitation = mu + alpha_scaled' * Rg_u;
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
    vector mu,
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
    
    compensator += t[1] * sum(mu);
    
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
    
      real excitation = mu[k_now];
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
      if (i < N) {
        D = t[i + 1] - t[i];
      } else {
        D = T_end - t[i];
      }
      real height = 0;
      
      for(q in 1:Q){
        
        real u_q = 0.5 * D * (x[q] + 1);
        
        height += w[q] * lambda_total_u(
          u_q, R_g, R_h, alpha_scaled, eta, gamma, phi, mu, K
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

  // Regularized horseshoe hyperparameters
  real<lower=0> tau0_effect;
  real<lower=0> slab_scale_effect;
  real<lower=0> slab_df_effect;

  // 1이면 Cauchy, 3~4이면 더 안정적인 Student-t
  real<lower=1> hs_df_local;
  real<lower=1> hs_df_global;
}

transformed data {
  int Q = 8;
  vector[8] x;
  vector[8] w;

  // 8-point Gauss-Legendre nodes
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
}

parameters {
  // Background intensity
  vector[K] lmu;

  // Signed regularized horseshoe
  matrix[K, K] z_effect;
  matrix<lower=0>[K, K] lambda_effect;
  real<lower=0> tau_effect;
  real<lower=0> c2_effect;

  // Decay parameters
  vector[K] leta;
  vector[K] lphi;
}

transformed parameters {
  vector<lower=0>[K] mu;
  vector<lower=0>[K] eta;
  vector<lower=0>[K] phi;

  // Latent signed horseshoe coefficient
  matrix[K, K] hs_effect;

  // Actual bounded interaction coefficient
  matrix[K, K] a_star;

  matrix<lower=0, upper=1>[K, K] alpha;
  matrix<lower=0, upper=1>[K, K] gamma;

  mu = exp(lmu);
  eta = exp(leta);
  phi = exp(lphi);

  for (l in 1:K) {
    for (k in 1:K) {
      real lambda2;
      real tau2;
      real lambda_tilde;

      lambda2 = square(lambda_effect[l, k]);
      tau2 = square(tau_effect);

      // Regularized horseshoe local scale
      lambda_tilde =
        sqrt(
          c2_effect * lambda2
          /
          (c2_effect + tau2 * lambda2)
        );

      // Symmetric signed coefficient
      hs_effect[l, k] =
        tau_effect
        * lambda_tilde
        * z_effect[l, k];

      // Preserve the original [-1, 1] restriction
      a_star[l, k] = tanh(hs_effect[l, k]);

      // Positive part: excitation
      // Negative part: inhibition
      alpha[l, k] =
        a_star[l, k] > 0
          ? a_star[l, k]
          : 0;

      gamma[l, k] =
        a_star[l, k] < 0
          ? -a_star[l, k]
          : 0;
    }
  }
}

model {
  // Background and decay priors
  lmu ~ normal(0, 1);
  leta ~ normal(0, 1);
  lphi ~ normal(0, 1);

  // Signed regularized horseshoe
  to_vector(z_effect) ~ std_normal();

  // Half-Student-t because lambda_effect is positive
  to_vector(lambda_effect) ~
    student_t(hs_df_local, 0, 1);

  // Half-Student-t global shrinkage
  tau_effect ~
    student_t(hs_df_global, 0, tau0_effect);

  // Regularized slab
  c2_effect ~ inv_gamma(
    0.5 * slab_df_effect,
    0.5 * slab_df_effect
      * square(slab_scale_effect)
  );

  // Likelihood
  target += MHPWI_loglike_lpdf(
    t |
    m,
    K,
    T_end,
    alpha,
    eta,
    gamma,
    phi,
    mu,
    Q,
    x,
    w
  );
}

"

file_Shoe = write_stan_file(Shoe_model_code)
Shoe_model = cmdstan_model(file_Shoe)

## load data ----

### simulated data with excitation + inhibition
sim_dat1 = readRDS("sim_dat1.rds")

## construct Stan data ----
shoe_data1 <- list(
  N     = nrow(sim_dat1),  # number of events
  t     = sim_dat1$time,   # increasing time sequence
  m     = sim_dat1$mark,   # integer type mark
  T_end = 86400,           # observation end time
  K     = 3,                # number of mark type 
  tau0_effect = 0.05,
  slab_scale_effect = 0.5,
  slab_df_effect = 4,
  
  hs_df_local = 4,
  hs_df_global = 4
)


## fit the model ---- 
fit_shoe1 = Shoe_model$sample(
  data = shoe_data1,      # data
  chains = 4,              # the number of chains
  parallel_chains = 4,     # the number of parallel chains (# of available cpu cores)
  iter_warmup = 10000,      # burn-in iterations
  iter_sampling = 5000,    # sampling iterations
  seed = 123,             # seed
  refresh = 10        # update
)

fit_shoe1$diagnostic_summary()

shoe_summary <- fit_shoe1$summary()
shoe_summary %>% filter(str_detect(variable, "phi"))
fit_naive_summary %>% filter(str_detect(variable, "phi"))
