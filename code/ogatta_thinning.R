ogata_thinning <- function(log_mu_grid, cseg, trpar, maxT){
  alpha <- trpar$alpha
  gamma <- trpar$gamma
  eta <- trpar$eta
  phi <- trpar$phi
  
  K <- nrow(alpha)
  
  alpha_eta_ratio <- sweep(alpha, 1, eta, "/") # alpha[l,k]를 eta[l]로 나눔
   
  baseline_at <- function(tt){ # 주어진 점 tt에서의 logmu(tt)를 격자점 기준 선형보간, mu(tt) 추정점 반환환
    left <- findInterval(tt, cseg, all.inside = T)
    
    if(left >= length(cseg)){
      return(exp(log_mu_grid[nrow(log_mu_grid), ]))
    }
    
    weight <- (tt - cseg[left]) / (cseg[left + 1] - cseg[left])
    lmu_tt <- (1-weight) * log_mu_grid[left, ] + weight * log_mu_grid[left + 1,]
    
    exp(lmu_tt)
  }
  
  
  maxmu <- exp(apply(log_mu_grid, 2, max)) # 각 mark K 별로 격자점에서의 mu 최댓값 구함
  alpha_row_total <- rowSums(alpha_eta_ratio) # 각 mark l 별로 alpha/eta 더함
  
  t_now <- 0
  events <- vector("list", K)
  for(k in 1 : K) events[[k]] <- numeric(0)
  
  all_time <- numeric(0)
  all_mark <- integer(0)
  
  while(t_now < maxT){
    Rg_now <- numeric(K)
    
  ################# t 이후에서의 intensity upper bound 설정
    
    for(l in 1 : K){ # 만약 이전에 l이 유발시킨 사건이 있어, remaining excitation 더할게 있는 경우
      if(length(events[[l]])  > 0){
        ages <- t_now - events[[l]] # 각 이벤트의 elapsed time
        Rg_now[l] <- sum(exp(-ages/eta[l])) # 각 이전 사건의 남은 exp항 excitation
      }
    }
    
    maxintensity <- sum(maxmu) + sum(alpha_row_total * Rg_now) 
    # Inhibition은 Intensity 감소시키므로 Upper bound 여부에 영향 주지 않음
    # mu(t)가 cseg 점 기준으로 선형보간 되므로, log_mu_grid 점 들 중 하나에서 최댓값을 가짐
    # Excitation 항은 alpha[l,k] /eta[l] * sum(이전 exp)이므로, alpha[l, k] / eta[l] <= sum(alpha[l,k]) / eta[l]
    
    ################## exp(max) 샘플링
    delta <- rexp(1, rate = maxintensity)
    t_candidate <- t_now + delta
    
    if(t_candidate > maxT) break
    
    ################## intensity at t_candidate과 max_intensity 비교교
    Rg <- numeric(K)
    Rh <- numeric(K)
    
    for(l in 1 : K){ # 각 l 별로 
      if(length(events[[l]]) > 0){
        ages <- t_candidate - events[[l]]
        Rg[l] <- sum(exp(-ages / eta[l]))
        Rh[l] <- sum(exp(-ages / phi[l]))
      }
    }
    
    mu_tt <- baseline_at(t_candidate)
    excitation <- mu_tt + as.vector(crossprod(alpha_eta_ratio, Rg)) 
    # [K x 1] 벡터, crossprod 결과의 각 entry는 lambda_k(t_candidate)의 Excitation 부분에서서 background intensity를 제외한 형태
    inhibition <- as.vector(crossprod(gamma, Rh))
    
    lambda <- excitation * exp(-inhibition) # 각 lambda_k(t_candidate)
    
    total_lambda <- sum(lambda) # 시각 t_candidate에서의 전체 intensity
    
    if (total_lambda > maxintensity * (1 + 1e-8)) {
      stop("Thinning upper bound violated.")
    }
    
    if(total_lambda <= 0){
      stop("intensity should be positive number.")
    }
    
    acceptance_threshold <- total_lambda / maxintensity
    unif <- runif(1)
    
    if(unif < acceptance_threshold){
      all_time <- c(all_time, t_candidate)
      
      candidate_mark <- sample(1:K, size = 1, prob = lambda / total_lambda)
      
      events[[candidate_mark]] <- c(events[[candidate_mark]], t_candidate)
      all_mark <- c(all_mark, candidate_mark)
      
      if(length(all_time) %% 100 == 0){
        message(length(all_time), " points generated so far.")
      }
    }
    
    t_now <- t_candidate
  }
  
  data.frame(
    time = all_time,
    mark = all_mark
  )
  
}





















