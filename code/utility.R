library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)

plot_trace_with_true <- function(fit, vars, true_vals){
  if(length(vars) != length(true_vals)){
    stop("vars와 true_vals의 길이가 같아야 합니다.")
  }
  
  draw <- fit$draws(variables = vars)
  
  true_df <- data.frame(
    parameter = vars,
    true_value = as.numeric(true_vals)
  )
  
  mcmc_trace(draws, pars = vars) +
    geom_hline(
      data = true_df,
      aes(yintercept = true_value),
      linetype = "dashed",
      linewidth = 0.4,
      inherit.aes = FALSE
    )
}

summarise_posterior_with_truth <- function(
    fit,
    vars,
    true_vals,
    probs = c(0.025, 0.975),
    digits = 4
) {
  if (length(vars) != length(true_vals)) {
    stop("vars와 true_vals의 길이가 같아야 합니다.")
  }
  
  if (anyDuplicated(vars)) {
    stop("vars에 중복된 파라미터 이름이 있습니다.")
  }
  
  if (
    length(probs) != 2 ||
    any(probs < 0) ||
    any(probs > 1) ||
    probs[1] >= probs[2]
  ) {
    stop("probs는 0과 1 사이의 오름차순 길이 2 벡터여야 합니다.")
  }
  
  true_df <- data.frame(
    Parameter = vars,
    true_value = as.numeric(true_vals)
  )
  
  fit$draws(
    variables = vars,
    format = "draws_df"
  ) %>%
    as.data.frame() %>%
    select(all_of(vars)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "Parameter",
      values_to = "draw"
    ) %>%
    group_by(Parameter) %>%
    summarise(
      mean_estimate = mean(draw),
      ci_lower = quantile(draw, probs[1]),
      ci_upper = quantile(draw, probs[2]),
      .groups = "drop"
    ) %>%
    left_join(true_df, by = "Parameter") %>%
    mutate(
      coverage = true_value >= ci_lower &
        true_value <= ci_upper,
      
      `True value` = round(true_value, digits),
      
      `Mean estimate` = round(mean_estimate, digits),
      
      `credible interval` = sprintf(
        paste0("[%.", digits, "f, %.", digits, "f]"),
        ci_lower,
        ci_upper
      )
    ) %>%
    select(
      Parameter,
      `True value`,
      `Mean estimate`,
      `credible interval`,
      coverage
    )
}