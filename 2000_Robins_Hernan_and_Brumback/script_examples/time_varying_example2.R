# Testing if Adjusting for Time-Varying Confounders via IPTW Works
# 
# Consider the fundamental difficulty that if there is time-varying confounding
# in our simulated data generating mechanism, we do not have a way of a priori
# telling if our causal estimate is correct or not, at least to my knowledge. 
# 
# What should be true, however, is that IPTW should produce a population that 
# is as if the time-varying confounding had not occured. 
# 
# If we were to have a parameter "amount_of_confounding" that varied, we would 
# expect to see that when we ramp this up in the positive direction we would see
# bias in one direction in a naive analysis, and similarly we would expect to
# see bias in the naive analysis in the opposite direction if this parameter
# is ramped up in the negative direction. If the "amount_of_confounding" is 0,
# then we should expect to see that the naive analysis matches that of a 
# IPTW analysis. 

# Our DAG, drawn here for reference
# This is the more complicated scenario:
# 
#   _ _ _ _ _ _
#  /           \
# L0    →  L1   \
#   ↘   ↗↘  ↓  ↘ ↘ 
#     A0  →  A1  →  Y
#      \          ↗
#       \ _ _ _ _/
# 
# Though I think we can simplify to a less complicated but still 
# time varying confounded scenario: 
# 
# L0    →  L1
#   ↘    ↗↘  ↓  ↘ 
#     A0  →  A1  →  Y
#      \          ↗
#       \ _ _ _ _/

# Shu described to me that the way we should think about IPTW is as 
# removing arrows from treatment to exposures, 
# That would imply we would think of the limiting case as: 
# 
# L0    →  L1
#        ↗      ↘ 
#     A0  →  A1  →  Y
#      \          ↗
#       \ _ _ _ _/
# 
# So I should be able to do simpler causal/mediation analysis in this limiting
# setting to derive the causal effects I expect to get after IPTW
# 
# and Shu suggests that I want to toggle/ramp-up/down is the strength of the 
# arrows from the confounders into the treatment as my "amount_of_confounding"


# Shu points out we could consider a deterministic scenario: 
# 
# A0 increases Y by 1 
# A1 increases Y by 1 
# A0 increases A1 by 1 
# and then A0 also increases L1 by 1 and L1 increases A1 by 1 and Y by 1 
# and then L0 could increase A0 and L1 by 1 
# 
# In the unconfounded setting, we should have that a 1-unit increase in A0 has: 
# - 1 unit increase directly on Y, unmediated 
# - 1 unit increase on L1 
# - 1 unit increase directly on A1, unmediated 
# - the 1 unit increase on L1 cascades to another 1 unit increase on Y 
# - and the 1-unit increase of A1 increases Y directly, unmediated by 1
#
# So we should get an increase of Y by 3 for every 1-unit increase in A0.
# 
# In the confounded setting, we should have that a 1-unit increase in A0 has:
# - 1 unit increase directly on Y, unmediated 
# - 1 unit increase on L1 
# - 1 unit increase directly on A1, unmediated 
# - the 1 unit increase on L1 cascades to another 1 unit increase on A1 
# - the 1 unit increase on L1 cascades to another 1 unit increase on Y 
# - and the 2-unit increase of A1 increases Y directly, unmediated by 2 
#
# So we should get an increase of Y by 4 for every 1-unit increase in A0.
# 
# Jarvis points out that the reason we're getting an estimate of 5.5 below 
# is because of the back-door paths: 
#   A0 -- L0 -- L1 -- Y
# and  A0 -- L0 -- L1 -- A1 -- Y 
# 
# Jarvis proposes renaming "amount of confounding" -> "confounding_path_strength" 
# where the amount of confounding should be a measure of 
# the bias in the true A0 effect on Y vs. the observed association 
# or otherwise some function of all of the paths. 

# -----------------------------------------------------------------------

library(tidyverse)
library(patchwork)

# Function to generate simulated data
generate_sim_data <- function(N, confounding_path_strength) {
  L0 <- rnorm(N)
  A0 <- rnorm(N, mean = L0 * confounding_path_strength)
  L1 <- rnorm(N, mean = L0 + A0)
  A1 <- rnorm(N, mean = A0 + L0 * confounding_path_strength + L1 * confounding_path_strength)
  Y <- rnorm(N, mean = A0 + A1 + L1)
  
  data.frame(L0 = L0, A0 = A0, L1 = L1, A1 = A1, Y = Y)
}

# Function to calculate IPTW
calculate_weights <- function(data) {
  model_A0 <- lm(A0 ~ L0, data = data)
  model_A1 <- lm(A1 ~ A0 + L0 + L1, data = data)
  
  weights_A0 <- 1 / dnorm(x = data$A0, mean = coef(model_A0)[1], sd = summary(model_A0)$sigma)
  weights_A1 <- 1 / dnorm(
    x = data$A1,
    mean = coef(model_A1)[1] + coef(model_A1)[2] * data$A0 +
      coef(model_A1)[3] * data$L0 + coef(model_A1)[4] * data$L1,
    sd = summary(model_A1)$sigma
  )
  
  weights_A0 * weights_A1
}

# Function to calculate stabilized IPTW
calculate_stable_weights <- function(data, use_log = TRUE) {
  model_A0 <- lm(A0 ~ L0, data = data)
  model_A1 <- lm(A1 ~ A0 + L0 + L1, data = data)
  
  model_A0_alone <- lm(A0 ~ 1, data = data)
  model_A1_stabilizer <- lm(A1 ~ A0, data = data)

  if (use_log) {
    inv_log_weight_A0 = dnorm(x = data$A0, mean = coef(model_A0)[1], sd = summary(model_A0)$sigma, log = TRUE)
    inv_log_weight_A1 = dnorm(
      x = data$A1,
      mean = coef(model_A1)[1] + coef(model_A1)[2] * data$A0 +
        coef(model_A1)[3] * data$L0 + coef(model_A1)[4] * data$L1,
      sd = summary(model_A1)$sigma,
      log = TRUE
    )
    
    log_A0_stabilizer <- dnorm(x = data$A0, mean = coef(model_A0_alone)[1], sd = summary(model_A0_alone)$sigma, log=TRUE)
    log_A1_stabilizer <- dnorm(x = data$A1, mean = coef(model_A0_alone)[1] + data$A0*coef(model_A1_stabilizer)[2], sd = summary(model_A1_stabilizer)$sigma, log=TRUE)

    log_stabilized_weight <- log_A0_stabilizer + log_A1_stabilizer - inv_log_weight_A0 - inv_log_weight_A1
    return(exp(log_stabilized_weight))
  } else {
    weights_A0 <- 1 / dnorm(x = data$A0, mean = coef(model_A0)[1], sd = summary(model_A0)$sigma)
    weights_A1 <- 1 / dnorm(
      x = data$A1,
      mean = coef(model_A1)[1] + coef(model_A1)[2] * data$A0 +
        coef(model_A1)[3] * data$L0 + coef(model_A1)[4] * data$L1,
      sd = summary(model_A1)$sigma
    )
    
    model_A0_alone <- lm(A0 ~ 1, data = data)
    model_A1_stabilizer <- lm(A1 ~ A0, data = data)
    A0_stabilizer <- dnorm(x = data$A0, mean = coef(model_A0_alone)[1], sd = summary(model_A0_alone)$sigma)
    A1_stabilizer <- dnorm(x = data$A1, mean = coef(model_A0_alone)[1] + data$A0*coef(model_A1_stabilizer)[2], sd = summary(model_A1_stabilizer)$sigma)
    
    return(A0_stabilizer * A1_stabilizer * weights_A0 * weights_A1)
  }
}

# Generate simulated data
N <- 1000
confounding_path_strength <- .3
sim_data <- generate_sim_data(N, confounding_path_strength)

# Calculate weights
sim_data$weights <- calculate_weights(sim_data)
sim_data$stable_weights <- calculate_stable_weights(sim_data)
sim_data$stable_weights_no_log_precision <- calculate_stable_weights(sim_data, use_log = FALSE)

# check if log precision improves things:
# answer: not really
# ggplot(sim_data, aes(x = weights, y = stable_weights_no_log_precision)) + 
#   geom_point() + 
#   geom_smooth(method = 'lm')

# Unadjusted model
model_unadjusted <- lm(Y ~ A0 + A1, data = sim_data)
summary(model_unadjusted)

# Adjust for confounders using MSM
model_msm <- lm(Y ~ A0 + A1, data = sim_data, weights = sim_data$stable_weights)
summary(model_msm)

make_exposure_confounder_bivariate_plt_weighted <- function(var1, var2) {
  ggplot(sim_data, aes(x = {{var1}}, y = {{var2}}, color = stable_weights, size = stable_weights)) + 
  geom_point(alpha = .5) + 
  scale_color_viridis_c(option = 'magma', end = .8) + 
  theme_bw() +
  guides(color= guide_legend(), size=guide_legend())
}

plt1 <- make_exposure_confounder_bivariate_plt_weighted(A0, L0)
plt2 <- make_exposure_confounder_bivariate_plt_weighted(A1, L0)
plt3 <- make_exposure_confounder_bivariate_plt_weighted(A0, L1)
plt4 <- make_exposure_confounder_bivariate_plt_weighted(A0, A1)
 
(plt1 + plt2 )/ (plt3 + plt4) 

hist(sim_data$stable_weights, breaks = 50)

library(scatterplot3d)
library(rgl)
library(car)

myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

cols <- myColorRamp(c("red", "blue"), sim_data$stable_weights) 

scatter3d(x = sim_data$A0, y = sim_data$A1, z = sim_data$L1, 
  point.col = cols,
  radius = sim_data$stable_weights)

library(plotly)

plot_ly(sim_data, 
  x = ~A0,
  y = ~A1, 
  z = ~L0,
  size = ~stable_weights,
  color = ~stable_weights,
  alpha = .5) %>% 
  add_markers(marker = list(sizeref = .1))

# Look at estimates when confounding_path_strength is ramped up or down 
results <- list()
for (confounding_path_strength in seq(-1,1,length.out = 7)) {
  for (i in 1:100) {
  # Generate simulated data
  N <- 1000
  sim_data <- generate_sim_data(N, confounding_path_strength)
  
  # Calculate weights
  sim_data$stable_weights <- calculate_stable_weights(sim_data)
  
  # Unadjusted model
  model_unadjusted <- lm(Y ~ A0 + A1, data = sim_data)

  # Adjust for confounders using MSM
  model_msm <- lm(Y ~ A0 + A1, data = sim_data, weights = sim_data$stable_weights)

  # store results from MSM model
  results[[length(results)+1]] <- 
    list(confounding_path_strength = confounding_path_strength,
      i = i,
      model = 'marginal structural model (with IPTW)',
      intercept = unname(coef(model_msm)[1]),
      A0_effect = unname(coef(model_msm)[2]),
      A1_effect = unname(coef(model_msm)[3]))
  
  # store results from unadjusted model
  results[[length(results)+1]] <- 
    list(confounding_path_strength = confounding_path_strength,
      i = i,
      model = 'unadjusted',
      intercept = unname(coef(model_unadjusted)[1]),
      A0_effect = unname(coef(model_unadjusted)[2]),
      A1_effect = unname(coef(model_unadjusted)[3]))
  }
}

# turn into a dataframe
results <- bind_rows(results)

labelling_helper <- function(x) { 
  paste0("confounding path coefficient: ", round(as.numeric(x), 1))
}

results |> # filter(confounding_path_strength %in% c(3, 2, 1, 0)) |> 
  tidyr::pivot_longer(cols = c(intercept, A0_effect, A1_effect), names_to = 'coefficient', values_to = 'estimate') |> 
  # mutate(coefficient = factor(coefficient, c('intercept', 'A0_effect', 'A1_effect'))) |> 
  ggplot(aes(x = coefficient, y = estimate, color = model, group = model)) + 
  geom_jitter(height = 0, width=.15, alpha = .15) + 
  stat_summary(fun = mean, geom = 'point', aes(group = model, shape = model), color = 'grey30') + 
  stat_summary(fun = mean, geom = 'line', aes(group = model)) + 
  facet_wrap(~confounding_path_strength, nrow = 1,
             labeller = as_labeller(labelling_helper)) + 
  theme_bw() + 
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 15, hjust = 1),
        plot.caption = element_text(hjust = .5)) + 
  ggtitle("Variation in Coefficient Estimates with/without IPTW",
          subtitle = "Datasets (N obs=1000) were simulated for varying coefficients for the confounding paths (L0→A0, L0→A1, L1→A1) 100 times") + 
  labs(caption = "Mean coefficient estimates from the 100 simulations are indicated by the black marks that are connected by colored lines for each model type.")

ggsave("2000_Robins_Hernan_and_Brumback/script_examples/time_varying_confounding/time_varying_confounding_lt1.png",
       width = 14, height = 3.5)




# Look at estimates when confounding_path_strength is ramped up or down 
results <- list()
for (confounding_path_strength in seq(-2,2,length.out = 7)) {
  for (i in 1:100) {
    # Generate simulated data
    N <- 1000
    sim_data <- generate_sim_data(N, confounding_path_strength)
    
    # Calculate weights
    sim_data$stable_weights <- calculate_stable_weights(sim_data)
    
    # Unadjusted model
    model_unadjusted <- lm(Y ~ A0 + A1, data = sim_data)
    
    # Adjust for confounders using MSM
    model_msm <- lm(Y ~ A0 + A1, data = sim_data, weights = sim_data$stable_weights)
    
    # store results from MSM model
    results[[length(results)+1]] <- 
      list(confounding_path_strength = confounding_path_strength,
           i = i,
           model = 'marginal structural model (with IPTW)',
           intercept = unname(coef(model_msm)[1]),
           A0_effect = unname(coef(model_msm)[2]),
           A1_effect = unname(coef(model_msm)[3]))
    
    # store results from unadjusted model
    results[[length(results)+1]] <- 
      list(confounding_path_strength = confounding_path_strength,
           i = i,
           model = 'unadjusted',
           intercept = unname(coef(model_unadjusted)[1]),
           A0_effect = unname(coef(model_unadjusted)[2]),
           A1_effect = unname(coef(model_unadjusted)[3]))
  }
}

# turn into a dataframe
results <- bind_rows(results)

labelling_helper <- function(x) { 
  paste0("confounding path coefficient: ", round(as.numeric(x), 1))
}

results |> # filter(confounding_path_strength %in% c(3, 2, 1, 0)) |> 
  tidyr::pivot_longer(cols = c(intercept, A0_effect, A1_effect), names_to = 'coefficient', values_to = 'estimate') |> 
  # mutate(coefficient = factor(coefficient, c('intercept', 'A0_effect', 'A1_effect'))) |> 
  ggplot(aes(x = coefficient, y = estimate, color = model, group = model)) + 
  geom_jitter(height = 0, width=.15, alpha = .15) + 
  stat_summary(fun = mean, geom = 'point', aes(group = model, shape = model), color = 'grey30') + 
  stat_summary(fun = mean, geom = 'line', aes(group = model)) + 
  facet_wrap(~confounding_path_strength, nrow = 1,
             labeller = as_labeller(labelling_helper)) + 
  theme_bw() + 
  theme(legend.position = 'bottom', axis.text.x = element_text(angle = 15, hjust = 1)) + 
  ggtitle("Variation in Coefficient Estimates with/without IPTW",
          subtitle = "Datasets (N obs=1000) were simulated for varying coefficients for the confounding paths (L0→A0, L0→A1, L1→A1) 100 times")

ggsave("2000_Robins_Hernan_and_Brumback/script_examples/time_varying_confounding/time_varying_confounding_wider_scenario.png",
       width = 14, height = 3.5)

