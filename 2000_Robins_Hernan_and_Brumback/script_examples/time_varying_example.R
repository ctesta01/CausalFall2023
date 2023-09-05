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

# dependencies ------------------------------------------------------------

library(tidyverse)

mostly_deterministic_simulation <- function(N = 1000, L0 = 0, A0_beyond_L0 = 0, confounding_path_strength = 0) {
  A0 <- L0*confounding_path_strength + A0_beyond_L0 + rnorm(n = N)
  L1 <- L0 + A0
  A1 <- A0 + L0*confounding_path_strength + L1*confounding_path_strength
  Y <- A0 + A1 + L1
  
  return(data.frame(L0 = L0, A0 = A0, L1 = L1, A1 = A1, Y = Y))
}
lm(Y ~ A0, data = mostly_deterministic_simulation(L0 = 0, A0_beyond_L0 = 0))
lm(Y ~ A0, data = mostly_deterministic_simulation(L0 = 0, A0_beyond_L0 = 1))
lm(Y ~ A0, mostly_deterministic_simulation(L0 = 0, A0_beyond_L0 = 0, confounding_path_strength = 1))
lm(Y ~ A0, mostly_deterministic_simulation(L0 = 0, A0_beyond_L0 = 1, confounding_path_strength = 1))

# create simulated data with a time-varying confounder affected by the exposure
simulate_data <- function(
    N = 1000, # population size 
    confounding_path_strength = 0,
    sd = 1
    ) {
  L0 <- rnorm(n = N, sd = sd)
  A0 <- rnorm(n = N, mean = L0*confounding_path_strength, sd = sd)
  L1 <- rnorm(n = N, mean = L0 + A0, sd = sd)
  A1 <- rnorm(n = N, mean = A0 + L0*confounding_path_strength + L1*confounding_path_strength, sd = sd)
  Y <- rnorm(n = N, mean = A0 + A1 + L1, sd = sd)
  cumulative_A <- A0 + A1
  
  tibble(L0 = L0, A0 = A0, L1 = L1, A1 = A1, cumulative_A, Y = Y)
}

df <- simulate_data(N = 10000, confounding_path_strength = 0, sd = 1)

lm(formula = Y ~ cumulative_A, data = df)

df <- simulate_data(N = 10000, confounding_path_strength = 1, sd = 1)

lm(formula = Y ~ A0, data = df)
lm(formula = Y ~ cumulative_A, data = df)

A0.L0.model <- lm(A0 ~ L0, data = df)
weights <- 1 / dnorm(
  x = df$A0, 
  mean = coef(A0.L0.model)[1] + coef(A0.L0.model)[2]*df$L0,
  sd = summary(A0.L0.model)$sigma)

lm(formula = Y ~ A0, data = df, weights = weights)

A0.model <- lm(A0 ~ 1, data = df)
stabilized_weights <- 
  dnorm(x = df$A0,
        mean = coef(A0.model)[1],
        sd = summary(A0.model)$sigma) / 
  dnorm(
    x = df$A0, 
    mean = coef(A0.L0.model)[1] + coef(A0.L0.model)[2]*df$L0,
    sd = summary(A0.L0.model)$sigma)

lm(formula = Y ~ A0, data = df, weights = stabilized_weights)
lm(formula = Y ~ cumulative_A, data = df, weights = stabilized_weights)

A1.A0.L0.L1.model <- lm(A1 ~ A0 + L0 + L1, data = df)

multi_time_weights <- 
  1 / (dnorm(
    x = df$A0, 
    mean = coef(A0.L0.model)[1] + coef(A0.L0.model)[2]*df$L0,
    sd = summary(A0.L0.model)$sigma) * 
      dnorm(
        x = df$A1,
        mean = coef(A1.A0.L0.L1.model)[1] + coef(A1.A0.L0.L1.model)[2]*df$A0 + 
          coef(A1.A0.L0.L1.model)[3]*df$L0 + coef(A1.A0.L0.L1.model)[4]*df$L1,
        sd = summary(A1.A0.L0.L1.model)$sigma
      ))

lm(formula = Y ~ A0, data = df, weights = multi_time_weights)
lm(formula = Y ~ cumulative_A, data = df, weights = multi_time_weights)


A1.A0.model <- lm(A1 ~ A0, data = df)

multi_time_stabilized_weights <- 
  (dnorm(x = df$A0,
        mean = coef(A0.model)[1],
        sd = summary(A0.model)$sigma) * 
     dnorm(x = df$A1,
           mean = coef(A1.A0.model)[1] + coef(A1.A0.model)[2]*df$A0,
           sd = summary(A1.A0.model)$sigma)
     )/ (dnorm(
    x = df$A0, 
    mean = coef(A0.L0.model)[1] + coef(A0.L0.model)[2]*df$L0,
    sd = summary(A0.L0.model)$sigma) * 
      dnorm(
        x = df$A1,
        mean = coef(A1.A0.L0.L1.model)[1] + coef(A1.A0.L0.L1.model)[2]*df$A0 + 
          coef(A1.A0.L0.L1.model)[3]*df$L0 + coef(A1.A0.L0.L1.model)[4]*df$L1,
        sd = summary(A1.A0.L0.L1.model)$sigma
      ))

lm(formula = Y ~ A0, data = df, weights = multi_time_stabilized_weights)
lm(formula = Y ~ cumulative_A, data = df, weights = multi_time_stabilized_weights)

