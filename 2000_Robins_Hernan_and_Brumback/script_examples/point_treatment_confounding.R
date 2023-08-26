

# dependencies ------------------------------------------------------------

library(tidyverse)
library(magrittr)

set.seed(1234)

# Let A = treatment with a drug
# Let L = a predictor of treatment and a predictor of outcome Y 
# Let Y = some successful outcome induced by drug A and covariate L 
# and N is the size of our study population 

# data generation: 
N <- 1000
L <- rbinom(n = N, prob = 0.5, size = 1)
A <- rbinom(n = N, prob = .5 + .3 * L, size = 1) # L affects A

logistic <- function(x) { exp(x) / (1 + exp(x)) } 

beta_0 <- -1
beta_1 <- 0.1 
beta_2 <- 0.5 

# Y is created from both A and L plus random noise
Y <- rbinom(n = N, prob = logistic(beta_0 + beta_1 * L + beta_2 * A), size = 1)

# naive methods -----------------------------------------------------------

summary(glm(Y ~ A, family = binomial(link = 'logit')))
summary(glm(Y ~ A + L, family = binomial(link = 'logit')))


# inverse probability of treatment weighting ------------------------------

treatment_model <- glm(A ~ L, family = binomial(link = 'logit'))
probs_of_treatment <- predict(treatment_model, type='response')
probs_of_treatment <- ifelse(A == 0, 1 - probs_of_treatment, probs_of_treatment) # use Pr(A=0|L=l) for A=0 situations

summary(glm(Y ~ A, family = binomial(link = 'logit'), 
            weights = 1 / probs_of_treatment))



# simulation study --------------------------------------------------------

df <- tibble(i = 1:1000)

beta_0 <- -1
beta_1 <- 0.5
beta_2 <- 0.5

df %<>% rowwise() %>% mutate(
  data = list(
    tibble(
      L = rbinom(n = N, prob = 0.5, size = 1),
      A = rbinom(n = N, prob = .5 + .3*L, size = 1), # L affects A
      Y = rbinom(n = N, prob = logistic(beta_0 + beta_1 * L + beta_2 * A), size = 1))))

df %<>% mutate(
  prob_of_treatment_model = list(glm(A ~ L, family = binomial(link='logit'), data = data)))

df %<>% mutate(
  data = list(bind_cols(
    data, 
    w = ifelse(
      data$A == 1,
      1 / predict(prob_of_treatment_model, type='response'),
      1 / (1 - predict(prob_of_treatment_model, type='response'))
    ))))

df %<>% mutate(
  naive_model_no_L = list(glm(Y ~ A, family = binomial(link = 'logit'), data = data)),
  naive_model_w_L = list(glm(Y ~ A + L, family = binomial(link = 'logit'), data = data)),
  iptw_model = list(glm(Y ~ A, family = binomial(link = 'logit'), data = data, weights = data$w))
  )

df %<>% mutate(
  naive_model_no_L_beta0 = list(coef(naive_model_no_L)[1]),
  naive_model_no_L_beta1 = list(coef(naive_model_no_L)[2]),
  naive_model_w_L_beta0 = list(coef(naive_model_w_L)[1]),
  naive_model_w_L_beta1 = list(coef(naive_model_w_L)[2]),
  iptw_model_beta0 = list(coef(iptw_model)[1]),
  iptw_model_beta1 = list(coef(iptw_model)[2]))

df %<>% ungroup() %>% mutate(across(contains("beta0"), unlist))
df %<>% ungroup() %>% mutate(across(contains("beta1"), unlist))

df %>% 
  select(contains("beta0"), contains("beta1")) %>% 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "model",
                      values_to = "estimate") %>% 
  mutate(intercept_or_slope = ifelse(stringr::str_detect(model, "beta0"), "intercept", "slope")) %>% 
  mutate(model = stringr::str_remove_all(model, "_beta0|_beta1")) %>% 
  ggplot(aes(x = model, y = estimate)) + 
  facet_grid(.~intercept_or_slope) + 
  ggforce::geom_sina(alpha = .05) + 
  geom_boxplot(width = .25, outlier.colour = NA) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# test if claim that Pr(Y=1|A=a) remains unchanged ------------------------

# data generation: 
N <- 1000
L <- rbinom(n = N, prob = 0.5, size = 1)
A <- rbinom(n = N, prob = .5 + .3 * L, size = 1) # L affects A

logistic <- function(x) { exp(x) / (1 + exp(x)) } 

beta_0 <- -1
beta_1 <- 0.1 
beta_2 <- 0.5 

# Y is created from both A and L plus random noise
Y <- rbinom(n = N, prob = logistic(beta_0 + beta_1 * L + beta_2 * A), size = 1)


mean(Y[which(A == 1)]) # Pr(Y=1|A=1)

# now if they are reweighted
weighted.mean(Y[which(A==1)], w = (1 / probs_of_treatment)[which(A==1)])

# and for A == 0?
mean(Y[which(A == 0)]) # Pr(Y=1|A=0)

# reweighted:
weighted.mean(Y[which(A==0)], w = (1 / probs_of_treatment)[which(A==0)])

ALw <- expand.grid(A = c(0,1), L = c(0,1)) 
ALw$prob_of_treatment <- predict(treatment_model, newdata = ALw, type = 'response')
ALw$prob_of_treatment <- ifelse(ALw$A == 0, 1-ALw$prob_of_treatment, ALw$prob_of_treatment) # use Pr(A=0|L=l) for A=0 situations
ALw$w <- 1/ALw$prob_of_treatment

library(tidyverse)
library(showtext)

data.frame(A = A, L = L, Y = Y) |> 
  group_by(A, L, Y) |> 
  tally() |> 
  ggplot(aes(x = L, y = A, label = scales::comma_format()(n))) + 
  geom_tile(alpha = 1, color = 'black', fill = 'white') + 
  geom_text(family = "Georgia") + 
  facet_grid(.~Y, labeller = labeller(Y = \(x) { paste0("Y = ", x) })) + 
  scale_y_continuous(breaks = c(0,1), labels = c("A = 0", "A = 1")) + 
  scale_x_continuous(breaks = c(0,1), labels = c("L = 0", "L = 1")) + 
  labs(x = "", y = "") + 
  ggtitle("Original Population") + 
  theme_bw() + 
  theme(
    legend.position = 'none', panel.grid = element_blank(), 
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    text=element_text(family="Georgia"))

ggsave("2000_Robins_Hernan_and_Brumback/script_examples/point_treatment_confounding/original_table.svg",
       height = 1.75, width=3.5)


data.frame(A = A, L = L, Y = Y) |> 
  group_by(A, L, Y) |> 
  tally() |> 
  left_join(ALw, by = c('A' = 'A', 'L' = 'L')) |> 
  ggplot(aes(x = L, y = A, label = paste0(scales::comma_format()(n), " × ", round(w, 1)))) + 
  geom_tile(alpha = 1, color = 'black', fill = 'white') + 
  geom_text(family = "Georgia") + 
  facet_grid(.~Y, labeller = labeller(Y = \(x) { paste0("Y = ", x) })) + 
  scale_y_continuous(breaks = c(0,1), labels = c("A = 0", "A = 1")) + 
  scale_x_continuous(breaks = c(0,1), labels = c("L = 0", "L = 1")) + 
  labs(x = "", y = "") + 
  ggtitle("Original Population with Weights") + 
  theme_bw() + 
  theme(
    legend.position = 'none', panel.grid = element_blank(), 
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    text=element_text(family="Georgia"))

ggsave("2000_Robins_Hernan_and_Brumback/script_examples/point_treatment_confounding/table_w_weights.svg",
       height = 1.75, width=4)

data.frame(A = A, L = L, Y = Y) |> 
  group_by(A, L, Y) |> 
  tally() |> 
  left_join(ALw, by = c('A' = 'A', 'L' = 'L')) |> 
  mutate(n = n * w) |> 
  ggplot(aes(x = L, y = A, fill = n, label = scales::comma_format(accuracy = .1)(round(n, 1)))) + 
  geom_tile(alpha = 1, color = 'black', fill = 'white') + 
  geom_text(family = "Georgia") + 
  facet_grid(.~Y, labeller = labeller(Y = \(x) { paste0("Y = ", x) })) + 
  scale_y_continuous(breaks = c(0,1), labels = c("A = 0", "A = 1")) + 
  scale_x_continuous(breaks = c(0,1), labels = c("L = 0", "L = 1")) + 
  labs(x = "", y = "") + 
  ggtitle("Reweighted Population") + 
  theme_bw() + 
  theme(
    legend.position = 'none', panel.grid = element_blank(), 
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    text=element_text(family="Georgia"))

ggsave("2000_Robins_Hernan_and_Brumback/script_examples/point_treatment_confounding/reweighted_table.svg",
       height = 1.75, width=3.5)

