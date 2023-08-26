
# Let A = treatment with a drug
# Let L = a predictor of treatment and a predictor of outcome Y 
# Let Y = some successful outcome induced by drug A and covariate L 
# and N is the size of our study population 

# data generation: 
N <- 10000
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

# test if claim that Pr(Y=1|A=a) remains unchanged ------------------------

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

data.frame(A = A, L = L, Y = Y) |> 
  group_by(A, L, Y) |> 
  tally() |> 
  left_join(ALw, by = c('A' = 'A', 'L' = 'L')) |> 
  ggplot(aes(x = L, y = A, label = paste0(scales::comma_format()(n), " * ", round(w, 1)))) + 
  geom_tile(alpha = 1, color = 'black', fill = 'white') + 
  geom_text() + 
  facet_grid(.~Y, labeller = labeller(Y = \(x) { paste0("Y = ", x) })) + 
  scale_y_continuous(breaks = c(0,1), labels = c("A = 0", "A = 1")) + 
  scale_x_continuous(breaks = c(0,1), labels = c("L = 0", "L = 1")) + 
  labs(x = "", y = "") + 
  ggtitle("Original Population") + 
  theme_bw() + 
  theme(
    legend.position = 'none', panel.grid = element_blank(), 
    axis.text.x = element_text(),
    axis.text.y = element_text())


data.frame(A = A, L = L, Y = Y, w = 1/probs_of_treatment) |> 
  group_by(A, L, Y) |> 
  summarize(n = sum(w)) |> 
  ggplot(aes(x = A, y = L, fill = n, label = scales::comma_format(accuracy = .1)(round(n, 1)))) + 
  geom_tile(alpha = 1, color = 'black', fill = 'white') + 
  geom_text() + 
  facet_grid(.~Y, labeller = labeller(Y = \(x) { paste0("Y = ", x) })) + 
  scale_x_continuous(breaks = c(0,1), labels = c("A = 0", "A = 1")) + 
  scale_y_continuous(breaks = c(0,1), labels = c("L = 0", "L = 1")) + 
  labs(x = "", y = "") + 
  ggtitle("Reweighted Population") + 
  theme_bw() + 
  theme(
    legend.position = 'none', panel.grid = element_blank(), 
    axis.text.x = element_text(),
    axis.text.y = element_text())

