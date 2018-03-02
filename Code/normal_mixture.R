
library(ggplot2)
library(latex2exp)
library(mcclust)
library(mcclust.ext)
library(sdols)
library(plyr)
theme_set(theme_bw(base_size = 14))


rm(list=ls())

setwd('~/Desktop/Discussion Paper/')
source('functions.R')


# SIMULATE DATA -----------------------------------------------------------
set.seed(124)
n <- 1000 # number of data points
J <- 7    # number of components 

pi_true <- rep(1/J, J) # true proportions
mu_true <- c(-3, -3.5, -2.6, 0, 1.8, 2.4, 7.1) # true means
sd_true <- rep(1, J)


# We generate the unnoised X_i from a mixture of two normals
# The goal is to infer this unobserved density, by only observing the noisy
# proxies W_{ij}
cl <- sample(1:J, n, TRUE, prob = pi_true)
X <- array(NA, dim = n)
for (i in 1:n){
  X[i] <- rnorm(1, mu_true[cl[i]], sd_true[cl[i]])
}

xgrid <- seq(min(X) - 2, max(X) + 2, length.out = 100)

# Inspect the predictive density
ggplot() + 
  geom_histogram(aes(x = X, y = ..density..), bins = 40, fill = 'gray', col = 'white') + 
  geom_line(aes(x = xgrid, y = dnormix(xgrid, mu_true, sd_true, pi_true), 
                col = 'True'), size = 0.8) + 
  labs(x = 'X', y = 'Density', title = 'Histogram') + 
  scale_color_brewer(name="", palette = "Set1")

hyperpars <- NULL
hyperpars$mu_th <- 0; hyperpars$tau_th <- 200
hyperpars$a_sig <- 1; hyperpars$b_sig <- 1
hyperpars$a_alpha <- 2; hyperpars$b_alpha <- 2


# Markov Chain Parameters
Niter <- 12000
burnin <- 2000
thin <- 2
N_tr <- 40

fit <- blocked_DP(X, xgrid, hyperpars, N_tr, Niter, burnin, thin)

samp_size <- (Niter - burnin)/thin


# INSPECT POSTERIOR -------------------------------------------------------


pred_quant <- apply(fit$pred, 2, quantile, prob = c(0.025, 0.975))
# Inspect the predictive density
ggplot() + 
  geom_histogram(aes(x = X, y = ..density..), bins = 40, fill = 'gray', col = 'white') + 
  geom_ribbon(aes(x = xgrid, ymin = pred_quant[1,], ymax = pred_quant[2,]), 
              fill = 2, alpha = 0.3) + 
  geom_line(aes(x = xgrid, y = colMeans(fit$pred), col = 'Estimated'), 
            size = 0.8) + 
  geom_line(aes(x = xgrid, y = dnormix(xgrid, mu_true, sd_true, pi_true), 
                col = 'True'), size = 0.8) + 
  labs(x = 'X', y = 'Density', title = 'Histogram') + 
  scale_color_brewer(name="", palette = "Set1")


# Number of groups
ggplot(data.frame(num = min(fit$M):max(fit$M), prob = as.numeric(table(factor(fit$M, levels = min(fit$M):max(fit$M)))/samp_size)), aes(num, prob)) + 
  geom_bar(stat = "identity", width = 0.3, fill = 2) + 
  labs(x = 'Number of groups', y = "Probability") + 
  scale_x_continuous(breaks = min(fit$M):max(fit$M)) + 
  scale_color_brewer(name="", palette = "Set1")



# OPTIMAL PARTITION -------------------------------------------------------


# (1) Posterior mode:
post_clust <- post_mode(fit$K)
opt_mode <- as.numeric(post_clust[1,1:n])
opt_mode # overfit


# (2) Optimal according to Dahl (2006): L2 distance
P <- expectedPairwiseAllocationMatrix(fit$K)
opt_Dahl <- salso(P, loss = c("squaredError"))
opt_Dahl


# (3) Optimal according to Lau and Green (2007): Binder's loss
opt_Binder <- salso(P, loss = c("binder"))
opt_Binder


# (4) Optimal according to Wade (2017): Variation of information
opt_VI <- minVI(comp.psm(fit$K), fit$K, method=("all"), include.greedy=TRUE)
summary(opt_VI)
opt_VI$cl[1,]

# Using the sdols package we should get the same result
opt_VI_2 <- salso(P, loss = c("lowerBoundVariationOfInformation"))
opt_VI_2


# Plot the results obtained via VI loss
cols <- c('dodgerblue2','firebrick1','darkorchid1','gold', 'seagreen3')
ggplot() + 
  geom_histogram(aes(x = X, y = ..density..), bins = 40, fill = 'gray', col = 'white') + 
  geom_ribbon(aes(x = xgrid, ymin = pred_quant[1,], ymax = pred_quant[2,]), 
              fill = 1, alpha = 0.3) + 
  geom_line(aes(x = xgrid, y = colMeans(fit$pred)), col = 'black', 
            size = 0.8) + 
  geom_point(aes(x = X, y = colMeans(fit$fit)), col = cols[factor(opt_VI$cl[1,])], 
             pch = 1, size = 1.8) + 
  labs(x = 'X', y = '') + 
  scale_color_brewer(name="", palette = "Set1") +
  theme(axis.text=element_text(size=11))


# (5) Optimal according to our method: repulsive loss function
lambda <- 100
b1 <- diff(range(fit$theta))/sqrt(10)
tau1sq <- 1

repuls <- min_repulsive_loss(fit$K, fit$theta, lambda, b1, tau1sq)

opt_part <- repuls$opt_part
opt_params <- repuls$opt_params

length(unique(opt_part)) # number of groups in the optimal partition

cols <- c('gold','darkorchid1','firebrick1','dodgerblue2', 'seagreen2', 'black')
# Plot the results obtained via our loss
ggplot() + 
  geom_histogram(aes(x = X, y = ..density..), bins = 40, fill = 'gray', col = 'white') + 
  geom_ribbon(aes(x = xgrid, ymin = pred_quant[1,], ymax = pred_quant[2,]), 
              fill = 1, alpha = 0.3) + 
  geom_line(aes(x = xgrid, y = colMeans(fit$pred)), col = 'black', 
            size = 0.8) + 
  geom_point(aes(x = X, y = colMeans(fit$fit)), col = cols[factor(opt_part)], 
             pch = 1, size = 1.8) + 
  labs(x = 'X', y = '') + 
  scale_color_brewer(name="", palette = "Set1")

