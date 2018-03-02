
library(ggplot2)
library(latex2exp)
library(mcclust)
library(mcclust.ext)
library(sdols)
library(plyr)
theme_set(theme_bw(base_size = 14))

setwd('Desktop/Discussion Paper/')

rm(list=ls())
source('functions.R')


# LOAD SARCOMA DATA -------------------------------------------------------

Ni <-  c(28, 29, 29, 26, 20, 15, 5, 12, 13, 2)
succ <- c(6, 7, 3, 5, 3, 2, 1, 1, 0, 0)

n <- length(succ)

Y <- array(NA, dim = c(n, max(Ni)))
for (i in 1:n){
  Y[i,1:succ[i]] <- 1
  if (succ[i] != max(Ni))
    Y[i,(succ[i] + 1):Ni[i]] <- 0
}


# Markov Chain Parameters
Niter <- 50000
burnin <- 20000
thin <- 6
N_tr <- 200

hyperpars <- NULL
hyperpars$a_tau <- hyperpars$b_tau <- .1
hyperpars$tau_mu <- sqrt(1000)
hyperpars$a_alpha <- 10; hyperpars$b_alpha <- 0.5


set.seed(1234)
fit <- probit_DP(Y, Ni, hyperpars, N_tr, Niter, burnin, thin)

samp_size <- (Niter - burnin)/thin
  

# INSPECT POSTERIOR -------------------------------------------------------

# Traceplot of theta_1
ggplot() + 
  geom_line(aes(x = 1:samp_size, 
                y = fit$thetas[,1]), 
            col = 'gray') + 
  labs(x = 'Iterations', y = TeX("$\\theta_{1}$"), 
       title = TeX('Traceplot of $\\theta_{1}')) + 
  scale_color_brewer(name="", palette = "Set1")


# Number of groups
ggplot(data.frame(num = min(fit$M):max(fit$M), 
                  prob = as.numeric(table(fit$M))/samp_size), aes(num, prob)) + 
  geom_bar(stat = "identity", width = 0.3, fill = 2) + 
  labs(x = 'Number of groups', y = "Probability") + 
  scale_x_continuous(breaks = min(fit$M):max(fit$M)) + 
  scale_color_brewer(name="", palette = "Set1")


# Compute the individual specific proportions from the cdf of the corresponding
# individual specific thetas
p_i <- array(NA, dim = c(samp_size, n))
for (i in 1:n){
  p_i[,i] <- pnorm(fit$thetas[,i])
}

post_quant <- apply(p_i, 2, quantile, prob = c(0.05, 0.5, 0.95))

# Plot
caterplot.df <- data.frame(group = factor(1:n), 
                             lo = post_quant[1,],
                             md = post_quant[2,],
                             hi = post_quant[3,])

# Compute the proportion of successes within each subtype
emp_succ <- succ/Ni

ggplot(caterplot.df, aes(ymin = lo, ymax = hi, x = group)) + 
  geom_linerange(lwd = 1) + 
  xlab("Observed Successes") +
  ylab("Success Rates") +
  # ggtitle("90% Credible intervals for the success probabilities") +
  geom_point(aes(y = md, col = "Estimated Median"), size = 2) +
  geom_point(aes(y = emp_succ, col = "True proportion"), size = 2) +
  theme_bw() + 
  geom_hline(yintercept = 0.1, col = 'black', lty = 2) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip() + 
  ylim(0, 0.4) + 
  scale_x_discrete(labels=c("1" = "6/28", "2" = "7/29", "3" = "3/29", 
                            "4" = "5/26", "5" = "3/20", "6" = "2/15", 
                            "7" = "1/5", "8" = "1/12", "9" = "0/13", 
                            "10" = "0/2")) + 
  annotate("text", x = 1, y = 0.32, label= "Leiomyo", hjust = 0) + 
  annotate("text", x = 2, y = 0.32, label= "Lipo", hjust = 0) + 
  annotate("text", x = 3, y = 0.32, label= "MFH", hjust = 0) + 
  annotate("text", x = 4, y = 0.32, label= "Osteo", hjust = 0) + 
  annotate("text", x = 5, y = 0.32, label= "Synovial", hjust = 0) + 
  annotate("text", x = 6, y = 0.32, label= "Angio", hjust = 0) + 
  annotate("text", x = 7, y = 0.32, label= "MPNST", hjust = 0) + 
  annotate("text", x = 8, y = 0.32, label= "Fibro", hjust = 0) + 
  annotate("text", x = 9, y = 0.32, label= "Ewings", hjust = 0) + 
  annotate("text", x = 10, y = 0.32, label= "Rhabdo", hjust = 0) + 
  theme(legend.title=element_blank(), legend.text=element_text(size=11), 
        axis.text=element_text(size=11))



# OPTIMAL PARTITION -------------------------------------------------------


# (1) Posterior mode
opt_mode <- as.numeric(post_mode(fit$K)[1,1:n])
opt_mode


# (2) Optimal according to Dahl (2006): L2 distance
P <- expectedPairwiseAllocationMatrix(fit$K)
opt_Dahl <- salso(P, loss = c("squaredError"))
length(opt_Dahl)


# (3) Optimal according to Wade (2017): Variation of information
opt_VI <- minVI(comp.psm(fit$K), fit$K, method=("all"), include.greedy=TRUE)
summary(opt_VI)


# (4) Optimal according to our method: repulsive loss function
lambda <- 100
b <- diff(range(fit$theta))/sqrt(10)
tau1sq <- 1

repuls <- min_repulsive_loss(fit$K, fit$theta, lambda, b, tau1sq)
  
opt_part <- repuls$opt_part
opt_params <- repuls$opt_params


ggplot(caterplot.df, aes(ymin = lo, ymax = hi, x = group)) + 
  geom_linerange(lwd = 1) + 
  xlab("Observed Successes") +
  ylab("Success Rates") +
  ggtitle("90% Credible intervals for the success probabilities") +
  geom_point(aes(y = md), col = 1, size = 2) +
  theme_bw() + 
  geom_hline(yintercept = 0.1, col = 'black', lty = 2) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip() + 
  ylim(0, 0.4) + 
  scale_x_discrete(labels=c("1" = "6/28", "2" = "7/29", "3" = "3/29", 
                            "4" = "5/26", "5" = "3/20", "6" = "2/15", 
                            "7" = "1/5", "8" = "1/12", "9" = "0/13", 
                            "10" = "0/2")) + 
  geom_point(aes(y = pnorm(opt_params)), col = opt_part + 1, 
             size = 2) +
  theme(legend.title=element_blank())

