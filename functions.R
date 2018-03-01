

# AUXILIARY FUNCTIONS -----------------------------------------------------


dnormix <- function(x, mu, sigma, w){
  # Computes the likelihood for the 1D mixture of K normals with same variance
  # --------------------------------------------------------------------------
  # Args:
  #   - x: point where to evaluate the likelihood
  #   - mu: vector of the means
  #   - sigma: scalar for the common standard deviation
  #   - w: vector of the weights
  # Returns:
  #   - out: value of the likelihood
  # ------------------------------------
  if (length(mu) != length(sigma)){
    stop("Error: the number of means and the number of variances should be 
         the same")
  }
  K = length(mu)
  out = rep(0, length(x))
  for (j in 1:K){
    out <- out + w[j] * dnorm(x, mu[j], sigma[j])
  }
  return (out)
}


rtrnorm <- function(n, mean, sd, a, b){
  # Sample from a truncated normal distribution on [a, b]
  # -----------------------------------------------------
  # Args:
  #   - n: number of samples needed
  #   - mean: mean of the normal distribution
  #   - sd: sd of the normal distribution
  #   - a: left bound of the truncation domain
  #   - b: right bound of the truncation domain
  # Returns:
  #   - x: sample from the truncated normal distribution
  # ----------------------------------------------------
  Fa <- pnorm(a, mean, sd)
  Fb <- pnorm(b, mean, sd)
  
  u <- runif(n, 0, 1)
  x <- qnorm(u * (Fb - Fa) + Fa, mean, sd)
  
  return(x)
}


gaps_removal <- function(mu, K){
  # Function to reorder the cluster indicators in a set 1:m, where m is the
  # number of groups
  # ----------------
  # Args:
  #   - mu: vector of group specific parameters
  #   - K: vector of cluster membership
  # Returns:
  #   - mu, sig, K ordered
  # -----------------
  N_tr <- length(mu)
  nj <- as.numeric(table(factor(K, levels = 1:N_tr)))
  empty <- which(nj == 0)
  nonempty <- which(nj != 0)
  mu_ord <- mu
  if (length(empty) > 0){
    for (k in nonempty){
      if (length(which(k > empty)) > 0){
        sub <- which(k > empty)[length(which(k > empty))]
        K[K == k] <- K[K == k] - sub
      }
    }
  }
  mu_ord[1:length(nonempty)] <- mu[nonempty]
  mu_ord[(length(nonempty) + 1):N_tr] <- mu[empty]
  return (list('mu_ord' = mu_ord, 'K_ord' = K))
}


gaps_removal_2 <- function(mu, sig, K){
  # Function to reorder the cluster indicators in a set 1:m, where m is the
  # number of groups
  # ----------------
  # Args:
  #   - mu: vector of group specific parameters
  #   - K: vector of cluster membership
  # Returns:
  #   - mu, sig, K ordered
  # -----------------
  N_tr <- length(mu)
  nj <- as.numeric(table(factor(K, levels = 1:N_tr)))
  empty <- which(nj == 0)
  nonempty <- which(nj != 0)
  mu_ord <- mu
  sig_ord <- sig
  if (length(empty) > 0){
    for (k in nonempty){
      if (length(which(k > empty)) > 0){
        sub <- which(k > empty)[length(which(k > empty))]
        K[K == k] <- K[K == k] - sub
      }
    }
  }
  mu_ord[1:length(nonempty)] <- mu[nonempty]
  mu_ord[(length(nonempty) + 1):N_tr] <- mu[empty]
  sig_ord[1:length(nonempty)] <- sig[nonempty]
  sig_ord[(length(nonempty) + 1):N_tr] <- sig[empty]
  return (list('mu_ord' = mu_ord, 'sig_ord' = sig_ord, 'K_ord' = K))
}



# BLOCKED GIBBS SAMPLER --------------------------------------------------


blocked_DP_fixvar <- function(X, xgrid, hypers, N_tr, Niter, burnin, thin){
  # Blocked Gibbs Sampler for a Dirichlet Process Mixture model 
  # as described in Ishwaran, James (2001) in 5.2.
  # Here only the means are group-dependent, whereas the variance is known
  # ----------------------------------------------------------------------
  # Args:
  #   - X: univariate data
  #   - xgrid: grid where to evaluate the predictive distribution
  #   - hypers: hyperparameters
  #   - N_tr: truncation level
  #   - Niter, burnin, thin: parameters of the Markov Chain
  # Returns:
  #   - MCMC output
  # ---------------
  
  n <- length(X)
  n_gr <- length(xgrid)
  samp_size <- (Niter - burnin)/thin
  
  
  # Define the data structures
  P <- array(NA, dim = c(n, N_tr))
  M <- array(NA, dim = samp_size)
  K_chain <- array(NA, dim = c(samp_size, n))
  theta_chain <- array(NA, dim = c(samp_size, n))
  pi_chain <- array(NA, dim = c(samp_size, N_tr))
  alpha_chain <- array(NA, dim = samp_size)
  pred_dens <- array(NA, dim = c(samp_size, n_gr))
  fit_dens <- array(NA, dim = c(samp_size, n))
  
  # Initialize the chain
  K_old <- sample(1:5, n, replace = T) # initialize the labels
  theta_old <- (rnorm(N_tr, hyperpars$mu_th, hyperpars$tau_th))
  pi_old <- runif(N_tr)
  pi_old <- pi_old/sum(pi_old)
  alpha_old <- 1
  
  
  # Auxiliary quantities
  xmat <- matrix(rep(xgrid, each=N_tr), n_gr, N_tr, byrow=T)
  xmat_fit <- matrix(rep(X, each=N_tr), n, N_tr, byrow=T)
  m <- length(unique(K_old)) # number of groups
  ni <- as.numeric(table(factor(K_old, levels = 1:N_tr))) # cluster sizes
  it <- 1
  
  # Gibbs Sampler
  pb <- txtProgressBar(style = 3)
  for (iter in 2:Niter){
    
    # --------------------------------------
    # (1) Cluster-specific parameters: means
    # --------------------------------------
    # group sufficient statistics
    Xsum <- aggregate(X, list(K_old), sum)$x # group sufficient statistics
    
    # Posterior updates
    var_post <- hypers$tau_th^2 / (ni[1:m] * hypers$tau_th^2 + 1)
    mean_post <- (hypers$tau_th^2 * Xsum + hypers$mu_th) / 
      (ni[1:m] * hypers$tau_th^2 + 1)
    
    # Draw the existing clusters' parameters
    theta_old[1:m] <- rnorm(m, mean = mean_post, sd = sqrt(var_post))
    # Draws from the prior
    theta_old[(m+1):N_tr] <- rnorm(N_tr-m, mean = hypers$mu_th, 
                                   sd = hypers$tau_th)
    
    
    # ---------------------
    # (2) Individual labels
    # ---------------------
    # K_i | rest ~ \sum_{k=1}^N.tr p_{k,i} \delta_k where
    # (p_{1,i}, ..., p_{N.tr,i}) \propto 
    #    (p_1 f(X_i | \mu_1), ..., p_N f(X_i | \mu_N))
    for (i in 1:n){
      P[i,] <- log(pi_old) + dnorm(X[i], theta_old, 1, log = T)
      K_old[i] <- sample(1:N_tr, 1, prob = exp(P[i,]))
    }
    ord <- gaps_removal(theta_old, K_old)
    theta_old <- ord$mu_ord 
    K_old <- ord$K_ord
    m <- length(unique(K_old))
    ni <- as.numeric(table(factor(K_old, levels = 1:N_tr)))
    
    
    # --------------------------
    # (3) Stick breaking weights
    # --------------------------
    # V_k \sim Beta(1 + M_k, \alpha + \sum_{l=k+1}^N.tr M_l)
    res <- n - cumsum(ni)
    V <- c(rbeta(N_tr - 1, 1 + ni[1:(N_tr-1)], alpha_old + res[1:(N_tr-1)]), 1)
    res_stick  <- c(0, cumsum(log(1 - V[1:(N_tr-1)])))
    log_pi_old <- log(V) + res_stick
    pi_old <- exp(log_pi_old)
    
    if ( (alpha_old/(alpha_old+1))^(N_tr-1) > 1E-2)
      stop('Incorrect Approximation: need larger N_tr')
    
    
    # ------------------------
    # (4) Total mass parameter
    # ------------------------
    alpha_old <- rgamma(1, hyperpars$a_alpha + N_tr - 1,
                        hyperpars$b_alpha - sum(log(1 + 1E-6 - V[1:(N_tr-1)])))
    
    
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # Save parameters in the chain
      M[it] <- m
      K_chain[it,] <- K_old
      pi_chain[it,] <- pi_old
      alpha_chain[it] <- alpha_old
      theta_chain[it,] <- theta_old[K_old]

      pars <- matrix(rep(theta_old, each=n_gr), n_gr, N_tr, byrow=F)
      pars_2 <- matrix(rep(theta_old, each=n), n, N_tr, byrow=F)
      pred_dens[it,] <- colSums(pi_chain[it,] * t(dnorm(xmat, pars, 1)))
      fit_dens[it,] <- colSums(pi_chain[it,] * t(dnorm(xmat_fit, pars_2, 1)))
      
      it <- it + 1
    }

  setTxtProgressBar(pb, iter/Niter)
  }
  
  
  return(list('K' = K_chain, 'theta' = theta_chain, 'pi' = pi_chain,
              'pred' = pred_dens, 'fit' = fit_dens, 'M' = M, 'alpha' = alpha_chain))
}


blocked_DP <- function(X, xgrid, hypers, N_tr, Niter, burnin, thin){
  # Blocked Gibbs Sampler for a Dirichlet Process Mixture model 
  # as described in Ishwaran, James (2001) in 5.2.
  # Here both the means and the standard deviations are group-dependent
  # (in the DP)
  # -----------
  # Args:
  #   - X: univariate data
  #   - xgrid: grid where to evaluate the predictive distribution
  #   - hypers: hyperparameters
  #   - N_tr: truncation level
  #   - Niter, burnin, thin: parameters of the Markov Chain
  # Returns:
  #   - MCMC output
  # ---------------
  
  n <- length(X)
  n_gr <- length(xgrid)
  samp_size <- (Niter - burnin)/thin
  
  
  # Define the data structures
  P <- array(NA, dim = c(n, N_tr))
  M <- array(NA, dim = samp_size)
  K_chain <- array(NA, dim = c(samp_size, n))
  theta_chain <- array(NA, dim = c(samp_size, n))
  sigma_chain <- array(NA, dim = c(samp_size, n))
  pi_chain <- array(NA, dim = c(samp_size, N_tr))
  alpha_chain <- array(NA, dim = samp_size)
  pred_dens <- array(NA, dim = c(samp_size, n_gr))
  fit_dens <- array(NA, dim = c(samp_size, n))
  
  # Initialize the chain
  K_old <- sample(1:5, n, replace = T) # initialize the labels
  theta_old <- (rnorm(N_tr, hyperpars$mu_th, hyperpars$tau_th))
  sigma_old <- (rnorm(N_tr, hyperpars$a_sig, hyperpars$b_sig))
  pi_old <- runif(N_tr)
  pi_old <- pi_old/sum(pi_old)
  alpha_old <- 1
  
  
  # Auxiliary quantities
  xmat <- matrix(rep(xgrid, each=N_tr), n_gr, N_tr, byrow=T)
  xmat_fit <- matrix(rep(X, each=N_tr), n, N_tr, byrow=T)
  m <- length(unique(K_old)) # number of groups
  ni <- as.numeric(table(factor(K_old, levels = 1:N_tr))) # cluster sizes
  it <- 1
  
  # Gibbs Sampler
  pb <- txtProgressBar(style = 3)
  for (iter in 2:Niter){
    
    # --------------------------------------
    # (1) Cluster-specific parameters: means
    # --------------------------------------
    # group sufficient statistics
    Xsum <- aggregate(X, list(K_old), sum)$x # group sufficient statistics
    
    # Posterior updates
    var_post <- sigma_old[1:m]^2 * hypers$tau_th^2 / 
      (ni[1:m] * hypers$tau_th^2 + 1)
    mean_post <- (hypers$tau_th^2 * Xsum + hypers$mu_th) / 
      (ni[1:m] * hypers$tau_th^2 + 1)
    
    # Draw the existing clusters' parameters
    theta_old[1:m] <- rnorm(m, mean = mean_post, sd = sqrt(var_post))
    # Draws from the prior
    theta_old[(m+1):N_tr] <- rnorm(N_tr-m, mean = hypers$mu_th, 
                                   sd = hypers$tau_th)
    
    
    # ------------------------------------------
    # (2) Cluster-specific parameters: variances
    # ------------------------------------------
    RSS1 <- 0.5 * (aggregate(X^2, list(K_old), sum)$x - 
                     2 * aggregate(X, list(K_old), sum)$x * theta_old[1:m] + 
                     ni[1:m] * theta_old[1:m]^2)
    
    RSS2 <- 0.5/hyperpars$tau_th^2 * (theta_old[1:m] - hyperpars$mu_th)^2
    
    # Draw the existing clusters' parameters
    sigma_old[1:m] <- 1 / sqrt(rgamma(m, hyperpars$a_sig + 0.5 + 0.5*ni[1:m], 
                                      hyperpars$b_sig + RSS1 + RSS2))
    # Draws from the prior
    sigma_old[(m+1):N_tr] <- 1 / 
      sqrt(rgamma(N_tr-m, hyperpars$a_sig, hyperpars$b_sig))
    

    # ---------------------
    # (3) Individual labels
    # ---------------------
    # K_i | rest ~ \sum_{k=1}^N.tr p_{k,i} \delta_k where
    # (p_{1,i}, ..., p_{N.tr,i}) \propto 
    #    (p_1 f(X_i | \mu_1), ..., p_N f(X_i | \mu_N))
    for (i in 1:n){
      P[i,] <- log(pi_old) + dnorm(X[i], theta_old, sigma_old, log = T)
      K_old[i] <- sample(1:N_tr, 1, prob = exp(P[i,]))
    }
    ord <- gaps_removal_2(theta_old, sigma_old, K_old)
    theta_old <- ord$mu_ord
    sigma_old <- ord$sig_ord
    K_old <- ord$K_ord
    m <- length(unique(K_old))
    ni <- as.numeric(table(factor(K_old, levels = 1:N_tr)))
    
    
    # --------------------------
    # (4) Stick breaking weights
    # --------------------------
    # V_k \sim Beta(1 + M_k, \alpha + \sum_{l=k+1}^N.tr M_l)
    res <- n - cumsum(ni)
    V <- c(rbeta(N_tr - 1, 1 + ni[1:(N_tr-1)], alpha_old + res[1:(N_tr-1)]), 1)
    res_stick  <- c(0, cumsum(log(1 - V[1:(N_tr-1)])))
    log_pi_old <- log(V) + res_stick
    pi_old <- exp(log_pi_old)
    
    if ( (alpha_old/(alpha_old+1))^(N_tr-1) > 1E-2)
      stop('Incorrect Approximation: need larger N_tr')
    
    
    # ------------------------
    # (5) Total mass parameter
    # ------------------------
    alpha_old <- rgamma(1, hyperpars$a_alpha + N_tr - 1,
                        hyperpars$b_alpha - sum(log(1 + 1E-6 - V[1:(N_tr-1)])))
    
    
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # Save parameters in the chain
      M[it] <- m
      K_chain[it,] <- K_old
      pi_chain[it,] <- pi_old
      alpha_chain[it] <- alpha_old
      theta_chain[it,] <- theta_old[K_old]
      sigma_chain[it,] <- sigma_old[K_old]
      
      pars <- matrix(rep(theta_old, each=n_gr), n_gr, N_tr, byrow=F)
      pars_2 <- matrix(rep(theta_old, each=n), n, N_tr, byrow=F)
      pars_var <- matrix(rep(sigma_old, each=n_gr), n_gr, N_tr, byrow=F)
      pars_var_2 <- matrix(rep(sigma_old, each=n), n, N_tr, byrow=F)
      pred_dens[it,] <- colSums(pi_chain[it,] * t(dnorm(xmat, pars, pars_var)))
      fit_dens[it,] <- colSums(pi_chain[it,] * t(dnorm(xmat_fit, pars_2, 
                                                       pars_var_2)))
      
      it <- it + 1
    }
    
    setTxtProgressBar(pb, iter/Niter)
  }
  
  
  return(list('K' = K_chain, 'theta' = theta_chain, 'sigma' = sigma_chain, 
              'pi' = pi_chain, 'pred' = pred_dens, 'fit' = fit_dens, 'M' = M, 
              'alpha' = alpha_chain))
}


blocked_DP_2 <- function(X, xgrid, hypers, N_tr, Niter, burnin, thin){
  # Blocked Gibbs Sampler for a Dirichlet Process Mixture model 
  # as described in Ishwaran, James (2001) in 5.2.
  # Here both the means and the standard deviations are group-dependent
  # (in the DP)
  # -----------
  # Args:
  #   - X: univariate data
  #   - xgrid: grid where to evaluate the predictive distribution
  #   - hypers: hyperparameters
  #   - N_tr: truncation level
  #   - Niter, burnin, thin: parameters of the Markov Chain
  # Returns:
  #   - MCMC output
  # ---------------
  
  n <- length(X)
  n_gr <- length(xgrid)
  samp_size <- (Niter - burnin)/thin
  
  
  # Define the data structures
  P <- array(NA, dim = c(n, N_tr))
  M <- array(NA, dim = samp_size)
  K_chain <- array(NA, dim = c(samp_size, n))
  theta_chain <- array(NA, dim = c(samp_size, n))
  sigma_chain <- array(NA, dim = c(samp_size, n))
  pi_chain <- array(NA, dim = c(samp_size, N_tr))
  alpha_chain <- array(NA, dim = samp_size)
  mu_th_chain <- array(NA, dim = samp_size)
  tau_th_chain <- array(NA, dim = samp_size)
  b_sig_chain <- array(NA, dim = samp_size)
  
  pred_dens <- array(NA, dim = c(samp_size, n_gr))
  fit_dens <- array(NA, dim = c(samp_size, n))
  
  # Initialize the chain
  K_old <- sample(1:5, n, replace = T) # initialize the labels
  pi_old <- runif(N_tr)
  pi_old <- pi_old/sum(pi_old)
  alpha_old <- 1
  mu_th_old <- 0
  tau_th_old <- 10
  b_sig_old <- 1
  theta_old <- (rnorm(N_tr, mu_th_old, tau_th_old))
  sigma_old <- (rnorm(N_tr, hyperpars$a_sig, b_sig_old))
  
  # Auxiliary quantities
  xmat <- matrix(rep(xgrid, each=N_tr), n_gr, N_tr, byrow=T)
  xmat_fit <- matrix(rep(X, each=N_tr), n, N_tr, byrow=T)
  RSS1 <- array(NA, dim = N_tr)
  m <- length(unique(K_old)) # number of groups
  ni <- as.numeric(table(factor(K_old, levels = 1:N_tr))) # cluster sizes
  it <- 1
  
  # Gibbs Sampler
  pb <- txtProgressBar(style = 3)
  for (iter in 2:Niter){
    
    # --------------------------------------
    # (1) Cluster-specific parameters: means
    # --------------------------------------
    # group sufficient statistics
    Xsum <- aggregate(X, list(K_old), sum)$x # group sufficient statistics
    
    # Posterior updates
    var_post <- sigma_old[1:m]^2 * tau_th_old^2 / 
      (ni[1:m] * tau_th_old^2 + 1)
    mean_post <- (tau_th_old^2 * Xsum + mu_th_old) / 
      (ni[1:m] * tau_th_old^2 + 1)
    
    # Draw the existing clusters' parameters
    theta_old[1:m] <- rnorm(m, mean = mean_post, sd = sqrt(var_post))
    # Draws from the prior
    theta_old[(m+1):N_tr] <- rnorm(N_tr-m, mean = mu_th_old, 
                                   sd = tau_th_old)
    
    
    # ------------------------------------------
    # (2) Cluster-specific parameters: variances
    # ------------------------------------------
    for (k in 1:m){
      RSS1[k] <- 0.5 * sum((X[K_old == k] - theta_old[k])^2)
    }

    RSS2 <- 0.5/tau_th_old^2 * (theta_old[1:m] - mu_th_old)^2

    
    # Draw the existing clusters' parameters
    sigma_old[1:m] <- 1 / sqrt(rgamma(m, hyperpars$a_sig + 0.5 + 0.5*ni[1:m],
                                      b_sig_old + RSS1[1:m] + RSS2))
    
    # Draws from the prior
    sigma_old[(m+1):N_tr] <- 1 / 
      sqrt(rgamma(N_tr-m, hyperpars$a_sig, b_sig_old))

    
    # ---------------------------------------
    # (3) Hyperparameters of the base measure
    # ---------------------------------------
    sum_prec <- sum(1/sigma_old[1:m]^2)
    den <- hypers$s1^2 * sum_prec + tau_th_old^2
    mu_th_old <- rnorm(1, (hypers$s1^2 * sum(theta_old[1:m]/sigma_old[1:m]^2) +
                             tau_th_old^2 * hypers$m1)/den,
                       sqrt(tau_th_old^2 * hypers$s1^2/den))
    tau_th_old <- 1/sqrt(rgamma(1, hypers$tau1 + m/2, hypers$tau2 + 0.5 *
                                  sum(((theta_old[1:m] - mu_th_old)/sigma_old[1:m])^2)))
    b_sig_old <- rgamma(1, hypers$nu1 + m*hypers$a_sig, hypers$nu2 + sum_prec)
    
    
    # ---------------------
    # (4) Individual labels
    # ---------------------
    # K_i | rest ~ \sum_{k=1}^N.tr p_{k,i} \delta_k where
    # (p_{1,i}, ..., p_{N.tr,i}) \propto 
    #    (p_1 f(X_i | \mu_1), ..., p_N f(X_i | \mu_N))
    for (i in 1:n){
      P[i,] <- log(pi_old) + dnorm(X[i], theta_old, sigma_old, log = T)
      K_old[i] <- sample(1:N_tr, 1, prob = exp(P[i,]))
    }
    ord <- gaps_removal_2(theta_old, sigma_old, K_old)
    theta_old <- ord$mu_ord
    sigma_old <- ord$sig_ord
    K_old <- ord$K_ord
    m <- length(unique(K_old))
    ni <- as.numeric(table(factor(K_old, levels = 1:N_tr)))
    
    
    # --------------------------
    # (5) Stick breaking weights
    # --------------------------
    # V_k \sim Beta(1 + M_k, \alpha + \sum_{l=k+1}^N.tr M_l)
    res <- n - cumsum(ni)
    V <- c(rbeta(N_tr - 1, 1 + ni[1:(N_tr-1)], alpha_old + res[1:(N_tr-1)]), 1)
    res_stick  <- c(0, cumsum(log(1 - V[1:(N_tr-1)])))
    log_pi_old <- log(V) + res_stick
    pi_old <- exp(log_pi_old)
    
    if ( (alpha_old/(alpha_old+1))^(N_tr-1) > 1E-2)
      stop('Incorrect Approximation: need larger N_tr')
    
    
    # ------------------------
    # (6) Total mass parameter
    # ------------------------
    alpha_old <- rgamma(1, hyperpars$a_alpha + N_tr - 1,
                        hyperpars$b_alpha - sum(log(1 + 1E-6 - V[1:(N_tr-1)])))
    
    
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # Save parameters in the chain
      M[it] <- m
      K_chain[it,] <- K_old
      pi_chain[it,] <- pi_old
      alpha_chain[it] <- alpha_old
      theta_chain[it,] <- theta_old[K_old]
      sigma_chain[it,] <- sigma_old[K_old]
      mu_th_chain[it] <- mu_th_old
      tau_th_chain[it] <- tau_th_old
      b_sig_chain[it] <- b_sig_old
      
      pars <- matrix(rep(theta_old, each=n_gr), n_gr, N_tr, byrow=F)
      pars_2 <- matrix(rep(theta_old, each=n), n, N_tr, byrow=F)
      pars_var <- matrix(rep(sigma_old, each=n_gr), n_gr, N_tr, byrow=F)
      pars_var_2 <- matrix(rep(sigma_old, each=n), n, N_tr, byrow=F)
      pred_dens[it,] <- colSums(pi_chain[it,] * t(dnorm(xmat, pars, pars_var)))
      fit_dens[it,] <- colSums(pi_chain[it,] * t(dnorm(xmat_fit, pars_2, 
                                                       pars_var_2)))
      
      it <- it + 1
    }
    
    setTxtProgressBar(pb, iter/Niter)
  }
  
  
  return(list('K' = K_chain, 'theta' = theta_chain, 'sigma' = sigma_chain, 
              'pi' = pi_chain, 'pred' = pred_dens, 'fit' = fit_dens, 'M' = M, 
              'alpha' = alpha_chain, 'mu_th' = mu_th_chain, 
              'tau_th' = tau_th_chain, 'b_sig' = b_sig_chain))
}


probit_DP <- function(Y, Ni, hyperpars, N_tr, Niter, burnin, thin){
  # Blocked Gibbs Sampler for a Binomial-Probit Dirichlet Process Mixture model 
  # where the normal inverse cdf of the probabilities of success are modelled 
  # as a nonparametric mixture model
  # --------------------------------
  # Args:
  #   - Y: number of successes of the Binomial trials
  #   - Ni: number of Binomial trials
  #   - hypers: hyperparameters
  #   - N_tr: truncation level
  #   - Niter, burnin, thin: parameters of the Markov Chain
  # Returns:
  #   - MCMC output
  # ---------------
  
  n <- nrow(Y)
  N <- max(Ni)
  samp_size <- (Niter - burnin)/thin
  

  # Define the data structures
  P <- array(NA, dim = c(n, N_tr))
  M <- array(NA, dim = samp_size)
  K_chain <- array(NA, dim = c(samp_size, n))
  thetas_chain <-  array(NA, dim = c(samp_size, n))
  mu_chain <- array(NA, dim = samp_size)
  tau_chain <- array(NA, dim = samp_size)
  pi_chain <- array(NA, dim = c(samp_size, N_tr))
  alpha_chain <- array(NA, dim = samp_size)
  
  
  # Initialize the chain
  K_old <- sample(1:2, n, replace = T) # initialize the labels
  thetas_old <- rep(0, n)
  mu_old <- 0
  tau_old <- 10
  zetas_old <- array(NA, dim = c(n, N))
  pi_old <- runif(N_tr)
  pi_old <- pi_old/sum(pi_old)
  alpha_old <- 1


  # Auxiliary quantities
  m <- M_old <- length(unique(K_old)) # number of groups
  ni <- as.numeric(table(factor(K_old, levels = 1:N_tr))) # cluster sizes
  it <- 1
  
  # Gibbs Sampler
  pb <- txtProgressBar(style = 3)
  for (iter in 2:Niter){
    
    # --------------------------------------
    # (1) Imputation of the latent variables
    # --------------------------------------
    for (i in 1:n){
      for (j in 1:Ni[i]){
        if (Y[i,j] == 1){
          zetas_old[i,j] <- rtrnorm(1, thetas_old[K_old[i]], 1, 0, +Inf)
        }
        else{
          zetas_old[i,j] <- rtrnorm(1, thetas_old[K_old[i]], 1, -Inf, 0)
        }
      }
    }
    
    
    # --------------------------------------
    # (2) Cluster-specific parameters: means
    # --------------------------------------
    # Posterior updates
    var_post <- mean_post <- array(NA, dim = m)
    for (k in 1:m){
      var_post[k] <- tau_old^2/
        (tau_old^2 * sum(Ni[K_old == k]) + 1)
      mean_post[k] <- (mu_old + tau_old^2 * 
                         sum(zetas_old[K_old == k,], na.rm = T)) / 
        (tau_old^2 * sum(Ni[K_old == k]) + 1)
    }
    
    # Draw the existing clusters' parameters
    thetas_old[1:m] <- rnorm(m, mean = mean_post, sd = sqrt(var_post))
    
    # Draws from the prior
    thetas_old[(m+1):N_tr] <- rnorm(N_tr - m, mean = mu_old, 
                                sd = tau_old)
    
    
    # ---------------------
    # (3) Individual labels
    # ---------------------
    # K_i | rest ~ \sum_{k=1}^N.tr p_{k,i} \delta_k where
    # (p_{1,i}, ..., p_{N.tr,i}) \propto 
    #    (p_1 f(X_i | \mu_1), ..., p_N f(X_i | \mu_N))
    for (i in 1:n){
      for (k in 1:N_tr){
        P[i,k] <- log(pi_old[k]) +
          sum(dnorm(zetas_old[i, 1:Ni[i]], thetas_old[k], 1, log = T))
      }
      K_old[i] <- sample(1:N_tr, 1, prob = exp(P[i,]))
    }
    ord <- gaps_removal(thetas_old, K_old)
    thetas_old <- ord$mu_ord
    K_old <- ord$K_ord
    m <- M_old <- length(unique(K_old))
    ni <- as.numeric(table(factor(K_old, levels = 1:N_tr)))
    
    
    # --------------------------
    # (4) Stick breaking weights
    # --------------------------
    # V_k \sim Beta(1 + M_k, \alpha + \sum_{l=k+1}^N.tr M_l)
    res <- n - cumsum(ni)
    V <- c(rbeta(N_tr - 1, 1 + ni[1:(N_tr-1)], alpha_old + res[1:(N_tr-1)]), 1)
    res_stick  <- c(0, cumsum(log(1 - V[1:(N_tr-1)])))
    log_pi_old <- log(V) + res_stick
    pi_old <- exp(log_pi_old)
    
    if ( (alpha_old/(alpha_old+1))^(N_tr-1) > 1E-2)
      stop('Incorrect Approximation: need larger N_tr')
    
    
    # ------------------------
    # (5) Total mass parameter
    # ------------------------
    alpha_old <- rgamma(1, hyperpars$a_alpha + N_tr - 1,
                        hyperpars$b_alpha - sum(log(1 + 1E-6 - V[1:(N_tr-1)])))
    #alpha_old <- 10
    
    # -------------------
    # (6) Hyperparameters
    # -------------------
    mean_post <- hyperpars$tau_mu^2 * sum(thetas_old) /
      (N_tr * hyperpars$tau_mu^2 + tau_old^2)
    var_post <- (hyperpars$tau_mu^2 * tau_old^2) /
      (N_tr * hyperpars$tau_mu^2 + tau_old^2)
    mu_old <- rnorm(1, mean_post, sqrt(var_post))

    RSS <- sum((thetas_old - mu_old)^2)
    a_post <- hyperpars$a_tau + N_tr/2
    b_post <- hyperpars$b_tau + 0.5 * RSS
    tau_old <- 1/sqrt(rgamma(1, a_post, b_post))
    
    
    if ( (iter > burnin) && (iter %% thin == 0) ){
      # Save parameters in the chain
      M[it] <- M_old
      K_chain[it,] <- K_old
      thetas_chain[it,] <- thetas_old[K_old]
      pi_chain[it,] <- pi_old
      alpha_chain[it] <- alpha_old
      mu_chain[it] <- mu_old
      tau_chain[it] <- tau_old

      it <- it + 1
    }
    
    setTxtProgressBar(pb, iter/Niter)
  }

  return(list('thetas' = thetas_chain, 'pi' = pi_chain, 'alpha' = alpha_chain,
              'mu' = mu_chain, 'tau' = tau_chain, 'K' = K_chain, 'M' = M))
}




# POSTPROCESSING FUNCTIONS ------------------------------------------------


post_mode <- function(s){
  # Function that computes the uniques rows of a matrix ordering them by their
  # frequencies
  # -----------
  # Args: 
  #   - s: matrix of cluster allocations
  # Returns: 
  #   - part_unique: uniques rows ordered by frequence
  # --------------------------------------------------
  n <- ncol(s)
  colnames(s) <- as.character(1:n)
  part_unique <- plyr::count(as.data.frame(s), vars = 1:n)
  part_unique <- part_unique[order(part_unique$freq,decreasing=TRUE),]
  
  return(part_unique)
}


post_prob <- function(s){
  # Function that computes the matrix of posterior probability of belonging 
  # to the same cluster
  # -------------------
  # Args: 
  #   - s: first vector of points
  # Returns: 
  #   - P: matrix of posterior probability of belonging to the same cluster
  # -----------------------------------------------------------------------
  samp_size <- nrow(s)
  n <- ncol(s)
  P <- matrix(NA, nrow = n, ncol = n)
  diag(P) <- rep(1,n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      P[i,j] <- P[j,i] <- sum(s[,i]==s[,j])/samp_size
    }
  }
  
  return(P)
}


adjacency_mat <- function(labels){
  # Function that computes the adjacency matrix, given a vector of 
  # cluster labels
  # --------------
  # Args: 
  #   - labels: vector of cluster allocation
  # Returns: 
  #   - P: incidence matrix
  # -----------------------
  n <- length(labels)
  P <- matrix(0, nrow = n, ncol = n)
  K <- length(unique(labels))
  
  diag(P) <- rep(1,n)

  for (k in 1:K){
    cl_k <- which(labels == k)
    for (i in 1:length(cl_k)){
      for (j in (i+1):length(cl_k)){
        P[cl_k[i], cl_k[j]] <- P[cl_k[j], cl_k[i]] <- 1
      }
    }
  }

  return(P)
}


my_dist <- function(x, y){
  # Function that computes the distance matrix D, where D_{ij} is x_i - y_j
  # -----------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points 
  #   - y: second vector of points 
  # Returns: 
  #   - D: matrix of pairwise distances
  # -----------------------------------
  n <- length(x)
  m <- length(y)
  
  D <- array(NA, dim = c(n, m))
  
  # In order to speed up the code, we vectorize the distance computations: 
  # we create two matrices, X and Y, by replicating the vectors x and y. 
  # The only operation we do is the difference between the two matrices.
  X <- matrix(rep(x, m), nrow = n, ncol = m, byrow = F)
  Y <- matrix(rep(y, n), nrow = n, ncol = m, byrow = T)
  
  D <- abs(X - Y)
  
  return (D)
}


C_squaredexp <- function(x, y, b, tau1sq){
  # Function that computes the squared exponential covariance function
  # ------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points
  #   - y: second vector of points
  #   - b, tau1sq: parameters of the SE covariance function
  # Returns: 
  #   - CSE: squared exponential covariance matrix
  # ----------------------------------------------
  n <- length(x)
  m <- length(y)
  CSE <- array(NA, dim = c(n, m))
  
  d <- my_dist(x, y)
  CSE <- tau1sq * exp(-(1/2)*(d/b)^2)
  
  return(CSE)
}


min_repulsive_loss <- function(s, theta, lambda, b, tau1sq){
  # Function that computes the incidence matrix, given a vector of 
  # cluster labels
  # --------------
  # Args: 
  #   - s: matrix of cluster allocations
  #   - theta: matrix of individual specific parameters
  #   - lambda: penalty to induce repulsion
  #   - b, tau1sq: parameters of the squared exponential
  # Returns: 
  #   - opt_part, opt_params: optimal partition and corresponding parameters
  # ------------------------------------------------------------------------
  n <- ncol(s)
  samp_size <- nrow(s)
  
  # Compute the posterior mean of the individual specific parameters (denoted 
  # by c_n in the draft)
  c_n <- array(NA, dim = n)
  for (i in 1:n){
    c_n[i] <- mean(theta[,i])
  }
  
  # Compute the loss function on the visited partitions
  l <- array(NA, dim = samp_size)
  pb <- txtProgressBar(style = 3)
  for (i in 1:samp_size){
    unique_pars <- unique(theta[i,])
    C_mat <- C_squaredexp(unique_pars, unique_pars, b, tau1sq)
    l[i] <- sum(theta[i,]^2 - 2 * c_n * theta[i,])/n - 
      lambda * length(unique_pars) * det(C_mat)
    setTxtProgressBar(pb, i/samp_size)
  }
  
  opt_part <- s[which.min(l),]
  opt_params <- theta[which.min(l),]
  
  return(list('opt_part' = opt_part, 'opt_params' = opt_params))
}



C_squaredexp_2D <- function(x, y, b, tau1sq){
  # Function that computes the squared exponential covariance function
  # ------------------------------------------------------------------
  # Args: 
  #   - x: first vector of points
  #   - y: second vector of points
  #   - b, tau1sq: parameters of the SE covariance function
  # Returns: 
  #   - CSE: squared exponential covariance matrix
  # ----------------------------------------------
  n <- nrow(x)
  m <- nrow(y)
  
  if (ncol(x) != ncol(y)){
    stop('')
  }
  p <- ncol(x)
  
  CSE <- array(1, dim = c(n, m))
  
  for (j in 1:p){
    d <- my_dist(x[,j], y[,j])
    CSE <- CSE * exp(-(1/2) * (d/b[j])^2)
  }

  CSE <- CSE * tau1sq
  
  return(CSE)
}


min_repulsive_loss_2D <- function(s, theta, sigma, lambda, b, tau1sq){
  # Function that computes the incidence matrix, given a vector of 
  # cluster labels
  # --------------
  # Args: 
  #   - s: matrix of cluster allocations
  #   - theta: matrix of individual specific parameters
  #   - lambda: penalty to induce repulsion
  #   - b, tau1sq: parameters of the squared exponential
  # Returns: 
  #   - opt_part, opt_params: optimal partition and corresponding parameters
  # ------------------------------------------------------------------------
  n <- ncol(s)
  samp_size <- nrow(s)
  p <- length(b)
  
  # Compute the posterior mean of the individual specific parameters (denoted 
  # by c_n in the draft)
  c_n <- array(NA, dim = c(n, p))
  for (i in 1:n){
    c_n[i,] <- colMeans(cbind(theta[,i], sigma[,i]))
  }
  
  # Compute the loss function on the visited partitions
  l <- array(NA, dim = samp_size)
  pb <- txtProgressBar(style = 3)
  for (i in 1:samp_size){
    pars <- cbind(theta[i,],sigma[i,])
    unique_pars <- cbind(unique(theta[i,]), unique(sigma[i,]))
    C_mat <- C_squaredexp_2D(unique_pars, unique_pars, b, tau1sq)
    l[i] <- sum(pars^2 - 2 * c_n * pars)/n - 
      lambda * nrow(unique_pars) * det(C_mat)
    setTxtProgressBar(pb, i/samp_size)
  }
  
  opt_part <- s[which.min(l),]
  opt_params <- cbind(theta[which.min(l),], sigma[which.min(l),])
  
  return(list('opt_part' = opt_part, 'opt_params' = opt_params))
}
