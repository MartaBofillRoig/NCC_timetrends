######################################################################################################################################

data_sim_block <- function(K=1, mu0=0, delta, p0, OR, lambda, sigma, N1, N2, N3, N_peak, trend, trend_param, endpoint){
  N1 <- round(N1)
  N2 <- round(N2)
  N3 <- round(N3)
  
  SS_matrix <- matrix(c(c(N1, NA), N2, c(N3[1], NA, N3[2])), nrow = 3)
  
  alloc_ratios <- ifelse(!is.na(SS_matrix), 1, 0)
  
  num_periods <- ncol(alloc_ratios) # total number of periods
  num_arms <- nrow(alloc_ratios)-1 # total number of treatment arms
  
  N_period <- colSums(SS_matrix, na.rm=T) # sample sizes per period
  N_arm <- rowSums(SS_matrix, na.rm=T) # sample sizes per arm
  n_total <- sum(SS_matrix, na.rm = T) # total sample size
  active_arms <- colSums(apply(SS_matrix, 2, is.na)==0) # active arms per period
  block_sizes <- 2*active_arms # block sizes per period
  
  
  t <- c()
  
  for (i in 1:num_periods){
    m_i <- t(replicate(trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i]),
                       sample(rep(rep(c(0:(num_arms)), alloc_ratios[,i]), block_sizes[i]/length(rep(c(0:(num_arms)), alloc_ratios[,i]))))))
    
    t_i <- c(t(m_i), sample(rep(c(0:(num_arms)), alloc_ratios[,i]),
                            size = sum(SS_matrix[,i], na.rm = T)-block_sizes[i]*trunc(sum(SS_matrix[,i], na.rm = T)/block_sizes[i])))
    
    t <- c(t, t_i)
  }
  
  
  j <- c(1:n_total)
  
  
  for (i in 0:num_arms) {
    assign(paste0("j", i), which(t==i)) # j0, j1, j2 ... position in time (order) of allocated patients in every arm
  }
  
  cj <- rep(1:num_periods, N_period) # period indicator
  
  
  
  
  # Simulation of individual trend
  
  if(trend=="linear"){
    ind_trend0 <- linear_trend(j=j0,
                               lambda=lambda[1],sample_size=n_total)
    ind_trend1 <- linear_trend(j=j1,
                               lambda=lambda[2],sample_size=n_total)
    ind_trend2 <- linear_trend(j=j2,
                               lambda=lambda[3],sample_size=n_total)
  }
  
  if(trend=="linear2"){
    ind_trend0 <- linear_trend2(j=j0,
                               lambda=lambda[1],sample_size=c(sum(N1),sum(N2)+sum(N3)))
    ind_trend1 <- linear_trend2(j=j1,
                               lambda=lambda[2],sample_size=c(sum(N1),sum(N2)+sum(N3)))
    ind_trend2 <- linear_trend2(j=j2,
                               lambda=lambda[3],sample_size=c(sum(N1),sum(N2)+sum(N3)))
  }
  
  if(trend=="stepwise"){
    ind_trend0 <- sw_trend(cj=cj[j0], lambda=lambda[1])
    ind_trend1 <- sw_trend(cj=cj[j1], lambda=lambda[2])
    ind_trend2 <- sw_trend(cj=cj[j2], lambda=lambda[3])
  }
  
  if(trend=="inv_u"){
    
    j0_1 <- which(j0<=N_peak)
    j0_2 <- which(j0>N_peak)
    j1_1 <- which(j1<=N_peak)
    j1_2 <- which(j1>N_peak)
    j2_1 <- which(j2<=N_peak)
    j2_2 <- which(j2>N_peak)
    
    ind_trend0_1 <- linear_trend(j=j0[j0_1],
                                 lambda=lambda[1],sample_size=n_total)
    ind_trend0_2 <- linear_trend(j=j0[j0_2]-2*N_peak+1,
                                 lambda=-lambda[1],sample_size=n_total)
    ind_trend0 <- c(ind_trend0_1, ind_trend0_2)
    
    
    ind_trend1_1 <- linear_trend(j=j1[j1_1],
                                 lambda=lambda[2],sample_size=n_total)
    ind_trend1_2 <- linear_trend(j=j1[j1_2]-2*N_peak+1,
                                 lambda=-lambda[2],sample_size=n_total)
    ind_trend1 <- c(ind_trend1_1, ind_trend1_2)
    
    
    ind_trend2_1 <- linear_trend(j=j2[j2_1],
                                 lambda=lambda[3],sample_size=n_total)
    ind_trend2_2 <- linear_trend(j=j2[j2_2]-2*N_peak+1,
                                 lambda=-lambda[3],sample_size=n_total)
    ind_trend2 <- c(ind_trend2_1, ind_trend2_2)
  }
  
  # Simulation of continuous endpoint
  
  if(endpoint=="continuous"){
    
    means <- c()
    means[j0] <- ind_trend0
    means[j1] <- ind_trend1 + delta[1]
    means[j2] <- ind_trend2 + delta[2]
    
    X <- rnorm(n=n_total, mean=mu0+means, sd=sigma)
    
    Data <- data.frame(response = X,
                       treatment = t,
                       stage = c(rep(1, sum(N1)), rep(2, sum(N2)), rep(3, sum(N3))),
                       j = c(1:n_total),
                       lambda0 = lambda[1],
                       lambda1 = lambda[2],
                       lambda2 = lambda[3],
                       means = mu0+means)
  }
  
  # Simulation of binary endpoint (with different possible parametrizations of the time trend)
  
  if(endpoint=="binary"){
    
    # multiplicative effect on the odds ratio
    if(trend_param=="mult"){
      eta0 = log(p0/(1-p0)) + ind_trend0
      eta1 = log(p0/(1-p0)) + log(OR[1]) + ind_trend1
      eta2 = log(p0/(1-p0)) + log(OR[2]) + ind_trend2
      
      p <- c()
      p[j0] <- 1 / (1 + exp(-eta0))
      p[j1] <- 1 / (1 + exp(-eta1))
      p[j2] <- 1 / (1 + exp(-eta2))
  
      OR_inf2 <- (p[j2][1]/(1-p[j2][1]))/(p[j0][N1[1]+1]/(1-p[j0][N1[1]+1])) # time dependent OR2
    }
    
    # additive effect on the probabilities
    if(trend_param=="add"){
      O1 = OR[1]*p0/(1-p0)
      O2 = OR[2]*p0/(1-p0)
      
      p <- c()
      p[j0] = p0 + (ind_trend0)/10
      p[j1] = O1/(1+O1) + (ind_trend1)/10
      p[j2] = O2/(1+O2) + (ind_trend2)/10
      
      OR_inf2 <- (p[j2][1]/(1-p[j2][1]))/(p[j0][N1[1]+1]/(1-p[j0][N1[1]+1])) # time dependent OR2
    }
    
    if(sum(p<0 | p>1)>0){ # check if all probabilities are between 0 and 1
      stop("p must be between 0 and 1")
    }
    
    X <- rbinom(n = n_total, size = 1, prob = p)
    
    Data <- data.frame(response = X,
                       treatment = t,
                       stage = c(rep(1, sum(N1)), rep(2, sum(N2)), rep(3, sum(N3))),
                       j = c(1:n_total),
                       p = p,
                       OR2 = OR[2],
                       OR_inf2 = OR_inf2,
                       lambda0 = lambda[1],
                       lambda1 = lambda[2],
                       lambda2 = lambda[3],
                       trend_param = trend_param)
  }
  
  return(Data)
}

