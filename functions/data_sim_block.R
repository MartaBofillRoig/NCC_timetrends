#######################################################################################################################################

data_sim_block <- function(K=1, mu0=0, delta, p0, OR, lambda, sigma, N1, N2, N_add, N_peak, trend, trend_param, endpoint){
  N1<-round(N1)
  N2<-round(N2)
  N_add<-round(N_add)
  
  N <- sum(N1+N2)
  
  ar1 <- N1/(min(N1))
  ar2 <- c(N2,N_add)/min(c(N2,N_add))
  
  if(sum(ar1==1)==2 & sum(ar2==c(1,1,2))==3){
    
    m <- t(replicate(trunc(sum(N1)/4),sample(rep(c(0,1),2))))
    t1 <- c(t(m),sample(c(0,0,1,1),size=sum(N1)-4*round(sum(N1)/4)))  
    
    m <- t(replicate(trunc(sum(N2,N_add)/12),sample( rep(c(0,1,2,2),3) )))
    t2 <- c(t(m),sample(rep(c(0,1,2,2),3),size=sum(N2, N_add)-12*trunc(sum(N2,N_add)/12)))  
    
  }else{stop("1:1 period 1, and 1:1:2 period 2")}
  
  t <- c(t1,t2)
  
  j0 <- which(t==0)
  j1 <- which(t==1)
  j2 <- which(t==2)
  
  cj <- c(rep(FALSE, sum(N1)), rep(TRUE, sum(N2)+N_add)) # indicator whether subject j is in stage 2
  
  
  # Simulation of individual trend
  
  if(trend=="linear"){
    ind_trend0 <- linear_trend(j=j0,
                               lambda=lambda[1],sample_size=sum(N+N_add))
    ind_trend1 <- linear_trend(j=j1,
                               lambda=lambda[2],sample_size=sum(N+N_add))
    ind_trend2 <- linear_trend(j=j2,
                               lambda=lambda[3],sample_size=sum(N+N_add))
  }
  
  if(trend=="linear2"){
    ind_trend0 <- linear_trend2(j=j0,
                               lambda=lambda[1],sample_size=c(sum(N1),sum(N2)+N_add))
    ind_trend1 <- linear_trend2(j=j1,
                               lambda=lambda[2],sample_size=c(sum(N1),sum(N2)+N_add))
    ind_trend2 <- linear_trend2(j=j2,
                               lambda=lambda[3],sample_size=c(sum(N1),sum(N2)+N_add))
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
                                 lambda=lambda[1],sample_size=sum(N+N_add))
    ind_trend0_2 <- linear_trend(j=j0[j0_2]-2*N_peak+1,
                                 lambda=-lambda[1],sample_size=sum(N+N_add))
    ind_trend0 <- c(ind_trend0_1, ind_trend0_2)
    
    
    ind_trend1_1 <- linear_trend(j=j1[j1_1],
                                 lambda=lambda[2],sample_size=sum(N+N_add))
    ind_trend1_2 <- linear_trend(j=j1[j1_2]-2*N_peak+1,
                                 lambda=-lambda[2],sample_size=sum(N+N_add))
    ind_trend1 <- c(ind_trend1_1, ind_trend1_2)
    
    
    ind_trend2_1 <- linear_trend(j=j2[j2_1],
                                 lambda=lambda[3],sample_size=sum(N+N_add))
    ind_trend2_2 <- linear_trend(j=j2[j2_2]-2*N_peak+1,
                                 lambda=-lambda[3],sample_size=sum(N+N_add))
    ind_trend2 <- c(ind_trend2_1, ind_trend2_2)
  }
  
  # Simulation of continuous endpoint
  
  if(endpoint=="continuous"){
    
    means <- c()
    means[j0] <- ind_trend0
    means[j1] <- ind_trend1 + delta[1]
    means[j2] <- ind_trend2 + delta[2]
    
    X <- rnorm(n=sum(N+N_add), mean=mu0+means, sd=sigma)
    
    Data <- data.frame(response = X,
                       treatment = t,
                       stage = c(rep(1, sum(N1)), rep(2, sum(N2)+N_add)),
                       j = c(1:sum(N+N_add)),
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
    
    X <- rbinom(n = sum(N+N_add), size = 1, prob = p)
    
    Data <- data.frame(response = X,
                       treatment = t,
                       stage = c(rep(1, sum(N1)), rep(2, sum(N2)+N_add)),
                       j = c(1:sum(N+N_add)),
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

######################################################################################################################################


