#######################################################################################################################################

data_sim <- function(K=1, mu0=0, delta, p0, OR, lambda, sigma, N1, N2, N_add, trend, trend_param, endpoint){
  N1<-round(N1)
  N2<-round(N2)
  N_add<-round(N_add)
  
  N <- sum(N1+N2)
  
  t1 <- c(rep(0,N1[1]),rep(1,N1[2])) # vec of treatments for stage 1
  t1 <- t1[sample(1:sum(N1))] # randomize
  t2 <- c(rep(0,N2[1]),rep(1,N2[2]),rep(2,N_add)) # vec of treatments for stage 2
  t2 <- t2[sample(1:sum(N2, N_add))] # randomize
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
  
  # Simulation of continuous endpoint
  
  if(endpoint=="continuous"){
    
    means <- c()
    means[j0] <- ind_trend0
    means[j1] <- ind_trend1 + delta[1]
    means[j2] <- ind_trend2 + delta[2]
    
    X <- rnorm(n=sum(N+N_add), mean=mu0+means, sd=sigma)
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
    }
    
    # additive effect on the probabilities
    if(trend_param=="add"){
      O1 = OR[1]*p0/(1-p0)
      O2 = OR[2]*p0/(1-p0)
      
      p <- c()
      p[j0] = p0 + ind_trend0
      p[j1]  = O1/(1+O1) + ind_trend1
      p[j2] = O2/(1+O2) + ind_trend2
    }
    
    if(sum(p<0 | p>1)>0){ # check if all probabilities are between 0 and 1
      stop("p must be between 0 and 1")
    }

    X <- rbinom(n = sum(N+N_add), size = 1, prob = p)
  }
  
  Data <- tibble(response = X,
                 treatment = t,
                 stage = c(rep(1, sum(N1)), rep(2, sum(N2)+N_add)),
                 j = c(1:sum(N+N_add)))
  return(Data)
}


#######################################################################################################################################




