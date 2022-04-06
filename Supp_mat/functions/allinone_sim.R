######################################################################################################################################

allinone_simsce <- function(nsim, scenarios, endpoint, alpha=0.025){
  
  if(endpoint=="continuous"){ 
    results <- data.frame(reject_h02_a1=c(0), reject_h02_a2=c(0),
                          reject_h02_a1_int=c(0), reject_h02_a2_int=c(0),
                          reject_h02_b1=c(0), reject_h02_b2=c(0),
                          reject_h02_tsep=c(0), reject_h02_tpol=c(0),
                          reject_h02_zsep=c(0), reject_h02_zpol=c(0),
                          rMSE_a1=c(0), rMSE_a2=c(0),
                          rMSE_a1_int=c(0), rMSE_a2_int=c(0),
                          rMSE_b1=c(0), rMSE_b2=c(0),
                          rMSE_tsep=c(0), rMSE_tpol=c(0),
                          rMSE_zsep=c(0), rMSE_zpol=c(0),
                          bias_a1=c(0), bias_a2=c(0),
                          bias_a1_int=c(0), bias_a2_int=c(0),
                          bias_b1=c(0), bias_b2=c(0),
                          bias_tsep=c(0), bias_tpol=c(0),
                          bias_zsep=c(0), bias_zpol=c(0))
    
    for(i in 1:dim(scenarios)[1]){
      res <- replicate(nsim, allinone_model(data=data_sim_block(K=scenarios$K[i], mu0=scenarios$mu0[i],
                                                                delta=c(scenarios$delta1[i], scenarios$delta2[i]),
                                                                lambda=c(scenarios$lambda0[i], scenarios$lambda1[i], scenarios$lambda2[i]),
                                                                sigma=scenarios$sigma[i],
                                                                N1=c(scenarios$n01[i], scenarios$n11[i]),
                                                                N2=c(scenarios$n02[i], scenarios$n12[i], scenarios$n22[i]),
                                                                N3=c(scenarios$n03[i], scenarios$n23[i]),
                                                                N_peak=scenarios$N_peak[i],
                                                                trend = scenarios$trend[i],
                                                                endpoint = scenarios$endpoint[i]),
                                            delta = scenarios$delta2[i],
                                            endpoint = scenarios$endpoint[i],
                                            alpha = alpha))
      
      db <- as.data.frame(t(res))
      prob_rej <- colMeans(do.call(rbind,db$reject_h0))
      rMSE <- sqrt(colMeans(do.call(rbind,db$bias)^2))
      bias <- colMeans(do.call(rbind,db$bias))
      results[i,] <- c(prob_rej, rMSE, bias)
    }
    return(cbind(scenarios,results))
  }
  
  if(endpoint=="binary"){
    results <- data.frame(reject_h02_a1=c(0), reject_h02_a2=c(0),
                          reject_h02_a1_int=c(0), reject_h02_a2_int=c(0),
                          reject_h02_b1=c(0), reject_h02_b2=c(0),
                          reject_h02_zsep=c(0), reject_h02_zpol=c(0),
                          reject_h02_log_sep=c(0), reject_h02_log_pol=c(0),
                          rMSE_a1=c(0), rMSE_a2=c(0),
                          rMSE_a1_int=c(0), rMSE_a2_int=c(0),
                          rMSE_b1=c(0), rMSE_b2=c(0),
                          rMSE_zsep=c(0), rMSE_zpol=c(0),
                          rMSE_log_sep=c(0), rMSE_log_pol=c(0),
                          bias_a1=c(0), bias_a2=c(0),
                          bias_a1_int=c(0), bias_a2_int=c(0),
                          bias_b1=c(0), bias_b2=c(0),
                          bias_zsep=c(0), bias_zpol=c(0),
                          bias_log_sep=c(0), bias_log_pol=c(0))
    
    for(i in 1:dim(scenarios)[1]){
      res <- replicate(nsim, allinone_model(data=db<-data_sim_block(K=scenarios$K[i], p0=scenarios$p0[i],
                                                                    OR=c(scenarios$OR1[i], scenarios$OR2[i]),
                                                                    lambda=c(scenarios$lambda0[i],
                                                                             scenarios$lambda1[i],
                                                                             scenarios$lambda2[i]),
                                                                    sigma=scenarios$sigma[i],
                                                                    N1=c(scenarios$n01[i], scenarios$n11[i]),
                                                                    N2=c(scenarios$n02[i], scenarios$n12[i], scenarios$n22[i]),
                                                                    N3=c(scenarios$n03[i], scenarios$n23[i]),
                                                                    N_peak=scenarios$N_peak[i],
                                                                    trend = scenarios$trend[i],
                                                                    endpoint = scenarios$endpoint[i],
                                                                    trend_param = scenarios$trend_param[i]),
                                            OR = db$OR_inf2[1],
                                            endpoint = scenarios$endpoint[i],
                                            alpha = alpha))
      
      db <- as.data.frame(t(res))
      prob_rej <- colMeans(do.call(rbind,db$reject_h0))
      rMSE <- sqrt(colMeans(do.call(rbind,db$bias)^2))
      bias <- colMeans(do.call(rbind,db$bias))
      results[i,] <- c(prob_rej, rMSE, bias)
    }
    return(cbind(scenarios,results))
  }
  
}




