#######################################################################################################################################

# Pooled t-test

t_test_pol <- function(data, delta, alpha=0.025){
  
  x0 <- data$response[data$treatment==0]
  x1 <- data$response[data$treatment==2]
  
  n_1 <- length(x1)
  n_0 <- length(x0)
  
  # calculate the t-test 
  t_stat <- t.test(x1, x0, alternative = "greater")
    
  # metrics
  bias <- (t_stat$estimate[1]-t_stat$estimate[2]) - delta
  reject_h0 <- (t_stat$p.value < alpha)
  
  return(list(t_test=t_stat, bias=bias, reject_h02=reject_h0))
}

#######################################################################################################################################

# Separate t-test

t_test_sep <- function(data, delta, alpha=0.025){
  
  x0 <- data$response[data$treatment==0 & data$stage==2]
  x1 <- data$response[data$treatment==2]
  
  n_1 <- length(x1)
  n_0 <- length(x0)
  
  # calculate the z-statistic
  t_stat <- t.test(x1, x0, alternative = "greater", var.equal = T)
  
  # metrics
  bias <- (t_stat$estimate[1]-t_stat$estimate[2]) - delta
  reject_h0 <- (t_stat$p.value < alpha)
  
  return(list(t_test=t_stat, bias=bias, reject_h02=reject_h0)) 
} 

#######################################################################################################################################