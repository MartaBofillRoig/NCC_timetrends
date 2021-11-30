# Pooled z-test

z_test_pol <- function(data, mu0=0, sigma=1, delta, alpha=0.025){
  
  x0 <- data$response[data$treatment==0]
  x1 <- data$response[data$treatment==2]
  
  z_alpha <- qnorm(1-alpha,0,1)
  n_1 <- length(x1)
  n_0 <- length(x0)
  
  # calculate the z-statistic
  z_stat <- (mean(x1) - mean(x0) - mu0) / 
    sqrt(sigma^2/n_1 + sigma^2/n_0)
  
  # metrics
  bias <- (mean(x1) - mean(x0)) - delta
  reject_h02 <- (z_stat > z_alpha)
  pvalue <- 1-pnorm(z_stat)
  
  return(list(z_stat=z_stat, bias=bias, reject_h02=reject_h02, p_val=pvalue))
} 

# Separate z-test

z_test_sep <- function(data, mu0=0, sigma=1, delta, alpha=0.025){
  
  x0 <- data$response[data$treatment==0 & data$stage==2]
  x1 <- data$response[data$treatment==2]
  
  z_alpha <- qnorm(1-alpha,0,1)
  n_1 <- length(x1)
  n_0 <- length(x0)
  
  # calculate the z-statistic
  z_stat <- (mean(x1) - mean(x0) - mu0) / 
    sqrt(sigma^2/n_1 + sigma^2/n_0)
  
  # metrics
  bias <- (mean(x1) - mean(x0)) - delta
  reject_h02 <- (z_stat > z_alpha)
  pvalue <- 1-pnorm(z_stat)
  
  return(list(z_stat=z_stat, bias=bias, reject_h02=reject_h02, p_val=pvalue))
}
