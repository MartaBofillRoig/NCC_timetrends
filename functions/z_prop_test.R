#######################################################################################################################################

# Pooled two-proportions z-test

z_prop_pol <- function(data, OR, alpha=0.025){
  
  x1 <- data$response[data$treatment==0]
  x2 <- data$response[data$treatment==2]
  
  z_alpha <- qnorm(1-alpha,0,1)
  n1 <- length(x1)
  n2 <- length(x2)
  
  p1 <- sum(x1)/n1 # proportion 1
  p2 <- sum(x2)/n2 # proportion 2
  
  p <- sum(x1, x2)/(n1+n2) # overall proportion
  
  # calculate the z-statistic
  z_stat <- log((p2/(1-p2))/(p1/(1-p1)))/sqrt((1/n1 + 1/n2)*(1/(p*(1-p))))
  
  # metrics
  bias <- log((p2/(1-p2))/(p1/(1-p1))) - log(OR)
  reject_h02 <- (z_stat > z_alpha)
  pvalue <- 1-pnorm(z_stat)
  
  return(list(z_stat=z_stat, bias=bias, reject_h02=reject_h02, p_val=pvalue))
}

#######################################################################################################################################

# Separate two-proportions z-test

z_prop_sep <- function(data, OR, alpha=0.025){
  
  x1 <- data$response[data$treatment==0 & data$stage==2]
  x2 <- data$response[data$treatment==2]
  
  z_alpha <- qnorm(1-alpha,0,1)
  n1 <- length(x1)
  n2 <- length(x2)
  
  p1 <- sum(x1)/n1 # proportion 1
  p2 <- sum(x2)/n2 # proportion 2
  
  p <- sum(x1, x2)/(n1+n2) # overall proportion
  
  # calculate the z-statistic
  z_stat <- log((p2/(1-p2))/(p1/(1-p1)))/sqrt((1/n1 + 1/n2)*(1/(p*(1-p))))
  
  # metrics
  bias <- log((p2/(1-p2))/(p1/(1-p1))) - log(OR)
  reject_h02 <- (z_stat > z_alpha)
  pvalue <- 1-pnorm(z_stat)
  
  return(list(z_stat=z_stat, bias=bias, reject_h02=reject_h02, p_val=pvalue))
}

#######################################################################################################################################
