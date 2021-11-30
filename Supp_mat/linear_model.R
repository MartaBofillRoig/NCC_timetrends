
######################################################################################################################################

# Linear model M_a1 - using all data to estimate effects of all treatments + effect of continuous recruitment

linear_model_a1 <- function(data, delta, alpha=0.025){
  
  # fit linear model
  mod <- lm(response ~ as.factor(treatment) + j, data)
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pt(coef(res)["as.factor(treatment)2","t value"], mod$df, lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - delta
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat,bias=bias,reject_h02=reject_h02))
}

#linear_model_a1(data, delta = 0.15)

######################################################################################################################################

# Linear model M_a2 - using all data to estimate effects of all treatments + stage effect

linear_model_a2 <- function(data, delta, alpha=0.025){
  
  # fit linear model
  mod <- lm(response ~ as.factor(treatment) + as.factor(stage), data)
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pt(coef(res)["as.factor(treatment)2","t value"], mod$df, lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - delta
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat,bias=bias,reject_h02=reject_h02))
}

#linear_model_a2(data, delta = 0.15)

######################################################################################################################################

# Linear model M_b1 - using only control data (both stages) and newly added arm to estimate the effect of new treatment only + effect of continuous recruitment

linear_model_b1 <- function(data, delta, alpha=0.025){
  data <- data %>% filter(treatment %in% c(0, max(treatment)))
  
  # fit linear model
  mod <- lm(response ~ as.factor(treatment) + j, data)
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pt(coef(res)["as.factor(treatment)2","t value"], mod$df, lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - delta
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat,bias=bias,reject_h02=reject_h02))
}

######################################################################################################################################

# Linear model M_b2 - using only control data (both stages) and newly added arm to estimate the effect of new treatment only + stage effect

linear_model_b2 <- function(data, delta, alpha=0.025){
  data <- data %>% filter(treatment %in% c(0, max(treatment)))
  
  # fit linear model
  mod <- lm(response ~ as.factor(treatment) + as.factor(stage), data)
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pt(coef(res)["as.factor(treatment)2","t value"], mod$df, lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - delta
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat,bias=bias,reject_h02=reject_h02))
}

