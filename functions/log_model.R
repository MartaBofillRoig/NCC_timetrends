#######################################################################################################################################

# Logistic model M_a1 (ALLTC-Linear) - using all data to estimate effects of all treatments + effect of continuous recruitment

logistic_model_a1 <- function(data, OR, alpha=0.025){
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment) + j, data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model M_a2 (ALLTC-Step) - using all data to estimate effects of all treatments + stage effect

logistic_model_a2 <- function(data, OR, alpha=0.025){
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment) + as.factor(stage), data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model M_b1 (TC-Linear) - using only control data (both stages) and newly added arm to estimate the effect of new treatment only + effect of continuous recruitment

logistic_model_b1 <- function(data, OR, alpha=0.025){
  
  data <- data[data$treatment %in% c(0, 2),]
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment) + j, data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model M_b2 (TC-Step) - using only control data (both stages) and newly added arm to estimate the effect of new treatment only + stage effect

logistic_model_b2 <- function(data, OR, alpha=0.025){
  
  data <- data[data$treatment %in% c(0, 2),]
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment) + as.factor(stage), data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model M_a1_int (ALLTCI-Linear) - using all data to estimate effects of all treatments + effect of continuous recruitment (with interaction)

logistic_model_a1_int <- function(data, OR, alpha=0.025){
  
  data$trt_1 <- ifelse(data$treatment==1, 1, 0)
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment) + j + trt_1:j, data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model M_a2_int (ALLTCI-Step) - using all data to estimate effects of all treatments + stage effect (with interaction)

logistic_model_a2_int <- function(data, OR, alpha=0.025){
  
  data$trt_1 <- ifelse(data$treatment==1, 1, 0)
  data$stage_2 <- ifelse(data$stage==2, 1, 0)
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment) + stage_2 + trt_1:stage_2, data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model sepatate

logistic_model_sep <- function(data, OR, alpha=0.025){
  
  data <- data[data$treatment %in% c(0, 2) & data$stage==2,]
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment), data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

######################################################################################################################################

# Logistic model pooled

logistic_model_pol <- function(data, OR, alpha=0.025){
  
  data <- data[data$treatment %in% c(0, 2),]
  
  # fit logistic model
  mod <- glm(response ~ as.factor(treatment), data, family = "binomial")
  res <- summary(mod)
  
  # one-sided p-value for new treatment
  p_val_treat <- pnorm(coef(res)["as.factor(treatment)2", "z value"], lower = FALSE)
  
  # metrics
  bias <- res$coefficients["as.factor(treatment)2", "Estimate"] - log(OR)
  reject_h02 <- (p_val_treat < alpha)
  
  return(list(p_val=p_val_treat, bias=bias, reject_h02=reject_h02))
}

#######################################################################################################################################