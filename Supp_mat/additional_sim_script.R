rm(list = ls())

source("functions/linear_model.R")
source("functions/log_model.R")
source("functions/z_test.R")
source("functions/t_test.R")
source("functions/z_prop_test.R")
source("functions/trend_functions.R")
source("functions/data_sim_block.R")
source("functions/allinone_model.R")
source("functions/allinone_sim.R")
source("functions/allinone_sim_par.R")
library(tidyverse)
library(foreach)
library(doParallel)
library(rpact)

############################################################################################################
# Treatment effects

#power.t.test(n=250, sd=1, sig.level=0.025, power=0.8, type="two.sample", alternative="one.sided") # delta2=0.25

#power.prop.test(n=250, p1=0.7, sig.level=0.025, power=0.8, alternative="one.sided") # p2=0.8076999

#(0.8076999/(1-0.8076999)) / (0.7/(1-0.7)) # OR2=1.8

#power.prop.test(n=230, p1=0.7, sig.level=0.025, power=0.8, alternative="one.sided") # p2=0.8119242

#(0.8119242/(1-0.8119242)) / (0.7/(1-0.7)) # OR2=1.85

#z.alpha <- qnorm(1-0.025,0,1)
#z.beta <-  qnorm(1-0.2,0,1)

#((z.alpha+z.beta)^2/(log(1.8)^2)) * ((1/(0.7*(1-0.7)))+(1/(0.8076999*(1-0.8076999))))


# Sample sizes for new trial scheme - continuous case

#getSampleSizeMeans(alternative = 0.25, allocationRatioPlanned = 1/1.5) # n0=315; n1=n2=210

# Sample sizes for new trial scheme - binary case

#getSampleSizeRates(pi1 = 0.8076999, pi2 = 0.7, allocationRatioPlanned = 1/1.5) # n0=316; n1=n2=211

############################################################################################################
# Number of replications
nsim=100000


############################################################################################################
# BINARY CASE - linear trend - OR1>1 - equal time trends - trial scheme 2 (N_1 = 2 x n/2)
############################################################################################################


# Binary scenario (equal time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_lin_mult_alpha_OR1_eq_2 <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 1.8,
                                                    OR2 = 1,
                                                    lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                    n01 = 210*(1/2),
                                                    n11 = 210*(1/2),
                                                    n02 = 210*(1/2),
                                                    n12 = 210*(1/2),
                                                    n22 = 210*(1/2),
                                                    n03 = 210*(1/2),
                                                    n23 = 210*(1/2),
                                                    trend = "linear",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = (n01+n02+n03)/(n22+n23))



set.seed(1)
results_bin_lin_mult_alpha_OR1_eq_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_OR1_eq_2, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_OR1_eq_2, "results/results_bin_lin_mult_alpha_OR1_eq_2.csv")
gc()


############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - power 

scenarios_bin_lin_mult_pow_OR1_eq_2 <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 1.8,
                                                  OR2 = 1.8,
                                                  lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                  n01 = 210*(1/2),
                                                  n11 = 210*(1/2),
                                                  n02 = 210*(1/2),
                                                  n12 = 210*(1/2),
                                                  n22 = 210*(1/2),
                                                  n03 = 210*(1/2),
                                                  n23 = 210*(1/2),
                                                  trend = "linear",
                                                  trend_param = "mult",
                                                  endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = (n01+n02+n03)/(n22+n23))


set.seed(2)
results_bin_lin_mult_pow_OR1_eq_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_OR1_eq_2, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_OR1_eq_2, "results/results_bin_lin_mult_pow_OR1_eq_2.csv")
gc()

############################################################################################################







############################################################################################################
# BINARY CASE - linear trend - OR1>1 - different time trends, negative trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_lin_mult_alpha_OR1_diff_neg_2 <- data.frame(K = 1,
                                                          p0 = 0.7,
                                                          sigma = 1,
                                                          OR1 = 1.8,
                                                          OR2 = 1,
                                                          lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                          lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                          lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                          n01 = 210*(1/2),
                                                          n11 = 210*(1/2),
                                                          n02 = 210*(1/2),
                                                          n12 = 210*(1/2),
                                                          n22 = 210*(1/2),
                                                          n03 = 210*(1/2),
                                                          n23 = 210*(1/2),
                                                          trend = "linear",
                                                          trend_param = "mult",
                                                          endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = (n01+n02+n03)/(n22+n23))



set.seed(3)
results_bin_lin_mult_alpha_OR1_diff_neg_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_OR1_diff_neg_2, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_OR1_diff_neg_2, "results/results_bin_lin_mult_alpha_OR1_diff_neg_2.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_lin_mult_pow_OR1_diff_neg_2 <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1.8,
                                                        lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        n01 = 210*(1/2),
                                                        n11 = 210*(1/2),
                                                        n02 = 210*(1/2),
                                                        n12 = 210*(1/2),
                                                        n22 = 210*(1/2),
                                                        n03 = 210*(1/2),
                                                        n23 = 210*(1/2),
                                                        trend = "linear",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = (n01+n02+n03)/(n22+n23))


set.seed(4)
results_bin_lin_mult_pow_OR1_diff_neg_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_OR1_diff_neg_2, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_OR1_diff_neg_2, "results/results_bin_lin_mult_pow_OR1_diff_neg_2.csv")
gc()

############################################################################################################






############################################################################################################
# BINARY CASE - linear trend - OR1>1 - different time trends, positive trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_lin_mult_alpha_OR1_diff_pos_2 <- data.frame(K = 1,
                                                          p0 = 0.7,
                                                          sigma = 1,
                                                          OR1 = 1.8,
                                                          OR2 = 1,
                                                          lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                          lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                          lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                          n01 = 210*(1/2),
                                                          n11 = 210*(1/2),
                                                          n02 = 210*(1/2),
                                                          n12 = 210*(1/2),
                                                          n22 = 210*(1/2),
                                                          n03 = 210*(1/2),
                                                          n23 = 210*(1/2),
                                                          trend = "linear",
                                                          trend_param = "mult",
                                                          endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = (n01+n02+n03)/(n22+n23))



set.seed(5)
results_bin_lin_mult_alpha_OR1_diff_pos_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_OR1_diff_pos_2, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_OR1_diff_pos_2, "results/results_bin_lin_mult_alpha_OR1_diff_pos_2.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_lin_mult_pow_OR1_diff_pos_2 <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1.8,
                                                        lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        n01 = 210*(1/2),
                                                        n11 = 210*(1/2),
                                                        n02 = 210*(1/2),
                                                        n12 = 210*(1/2),
                                                        n22 = 210*(1/2),
                                                        n03 = 210*(1/2),
                                                        n23 = 210*(1/2),
                                                        trend = "linear",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = (n01+n02+n03)/(n22+n23))


set.seed(6)
results_bin_lin_mult_pow_OR1_diff_pos_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_OR1_diff_pos_2, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_OR1_diff_pos_2, "results/results_bin_lin_mult_pow_OR1_diff_pos_2.csv")
gc()

############################################################################################################









############################################################################################################
# CONTINUOUS CASE - equal time trends - trial scheme 2 (N_1 = 2 x n/2)
############################################################################################################

# Continuous case (equal time trends) - linear trend - type one error


scenarios_cont_lin_alpha_eq_2 <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0,
                                            lambda0 = seq(-0.15, 0.15, length.out = 19),
                                            lambda1 = seq(-0.15, 0.15, length.out = 19),
                                            lambda2 = seq(-0.15, 0.15, length.out = 19),
                                            n01 = 210*(1/2),
                                            n11 = 210*(1/2),
                                            n02 = 210*(1/2),
                                            n12 = 210*(1/2),
                                            n22 = 210*(1/2),
                                            n03 = 210*(1/2),
                                            n23 = 210*(1/2),
                                            trend = "linear",
                                            endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = (n01+n02+n03)/(n22+n23))


set.seed(1)
results_cont_lin_alpha_eq_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_alpha_eq_2, endpoint = "continuous")

write_csv(results_cont_lin_alpha_eq_2, "results/results_cont_lin_alpha_eq_2.csv")
gc()


############################################################################################################

# Continuous case (equal time trends) - linear trend - power

scenarios_cont_lin_pow_eq_2 <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0.25,
                                          lambda0 = seq(-0.15, 0.15, length.out = 19),
                                          lambda1 = seq(-0.15, 0.15, length.out = 19),
                                          lambda2 = seq(-0.15, 0.15, length.out = 19),
                                          n01 = 210*(1/2),
                                          n11 = 210*(1/2),
                                          n02 = 210*(1/2),
                                          n12 = 210*(1/2),
                                          n22 = 210*(1/2),
                                          n03 = 210*(1/2),
                                          n23 = 210*(1/2),
                                          trend = "linear",
                                          endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = (n01+n02+n03)/(n22+n23))



set.seed(2)
results_cont_lin_pow_eq_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_pow_eq_2, endpoint = "continuous")

write_csv(results_cont_lin_pow_eq_2, "results/results_cont_lin_pow_eq_2.csv")
gc()


############################################################################################################





############################################################################################################
# CONTINUOUS CASE - different time trends
############################################################################################################

# Continuous case (different time trends) - linear trend - type one error


scenarios_cont_lin_alpha_diff_2 <- data.frame(K = 1,
                                              mu0 = 0,
                                              sigma = 1,
                                              delta1 = 0.25,
                                              delta2 = 0,
                                              lambda0 = 0.1,
                                              lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                              lambda2 = 0.1,
                                              n01 = 210*(1/2),
                                              n11 = 210*(1/2),
                                              n02 = 210*(1/2),
                                              n12 = 210*(1/2),
                                              n22 = 210*(1/2),
                                              n03 = 210*(1/2),
                                              n23 = 210*(1/2),
                                              trend = "linear",
                                              endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = (n01+n02+n03)/(n22+n23))


set.seed(3)
results_cont_lin_alpha_diff_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_alpha_diff_2, endpoint = "continuous")

write_csv(results_cont_lin_alpha_diff_2, "results/results_cont_lin_alpha_diff_2.csv")
gc()



############################################################################################################

# Continuous case (different time trends) - linear trend - power

scenarios_cont_lin_pow_diff_2 <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0.25,
                                            lambda0 = 0.1,
                                            lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                            lambda2 = 0.1,
                                            n01 = 210*(1/2),
                                            n11 = 210*(1/2),
                                            n02 = 210*(1/2),
                                            n12 = 210*(1/2),
                                            n22 = 210*(1/2),
                                            n03 = 210*(1/2),
                                            n23 = 210*(1/2),
                                            trend = "linear",
                                            endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = (n01+n02+n03)/(n22+n23))



set.seed(4)
results_cont_lin_pow_diff_2 <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_pow_diff_2, endpoint = "continuous")

write_csv(results_cont_lin_pow_diff_2, "results/results_cont_lin_pow_diff_2.csv")
gc()


############################################################################################################



