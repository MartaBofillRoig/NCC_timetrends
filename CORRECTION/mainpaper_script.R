rm(list = ls())

source("linear_model.R")
source("log_model.R")
source("z_test.R")
source("t_test.R")
source("z_prop_test.R")
source("trend_functions.R")
source("data_sim_block.R")
source("allinone_model.R")
source("allinone_sim.R")
source("allinone_sim_par.R")
library(tidyverse)
library(foreach)
library(doParallel)

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



############################################################################################################
# Number of replications
nsim=100000

############################################################################################################
# BINARY CASE - stepwise trend - OR1<1 - equal time trends
############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - type one error

scenarios_bin_step_add_alpha_eq <- data.frame(K = 1,
                                              p0 = 0.7,
                                              sigma = 1,
                                              OR1 = 0.4,
                                              OR2 = 1,
                                              lambda0 = seq(-0.5, 0.5, length.out = 19),
                                              lambda1 = seq(-0.5, 0.5, length.out = 19),
                                              lambda2 = seq(-0.5, 0.5, length.out = 19),
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "stepwise",
                                              trend_param = "add",
                                              endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(1)
results_bin_step_add_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_add_alpha_eq, endpoint = "binary")

write_csv(results_bin_step_add_alpha_eq, "results/results_bin_step_add_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_step_mult_alpha_eq <- data.frame(K = 1,
                                               p0 = 0.7,
                                               sigma = 1,
                                               OR1 = 0.4,
                                               OR2 = 1,
                                               lambda0 = seq(-0.5, 0.5, length.out = 19),
                                               lambda1 = seq(-0.5, 0.5, length.out = 19),
                                               lambda2 = seq(-0.5, 0.5, length.out = 19),
                                               n0 = 250,
                                               n1 = 250,
                                               n01 = 250*0.5,
                                               n02 = 250-250*0.5,
                                               n11 = 250*0.5,
                                               n12 = 250-250*0.5,
                                               n22 = 250,
                                               trend = "stepwise",
                                               trend_param = "mult",
                                               endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(2)
results_bin_step_mult_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_alpha_eq, endpoint = "binary")

write_csv(results_bin_step_mult_alpha_eq, "results/results_bin_step_mult_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - power

scenarios_bin_step_add_pow_eq <- data.frame(K = 1,
                                            p0 = 0.7,
                                            sigma = 1,
                                            OR1 = 0.4,
                                            OR2 = 1.8,
                                            lambda0 = seq(-0.5, 0.5, length.out = 19),
                                            lambda1 = seq(-0.5, 0.5, length.out = 19),
                                            lambda2 = seq(-0.5, 0.5, length.out = 19),
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "stepwise",
                                            trend_param = "add",
                                            endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(3)
results_bin_step_add_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_add_pow_eq, endpoint = "binary")

write_csv(results_bin_step_add_pow_eq, "results/results_bin_step_add_pow_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - power

scenarios_bin_step_mult_pow_eq <- data.frame(K = 1,
                                             p0 = 0.7,
                                             sigma = 1,
                                             OR1 = 0.4,
                                             OR2 = 1.8,
                                             lambda0 = seq(-0.5, 0.5, length.out = 19),
                                             lambda1 = seq(-0.5, 0.5, length.out = 19),
                                             lambda2 = seq(-0.5, 0.5, length.out = 19),
                                             n0 = 250,
                                             n1 = 250,
                                             n01 = 250*0.5,
                                             n02 = 250-250*0.5,
                                             n11 = 250*0.5,
                                             n12 = 250-250*0.5,
                                             n22 = 250,
                                             trend = "stepwise",
                                             trend_param = "mult",
                                             endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(4)
results_bin_step_mult_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_pow_eq, endpoint = "binary")

write_csv(results_bin_step_mult_pow_eq, "results/results_bin_step_mult_pow_eq.csv")
gc()




############################################################################################################
# BINARY CASE - stepwise trend - OR1>1 - equal time trends
############################################################################################################


# Binary scenario (equal time trends) (OR>1) - additive case - type one error

scenarios_bin_step_add_alpha_OR1_eq <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 1.8,
                                                  OR2 = 1,
                                                  lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "stepwise",
                                                  trend_param = "add",
                                                  endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(5)
results_bin_step_add_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_add_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_step_add_alpha_OR1_eq, "results/results_bin_step_add_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_step_mult_alpha_OR1_eq <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 1.8,
                                                   OR2 = 1,
                                                   lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "stepwise",
                                                   trend_param = "mult",
                                                   endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(6)
results_bin_step_mult_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_step_mult_alpha_OR1_eq, "results/results_bin_step_mult_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - additive case - power 

scenarios_bin_step_add_pow_OR1_eq <- data.frame(K = 1,
                                                p0 = 0.7,
                                                sigma = 1,
                                                OR1 = 1.8,
                                                OR2 = 1.8,
                                                lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                n0 = 250,
                                                n1 = 250,
                                                n01 = 250*0.5,
                                                n02 = 250-250*0.5,
                                                n11 = 250*0.5,
                                                n12 = 250-250*0.5,
                                                n22 = 250,
                                                trend = "stepwise",
                                                trend_param = "add",
                                                endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(7)
results_bin_step_add_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_add_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_step_add_pow_OR1_eq, "results/results_bin_step_add_pow_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - power 

scenarios_bin_step_mult_pow_OR1_eq <- data.frame(K = 1,
                                                 p0 = 0.7,
                                                 sigma = 1,
                                                 OR1 = 1.8,
                                                 OR2 = 1.8,
                                                 lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                 n0 = 250,
                                                 n1 = 250,
                                                 n01 = 250*0.5,
                                                 n02 = 250-250*0.5,
                                                 n11 = 250*0.5,
                                                 n12 = 250-250*0.5,
                                                 n22 = 250,
                                                 trend = "stepwise",
                                                 trend_param = "mult",
                                                 endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(8)
results_bin_step_mult_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_step_mult_pow_OR1_eq, "results/results_bin_step_mult_pow_OR1_eq.csv")
gc()

############################################################################################################






























############################################################################################################
# BINARY CASE - stepwise trend - OR1<1 - different time trends, negative trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_step_mult_alpha_diff_neg <- data.frame(K = 1,
                                                     p0 = 0.7,
                                                     sigma = 1,
                                                     OR1 = 0.4,
                                                     OR2 = 1,
                                                     lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                     lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                     lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                     n0 = 250,
                                                     n1 = 250,
                                                     n01 = 250*0.5,
                                                     n02 = 250-250*0.5,
                                                     n11 = 250*0.5,
                                                     n12 = 250-250*0.5,
                                                     n22 = 250,
                                                     trend = "stepwise",
                                                     trend_param = "mult",
                                                     endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(9)
results_bin_step_mult_alpha_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_alpha_diff_neg, endpoint = "binary")

write_csv(results_bin_step_mult_alpha_diff_neg, "results/results_bin_step_mult_alpha_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_step_mult_pow_diff_neg <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 0.4,
                                                   OR2 = 1.8,
                                                   lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "stepwise",
                                                   trend_param = "mult",
                                                   endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(10)
results_bin_step_mult_pow_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_pow_diff_neg, endpoint = "binary")

write_csv(results_bin_step_mult_pow_diff_neg, "results/results_bin_step_mult_pow_diff_neg.csv")
gc()




############################################################################################################
# BINARY CASE - stepwise trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_step_mult_alpha_OR1_diff_neg <- data.frame(K = 1,
                                                         p0 = 0.7,
                                                         sigma = 1,
                                                         OR1 = 1.8,
                                                         OR2 = 1,
                                                         lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                         lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                         lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                         n0 = 250,
                                                         n1 = 250,
                                                         n01 = 250*0.5,
                                                         n02 = 250-250*0.5,
                                                         n11 = 250*0.5,
                                                         n12 = 250-250*0.5,
                                                         n22 = 250,
                                                         trend = "stepwise",
                                                         trend_param = "mult",
                                                         endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(11)
results_bin_step_mult_alpha_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_alpha_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_step_mult_alpha_OR1_diff_neg, "results/results_bin_step_mult_alpha_OR1_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_step_mult_pow_OR1_diff_neg <- data.frame(K = 1,
                                                       p0 = 0.7,
                                                       sigma = 1,
                                                       OR1 = 1.8,
                                                       OR2 = 1.8,
                                                       lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                       lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                       lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                       n0 = 250,
                                                       n1 = 250,
                                                       n01 = 250*0.5,
                                                       n02 = 250-250*0.5,
                                                       n11 = 250*0.5,
                                                       n12 = 250-250*0.5,
                                                       n22 = 250,
                                                       trend = "stepwise",
                                                       trend_param = "mult",
                                                       endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(12)
results_bin_step_mult_pow_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_pow_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_step_mult_pow_OR1_diff_neg, "results/results_bin_step_mult_pow_OR1_diff_neg.csv")
gc()

############################################################################################################














############################################################################################################
# BINARY CASE - stepwise trend - OR1<1 - different time trends, positive trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_step_mult_alpha_diff_pos <- data.frame(K = 1,
                                                     p0 = 0.7,
                                                     sigma = 1,
                                                     OR1 = 0.4,
                                                     OR2 = 1,
                                                     lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                     lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                     lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                     n0 = 250,
                                                     n1 = 250,
                                                     n01 = 250*0.5,
                                                     n02 = 250-250*0.5,
                                                     n11 = 250*0.5,
                                                     n12 = 250-250*0.5,
                                                     n22 = 250,
                                                     trend = "stepwise",
                                                     trend_param = "mult",
                                                     endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(13)
results_bin_step_mult_alpha_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_alpha_diff_pos, endpoint = "binary")

write_csv(results_bin_step_mult_alpha_diff_pos, "results/results_bin_step_mult_alpha_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_step_mult_pow_diff_pos <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 0.4,
                                                   OR2 = 1.8,
                                                   lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "stepwise",
                                                   trend_param = "mult",
                                                   endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(14)
results_bin_step_mult_pow_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_pow_diff_pos, endpoint = "binary")

write_csv(results_bin_step_mult_pow_diff_pos, "results/results_bin_step_mult_pow_diff_pos.csv")
gc()




############################################################################################################
# BINARY CASE - stepwise trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_step_mult_alpha_OR1_diff_pos <- data.frame(K = 1,
                                                         p0 = 0.7,
                                                         sigma = 1,
                                                         OR1 = 1.8,
                                                         OR2 = 1,
                                                         lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                         lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                         lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                         n0 = 250,
                                                         n1 = 250,
                                                         n01 = 250*0.5,
                                                         n02 = 250-250*0.5,
                                                         n11 = 250*0.5,
                                                         n12 = 250-250*0.5,
                                                         n22 = 250,
                                                         trend = "stepwise",
                                                         trend_param = "mult",
                                                         endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(15)
results_bin_step_mult_alpha_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_alpha_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_step_mult_alpha_OR1_diff_pos, "results/results_bin_step_mult_alpha_OR1_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_step_mult_pow_OR1_diff_pos <- data.frame(K = 1,
                                                       p0 = 0.7,
                                                       sigma = 1,
                                                       OR1 = 1.8,
                                                       OR2 = 1.8,
                                                       lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                       lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                       lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                       n0 = 250,
                                                       n1 = 250,
                                                       n01 = 250*0.5,
                                                       n02 = 250-250*0.5,
                                                       n11 = 250*0.5,
                                                       n12 = 250-250*0.5,
                                                       n22 = 250,
                                                       trend = "stepwise",
                                                       trend_param = "mult",
                                                       endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(16)
results_bin_step_mult_pow_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_step_mult_pow_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_step_mult_pow_OR1_diff_pos, "results/results_bin_step_mult_pow_OR1_diff_pos.csv")
gc()

############################################################################################################









































############################################################################################################
# BINARY CASE - linear trend - OR1<1 - equal time trends
############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - type one error

scenarios_bin_lin_add_alpha_eq <- data.frame(K = 1,
                                              p0 = 0.7,
                                              sigma = 1,
                                              OR1 = 0.4,
                                              OR2 = 1,
                                              lambda0 = seq(-0.5, 0.5, length.out = 19),
                                              lambda1 = seq(-0.5, 0.5, length.out = 19),
                                              lambda2 = seq(-0.5, 0.5, length.out = 19),
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "linear",
                                              trend_param = "add",
                                              endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(1)
results_bin_lin_add_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_add_alpha_eq, endpoint = "binary")

write_csv(results_bin_lin_add_alpha_eq, "results/results_bin_lin_add_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_lin_mult_alpha_eq <- data.frame(K = 1,
                                               p0 = 0.7,
                                               sigma = 1,
                                               OR1 = 0.4,
                                               OR2 = 1,
                                               lambda0 = seq(-0.5, 0.5, length.out = 19),
                                               lambda1 = seq(-0.5, 0.5, length.out = 19),
                                               lambda2 = seq(-0.5, 0.5, length.out = 19),
                                               n0 = 250,
                                               n1 = 250,
                                               n01 = 250*0.5,
                                               n02 = 250-250*0.5,
                                               n11 = 250*0.5,
                                               n12 = 250-250*0.5,
                                               n22 = 250,
                                               trend = "linear",
                                               trend_param = "mult",
                                               endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(2)
results_bin_lin_mult_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_eq, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_eq, "results/results_bin_lin_mult_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - power

scenarios_bin_lin_add_pow_eq <- data.frame(K = 1,
                                            p0 = 0.7,
                                            sigma = 1,
                                            OR1 = 0.4,
                                            OR2 = 1.8,
                                            lambda0 = seq(-0.5, 0.5, length.out = 19),
                                            lambda1 = seq(-0.5, 0.5, length.out = 19),
                                            lambda2 = seq(-0.5, 0.5, length.out = 19),
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "linear",
                                            trend_param = "add",
                                            endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(3)
results_bin_lin_add_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_add_pow_eq, endpoint = "binary")

write_csv(results_bin_lin_add_pow_eq, "results/results_bin_lin_add_pow_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - power

scenarios_bin_lin_mult_pow_eq <- data.frame(K = 1,
                                             p0 = 0.7,
                                             sigma = 1,
                                             OR1 = 0.4,
                                             OR2 = 1.8,
                                             lambda0 = seq(-0.5, 0.5, length.out = 19),
                                             lambda1 = seq(-0.5, 0.5, length.out = 19),
                                             lambda2 = seq(-0.5, 0.5, length.out = 19),
                                             n0 = 250,
                                             n1 = 250,
                                             n01 = 250*0.5,
                                             n02 = 250-250*0.5,
                                             n11 = 250*0.5,
                                             n12 = 250-250*0.5,
                                             n22 = 250,
                                             trend = "linear",
                                             trend_param = "mult",
                                             endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(4)
results_bin_lin_mult_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_eq, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_eq, "results/results_bin_lin_mult_pow_eq.csv")
gc()




############################################################################################################
# BINARY CASE - linear trend - OR1>1 - equal time trends
############################################################################################################


# Binary scenario (equal time trends) (OR>1) - additive case - type one error

scenarios_bin_lin_add_alpha_OR1_eq <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 1.8,
                                                  OR2 = 1,
                                                  lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "linear",
                                                  trend_param = "add",
                                                  endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(5)
results_bin_lin_add_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_add_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_lin_add_alpha_OR1_eq, "results/results_bin_lin_add_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_lin_mult_alpha_OR1_eq <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 1.8,
                                                   OR2 = 1,
                                                   lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "linear",
                                                   trend_param = "mult",
                                                   endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(6)
results_bin_lin_mult_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_OR1_eq, "results/results_bin_lin_mult_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - additive case - power 

scenarios_bin_lin_add_pow_OR1_eq <- data.frame(K = 1,
                                                p0 = 0.7,
                                                sigma = 1,
                                                OR1 = 1.8,
                                                OR2 = 1.8,
                                                lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                n0 = 250,
                                                n1 = 250,
                                                n01 = 250*0.5,
                                                n02 = 250-250*0.5,
                                                n11 = 250*0.5,
                                                n12 = 250-250*0.5,
                                                n22 = 250,
                                                trend = "linear",
                                                trend_param = "add",
                                                endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(7)
results_bin_lin_add_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_add_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_lin_add_pow_OR1_eq, "results/results_bin_lin_add_pow_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - power 

scenarios_bin_lin_mult_pow_OR1_eq <- data.frame(K = 1,
                                                 p0 = 0.7,
                                                 sigma = 1,
                                                 OR1 = 1.8,
                                                 OR2 = 1.8,
                                                 lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                 n0 = 250,
                                                 n1 = 250,
                                                 n01 = 250*0.5,
                                                 n02 = 250-250*0.5,
                                                 n11 = 250*0.5,
                                                 n12 = 250-250*0.5,
                                                 n22 = 250,
                                                 trend = "linear",
                                                 trend_param = "mult",
                                                 endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(8)
results_bin_lin_mult_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_OR1_eq, "results/results_bin_lin_mult_pow_OR1_eq.csv")
gc()

############################################################################################################









############################################################################################################
# BINARY CASE - linear trend - OR1<1 - different time trends, negative trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_lin_mult_alpha_diff_neg <- data.frame(K = 1,
                                                     p0 = 0.7,
                                                     sigma = 1,
                                                     OR1 = 0.4,
                                                     OR2 = 1,
                                                     lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                     lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                     lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                     n0 = 250,
                                                     n1 = 250,
                                                     n01 = 250*0.5,
                                                     n02 = 250-250*0.5,
                                                     n11 = 250*0.5,
                                                     n12 = 250-250*0.5,
                                                     n22 = 250,
                                                     trend = "linear",
                                                     trend_param = "mult",
                                                     endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(9)
results_bin_lin_mult_alpha_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_diff_neg, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_diff_neg, "results/results_bin_lin_mult_alpha_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_lin_mult_pow_diff_neg <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 0.4,
                                                   OR2 = 1.8,
                                                   lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "linear",
                                                   trend_param = "mult",
                                                   endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(10)
results_bin_lin_mult_pow_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_diff_neg, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_diff_neg, "results/results_bin_lin_mult_pow_diff_neg.csv")
gc()




############################################################################################################
# BINARY CASE - linear trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_lin_mult_alpha_OR1_diff_neg <- data.frame(K = 1,
                                                         p0 = 0.7,
                                                         sigma = 1,
                                                         OR1 = 1.8,
                                                         OR2 = 1,
                                                         lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                         lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                         lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                         n0 = 250,
                                                         n1 = 250,
                                                         n01 = 250*0.5,
                                                         n02 = 250-250*0.5,
                                                         n11 = 250*0.5,
                                                         n12 = 250-250*0.5,
                                                         n22 = 250,
                                                         trend = "linear",
                                                         trend_param = "mult",
                                                         endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(11)
results_bin_lin_mult_alpha_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_OR1_diff_neg, "results/results_bin_lin_mult_alpha_OR1_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_lin_mult_pow_OR1_diff_neg <- data.frame(K = 1,
                                                       p0 = 0.7,
                                                       sigma = 1,
                                                       OR1 = 1.8,
                                                       OR2 = 1.8,
                                                       lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                       lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                       lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                       n0 = 250,
                                                       n1 = 250,
                                                       n01 = 250*0.5,
                                                       n02 = 250-250*0.5,
                                                       n11 = 250*0.5,
                                                       n12 = 250-250*0.5,
                                                       n22 = 250,
                                                       trend = "linear",
                                                       trend_param = "mult",
                                                       endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(12)
results_bin_lin_mult_pow_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_OR1_diff_neg, "results/results_bin_lin_mult_pow_OR1_diff_neg.csv")
gc()

############################################################################################################














############################################################################################################
# BINARY CASE - linear trend - OR1<1 - different time trends, positive trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_lin_mult_alpha_diff_pos <- data.frame(K = 1,
                                                     p0 = 0.7,
                                                     sigma = 1,
                                                     OR1 = 0.4,
                                                     OR2 = 1,
                                                     lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                     lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                     lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                     n0 = 250,
                                                     n1 = 250,
                                                     n01 = 250*0.5,
                                                     n02 = 250-250*0.5,
                                                     n11 = 250*0.5,
                                                     n12 = 250-250*0.5,
                                                     n22 = 250,
                                                     trend = "linear",
                                                     trend_param = "mult",
                                                     endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(13)
results_bin_lin_mult_alpha_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_diff_pos, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_diff_pos, "results/results_bin_lin_mult_alpha_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_lin_mult_pow_diff_pos <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 0.4,
                                                   OR2 = 1.8,
                                                   lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "linear",
                                                   trend_param = "mult",
                                                   endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(14)
results_bin_lin_mult_pow_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_diff_pos, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_diff_pos, "results/results_bin_lin_mult_pow_diff_pos.csv")
gc()




############################################################################################################
# BINARY CASE - linear trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_lin_mult_alpha_OR1_diff_pos <- data.frame(K = 1,
                                                         p0 = 0.7,
                                                         sigma = 1,
                                                         OR1 = 1.8,
                                                         OR2 = 1,
                                                         lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                         lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                         lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                         n0 = 250,
                                                         n1 = 250,
                                                         n01 = 250*0.5,
                                                         n02 = 250-250*0.5,
                                                         n11 = 250*0.5,
                                                         n12 = 250-250*0.5,
                                                         n22 = 250,
                                                         trend = "linear",
                                                         trend_param = "mult",
                                                         endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(15)
results_bin_lin_mult_alpha_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_alpha_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_lin_mult_alpha_OR1_diff_pos, "results/results_bin_lin_mult_alpha_OR1_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_lin_mult_pow_OR1_diff_pos <- data.frame(K = 1,
                                                       p0 = 0.7,
                                                       sigma = 1,
                                                       OR1 = 1.8,
                                                       OR2 = 1.8,
                                                       lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                       lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                       lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                       n0 = 250,
                                                       n1 = 250,
                                                       n01 = 250*0.5,
                                                       n02 = 250-250*0.5,
                                                       n11 = 250*0.5,
                                                       n12 = 250-250*0.5,
                                                       n22 = 250,
                                                       trend = "linear",
                                                       trend_param = "mult",
                                                       endpoint = "binary") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(16)
results_bin_lin_mult_pow_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_lin_mult_pow_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_lin_mult_pow_OR1_diff_pos, "results/results_bin_lin_mult_pow_OR1_diff_pos.csv")
gc()

############################################################################################################


























############################################################################################################
# BINARY CASE - inverted U-1 trend - OR1<1 - equal time trends
############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - type one error

scenarios_bin_inv_1_add_alpha_eq <- data.frame(K = 1,
                                             p0 = 0.7,
                                             sigma = 1,
                                             OR1 = 0.4,
                                             OR2 = 1,
                                             lambda0 = seq(-0.5, 0.5, length.out = 19),
                                             lambda1 = seq(-0.5, 0.5, length.out = 19),
                                             lambda2 = seq(-0.5, 0.5, length.out = 19),
                                             n0 = 250,
                                             n1 = 250,
                                             n01 = 250*0.5,
                                             n02 = 250-250*0.5,
                                             n11 = 250*0.5,
                                             n12 = 250-250*0.5,
                                             n22 = 250,
                                             trend = "inv_u",
                                             trend_param = "add",
                                             endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(1)
results_bin_inv_1_add_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_add_alpha_eq, endpoint = "binary")

write_csv(results_bin_inv_1_add_alpha_eq, "results/results_bin_inv_1_add_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_1_mult_alpha_eq <- data.frame(K = 1,
                                              p0 = 0.7,
                                              sigma = 1,
                                              OR1 = 0.4,
                                              OR2 = 1,
                                              lambda0 = seq(-0.5, 0.5, length.out = 19),
                                              lambda1 = seq(-0.5, 0.5, length.out = 19),
                                              lambda2 = seq(-0.5, 0.5, length.out = 19),
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "inv_u",
                                              trend_param = "mult",
                                              endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(2)
results_bin_inv_1_mult_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_alpha_eq, endpoint = "binary")

write_csv(results_bin_inv_1_mult_alpha_eq, "results/results_bin_inv_1_mult_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - power

scenarios_bin_inv_1_add_pow_eq <- data.frame(K = 1,
                                           p0 = 0.7,
                                           sigma = 1,
                                           OR1 = 0.4,
                                           OR2 = 1.8,
                                           lambda0 = seq(-0.5, 0.5, length.out = 19),
                                           lambda1 = seq(-0.5, 0.5, length.out = 19),
                                           lambda2 = seq(-0.5, 0.5, length.out = 19),
                                           n0 = 250,
                                           n1 = 250,
                                           n01 = 250*0.5,
                                           n02 = 250-250*0.5,
                                           n11 = 250*0.5,
                                           n12 = 250-250*0.5,
                                           n22 = 250,
                                           trend = "inv_u",
                                           trend_param = "add",
                                           endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(3)
results_bin_inv_1_add_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_add_pow_eq, endpoint = "binary")

write_csv(results_bin_inv_1_add_pow_eq, "results/results_bin_inv_1_add_pow_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_1_mult_pow_eq <- data.frame(K = 1,
                                            p0 = 0.7,
                                            sigma = 1,
                                            OR1 = 0.4,
                                            OR2 = 1.8,
                                            lambda0 = seq(-0.5, 0.5, length.out = 19),
                                            lambda1 = seq(-0.5, 0.5, length.out = 19),
                                            lambda2 = seq(-0.5, 0.5, length.out = 19),
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            trend_param = "mult",
                                            endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(4)
results_bin_inv_1_mult_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_pow_eq, endpoint = "binary")

write_csv(results_bin_inv_1_mult_pow_eq, "results/results_bin_inv_1_mult_pow_eq.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-1 trend - OR1>1 - equal time trends
############################################################################################################


# Binary scenario (equal time trends) (OR>1) - additive case - type one error

scenarios_bin_inv_1_add_alpha_OR1_eq <- data.frame(K = 1,
                                                 p0 = 0.7,
                                                 sigma = 1,
                                                 OR1 = 1.8,
                                                 OR2 = 1,
                                                 lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                 n0 = 250,
                                                 n1 = 250,
                                                 n01 = 250*0.5,
                                                 n02 = 250-250*0.5,
                                                 n11 = 250*0.5,
                                                 n12 = 250-250*0.5,
                                                 n22 = 250,
                                                 trend = "inv_u",
                                                 trend_param = "add",
                                                 endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(5)
results_bin_inv_1_add_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_add_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_1_add_alpha_OR1_eq, "results/results_bin_inv_1_add_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_1_mult_alpha_OR1_eq <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 1.8,
                                                  OR2 = 1,
                                                  lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "inv_u",
                                                  trend_param = "mult",
                                                  endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(6)
results_bin_inv_1_mult_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_1_mult_alpha_OR1_eq, "results/results_bin_inv_1_mult_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - additive case - power 

scenarios_bin_inv_1_add_pow_OR1_eq <- data.frame(K = 1,
                                               p0 = 0.7,
                                               sigma = 1,
                                               OR1 = 1.8,
                                               OR2 = 1.8,
                                               lambda0 = seq(-0.5, 0.5, length.out = 19),
                                               lambda1 = seq(-0.5, 0.5, length.out = 19),
                                               lambda2 = seq(-0.5, 0.5, length.out = 19),
                                               n0 = 250,
                                               n1 = 250,
                                               n01 = 250*0.5,
                                               n02 = 250-250*0.5,
                                               n11 = 250*0.5,
                                               n12 = 250-250*0.5,
                                               n22 = 250,
                                               trend = "inv_u",
                                               trend_param = "add",
                                               endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(7)
results_bin_inv_1_add_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_add_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_1_add_pow_OR1_eq, "results/results_bin_inv_1_add_pow_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_1_mult_pow_OR1_eq <- data.frame(K = 1,
                                                p0 = 0.7,
                                                sigma = 1,
                                                OR1 = 1.8,
                                                OR2 = 1.8,
                                                lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                n0 = 250,
                                                n1 = 250,
                                                n01 = 250*0.5,
                                                n02 = 250-250*0.5,
                                                n11 = 250*0.5,
                                                n12 = 250-250*0.5,
                                                n22 = 250,
                                                trend = "inv_u",
                                                trend_param = "mult",
                                                endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(8)
results_bin_inv_1_mult_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_1_mult_pow_OR1_eq, "results/results_bin_inv_1_mult_pow_OR1_eq.csv")
gc()

############################################################################################################









############################################################################################################
# BINARY CASE - inverted U-1 trend - OR1<1 - different time trends, negative trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_1_mult_alpha_diff_neg <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 0.4,
                                                    OR2 = 1,
                                                    lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(9)
results_bin_inv_1_mult_alpha_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_alpha_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_1_mult_alpha_diff_neg, "results/results_bin_inv_1_mult_alpha_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_1_mult_pow_diff_neg <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 0.4,
                                                  OR2 = 1.8,
                                                  lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "inv_u",
                                                  trend_param = "mult",
                                                  endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(10)
results_bin_inv_1_mult_pow_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_pow_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_1_mult_pow_diff_neg, "results/results_bin_inv_1_mult_pow_diff_neg.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-1 trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_1_mult_alpha_OR1_diff_neg <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1,
                                                        lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        n0 = 250,
                                                        n1 = 250,
                                                        n01 = 250*0.5,
                                                        n02 = 250-250*0.5,
                                                        n11 = 250*0.5,
                                                        n12 = 250-250*0.5,
                                                        n22 = 250,
                                                        trend = "inv_u",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(11)
results_bin_inv_1_mult_alpha_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_alpha_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_1_mult_alpha_OR1_diff_neg, "results/results_bin_inv_1_mult_alpha_OR1_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_1_mult_pow_OR1_diff_neg <- data.frame(K = 1,
                                                      p0 = 0.7,
                                                      sigma = 1,
                                                      OR1 = 1.8,
                                                      OR2 = 1.8,
                                                      lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                      lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                      lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                      n0 = 250,
                                                      n1 = 250,
                                                      n01 = 250*0.5,
                                                      n02 = 250-250*0.5,
                                                      n11 = 250*0.5,
                                                      n12 = 250-250*0.5,
                                                      n22 = 250,
                                                      trend = "inv_u",
                                                      trend_param = "mult",
                                                      endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(12)
results_bin_inv_1_mult_pow_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_pow_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_1_mult_pow_OR1_diff_neg, "results/results_bin_inv_1_mult_pow_OR1_diff_neg.csv")
gc()

############################################################################################################














############################################################################################################
# BINARY CASE - inverted U-1 trend - OR1<1 - different time trends, positive trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_1_mult_alpha_diff_pos <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 0.4,
                                                    OR2 = 1,
                                                    lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(13)
results_bin_inv_1_mult_alpha_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_alpha_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_1_mult_alpha_diff_pos, "results/results_bin_inv_1_mult_alpha_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_1_mult_pow_diff_pos <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 0.4,
                                                  OR2 = 1.8,
                                                  lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "inv_u",
                                                  trend_param = "mult",
                                                  endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(14)
results_bin_inv_1_mult_pow_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_pow_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_1_mult_pow_diff_pos, "results/results_bin_inv_1_mult_pow_diff_pos.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-1 trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_1_mult_alpha_OR1_diff_pos <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1,
                                                        lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        n0 = 250,
                                                        n1 = 250,
                                                        n01 = 250*0.5,
                                                        n02 = 250-250*0.5,
                                                        n11 = 250*0.5,
                                                        n12 = 250-250*0.5,
                                                        n22 = 250,
                                                        trend = "inv_u",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(15)
results_bin_inv_1_mult_alpha_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_alpha_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_1_mult_alpha_OR1_diff_pos, "results/results_bin_inv_1_mult_alpha_OR1_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_1_mult_pow_OR1_diff_pos <- data.frame(K = 1,
                                                      p0 = 0.7,
                                                      sigma = 1,
                                                      OR1 = 1.8,
                                                      OR2 = 1.8,
                                                      lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                      lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                      lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                      n0 = 250,
                                                      n1 = 250,
                                                      n01 = 250*0.5,
                                                      n02 = 250-250*0.5,
                                                      n11 = 250*0.5,
                                                      n12 = 250-250*0.5,
                                                      n22 = 250,
                                                      trend = "inv_u",
                                                      trend_param = "mult",
                                                      endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(16)
results_bin_inv_1_mult_pow_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_1_mult_pow_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_1_mult_pow_OR1_diff_pos, "results/results_bin_inv_1_mult_pow_OR1_diff_pos.csv")
gc()

############################################################################################################





























############################################################################################################
# BINARY CASE - inverted U-2 trend - OR1<1 - equal time trends
############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - type one error

scenarios_bin_inv_2_add_alpha_eq <- data.frame(K = 1,
                                               p0 = 0.7,
                                               sigma = 1,
                                               OR1 = 0.4,
                                               OR2 = 1,
                                               lambda0 = seq(-0.5, 0.5, length.out = 19),
                                               lambda1 = seq(-0.5, 0.5, length.out = 19),
                                               lambda2 = seq(-0.5, 0.5, length.out = 19),
                                               n0 = 250,
                                               n1 = 250,
                                               n01 = 250*0.5,
                                               n02 = 250-250*0.5,
                                               n11 = 250*0.5,
                                               n12 = 250-250*0.5,
                                               n22 = 250,
                                               trend = "inv_u",
                                               trend_param = "add",
                                               endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(1)
results_bin_inv_2_add_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_add_alpha_eq, endpoint = "binary")

write_csv(results_bin_inv_2_add_alpha_eq, "results/results_bin_inv_2_add_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_2_mult_alpha_eq <- data.frame(K = 1,
                                                p0 = 0.7,
                                                sigma = 1,
                                                OR1 = 0.4,
                                                OR2 = 1,
                                                lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                n0 = 250,
                                                n1 = 250,
                                                n01 = 250*0.5,
                                                n02 = 250-250*0.5,
                                                n11 = 250*0.5,
                                                n12 = 250-250*0.5,
                                                n22 = 250,
                                                trend = "inv_u",
                                                trend_param = "mult",
                                                endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(2)
results_bin_inv_2_mult_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_alpha_eq, endpoint = "binary")

write_csv(results_bin_inv_2_mult_alpha_eq, "results/results_bin_inv_2_mult_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - power

scenarios_bin_inv_2_add_pow_eq <- data.frame(K = 1,
                                             p0 = 0.7,
                                             sigma = 1,
                                             OR1 = 0.4,
                                             OR2 = 1.8,
                                             lambda0 = seq(-0.5, 0.5, length.out = 19),
                                             lambda1 = seq(-0.5, 0.5, length.out = 19),
                                             lambda2 = seq(-0.5, 0.5, length.out = 19),
                                             n0 = 250,
                                             n1 = 250,
                                             n01 = 250*0.5,
                                             n02 = 250-250*0.5,
                                             n11 = 250*0.5,
                                             n12 = 250-250*0.5,
                                             n22 = 250,
                                             trend = "inv_u",
                                             trend_param = "add",
                                             endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(3)
results_bin_inv_2_add_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_add_pow_eq, endpoint = "binary")

write_csv(results_bin_inv_2_add_pow_eq, "results/results_bin_inv_2_add_pow_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_2_mult_pow_eq <- data.frame(K = 1,
                                              p0 = 0.7,
                                              sigma = 1,
                                              OR1 = 0.4,
                                              OR2 = 1.8,
                                              lambda0 = seq(-0.5, 0.5, length.out = 19),
                                              lambda1 = seq(-0.5, 0.5, length.out = 19),
                                              lambda2 = seq(-0.5, 0.5, length.out = 19),
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "inv_u",
                                              trend_param = "mult",
                                              endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(4)
results_bin_inv_2_mult_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_pow_eq, endpoint = "binary")

write_csv(results_bin_inv_2_mult_pow_eq, "results/results_bin_inv_2_mult_pow_eq.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-2 trend - OR1>1 - equal time trends
############################################################################################################


# Binary scenario (equal time trends) (OR>1) - additive case - type one error

scenarios_bin_inv_2_add_alpha_OR1_eq <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 1.8,
                                                   OR2 = 1,
                                                   lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "inv_u",
                                                   trend_param = "add",
                                                   endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(5)
results_bin_inv_2_add_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_add_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_2_add_alpha_OR1_eq, "results/results_bin_inv_2_add_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_2_mult_alpha_OR1_eq <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 1.8,
                                                    OR2 = 1,
                                                    lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(6)
results_bin_inv_2_mult_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_2_mult_alpha_OR1_eq, "results/results_bin_inv_2_mult_alpha_OR1_eq.csv")
gc()
############################################################################################################

# Binary scenario (equal time trends) (OR>1) - additive case - power 

scenarios_bin_inv_2_add_pow_OR1_eq <- data.frame(K = 1,
                                                 p0 = 0.7,
                                                 sigma = 1,
                                                 OR1 = 1.8,
                                                 OR2 = 1.8,
                                                 lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                 n0 = 250,
                                                 n1 = 250,
                                                 n01 = 250*0.5,
                                                 n02 = 250-250*0.5,
                                                 n11 = 250*0.5,
                                                 n12 = 250-250*0.5,
                                                 n22 = 250,
                                                 trend = "inv_u",
                                                 trend_param = "add",
                                                 endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(7)
results_bin_inv_2_add_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_add_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_2_add_pow_OR1_eq, "results/results_bin_inv_2_add_pow_OR1_eq.csv")

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_2_mult_pow_OR1_eq <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 1.8,
                                                  OR2 = 1.8,
                                                  lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "inv_u",
                                                  trend_param = "mult",
                                                  endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(8)
results_bin_inv_2_mult_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_2_mult_pow_OR1_eq, "results/results_bin_inv_2_mult_pow_OR1_eq.csv")
gc()

############################################################################################################









############################################################################################################
# BINARY CASE - inverted U-2 trend - OR1<1 - different time trends, negative trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_2_mult_alpha_diff_neg <- data.frame(K = 1,
                                                      p0 = 0.7,
                                                      sigma = 1,
                                                      OR1 = 0.4,
                                                      OR2 = 1,
                                                      lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                      lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                      lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                      n0 = 250,
                                                      n1 = 250,
                                                      n01 = 250*0.5,
                                                      n02 = 250-250*0.5,
                                                      n11 = 250*0.5,
                                                      n12 = 250-250*0.5,
                                                      n22 = 250,
                                                      trend = "inv_u",
                                                      trend_param = "mult",
                                                      endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(9)
results_bin_inv_2_mult_alpha_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_alpha_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_2_mult_alpha_diff_neg, "results/results_bin_inv_2_mult_alpha_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_2_mult_pow_diff_neg <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 0.4,
                                                    OR2 = 1.8,
                                                    lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(10)
results_bin_inv_2_mult_pow_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_pow_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_2_mult_pow_diff_neg, "results/results_bin_inv_2_mult_pow_diff_neg.csv")
gc()



############################################################################################################
# BINARY CASE - inverted U-2 trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_2_mult_alpha_OR1_diff_neg <- data.frame(K = 1,
                                                          p0 = 0.7,
                                                          sigma = 1,
                                                          OR1 = 1.8,
                                                          OR2 = 1,
                                                          lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                          lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                          lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                          n0 = 250,
                                                          n1 = 250,
                                                          n01 = 250*0.5,
                                                          n02 = 250-250*0.5,
                                                          n11 = 250*0.5,
                                                          n12 = 250-250*0.5,
                                                          n22 = 250,
                                                          trend = "inv_u",
                                                          trend_param = "mult",
                                                          endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(11)
results_bin_inv_2_mult_alpha_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_alpha_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_2_mult_alpha_OR1_diff_neg, "results/results_bin_inv_2_mult_alpha_OR1_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_2_mult_pow_OR1_diff_neg <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1.8,
                                                        lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        n0 = 250,
                                                        n1 = 250,
                                                        n01 = 250*0.5,
                                                        n02 = 250-250*0.5,
                                                        n11 = 250*0.5,
                                                        n12 = 250-250*0.5,
                                                        n22 = 250,
                                                        trend = "inv_u",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(12)
results_bin_inv_2_mult_pow_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_pow_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_2_mult_pow_OR1_diff_neg, "results/results_bin_inv_2_mult_pow_OR1_diff_neg.csv")
gc()

############################################################################################################














############################################################################################################
# BINARY CASE - inverted U-2 trend - OR1<1 - different time trends, positive trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_2_mult_alpha_diff_pos <- data.frame(K = 1,
                                                      p0 = 0.7,
                                                      sigma = 1,
                                                      OR1 = 0.4,
                                                      OR2 = 1,
                                                      lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                      lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                      lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                      n0 = 250,
                                                      n1 = 250,
                                                      n01 = 250*0.5,
                                                      n02 = 250-250*0.5,
                                                      n11 = 250*0.5,
                                                      n12 = 250-250*0.5,
                                                      n22 = 250,
                                                      trend = "inv_u",
                                                      trend_param = "mult",
                                                      endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(13)
results_bin_inv_2_mult_alpha_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_alpha_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_2_mult_alpha_diff_pos, "results/results_bin_inv_2_mult_alpha_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_2_mult_pow_diff_pos <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 0.4,
                                                    OR2 = 1.8,
                                                    lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(14)
results_bin_inv_2_mult_pow_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_pow_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_2_mult_pow_diff_pos, "results/results_bin_inv_2_mult_pow_diff_pos.csv")
gc()



############################################################################################################
# BINARY CASE - inverted U-2 trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_2_mult_alpha_OR1_diff_pos <- data.frame(K = 1,
                                                          p0 = 0.7,
                                                          sigma = 1,
                                                          OR1 = 1.8,
                                                          OR2 = 1,
                                                          lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                          lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                          lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                          n0 = 250,
                                                          n1 = 250,
                                                          n01 = 250*0.5,
                                                          n02 = 250-250*0.5,
                                                          n11 = 250*0.5,
                                                          n12 = 250-250*0.5,
                                                          n22 = 250,
                                                          trend = "inv_u",
                                                          trend_param = "mult",
                                                          endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(15)
results_bin_inv_2_mult_alpha_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_alpha_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_2_mult_alpha_OR1_diff_pos, "results/results_bin_inv_2_mult_alpha_OR1_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_2_mult_pow_OR1_diff_pos <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1.8,
                                                        lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        n0 = 250,
                                                        n1 = 250,
                                                        n01 = 250*0.5,
                                                        n02 = 250-250*0.5,
                                                        n11 = 250*0.5,
                                                        n12 = 250-250*0.5,
                                                        n22 = 250,
                                                        trend = "inv_u",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(16)
results_bin_inv_2_mult_pow_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_2_mult_pow_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_2_mult_pow_OR1_diff_pos, "results/results_bin_inv_2_mult_pow_OR1_diff_pos.csv")
gc()

############################################################################################################












############################################################################################################
# BINARY CASE - inverted U-3 trend - OR1<1 - equal time trends
############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - type one error

scenarios_bin_inv_3_add_alpha_eq <- data.frame(K = 1,
                                               p0 = 0.7,
                                               sigma = 1,
                                               OR1 = 0.4,
                                               OR2 = 1,
                                               lambda0 = seq(-0.5, 0.5, length.out = 19),
                                               lambda1 = seq(-0.5, 0.5, length.out = 19),
                                               lambda2 = seq(-0.5, 0.5, length.out = 19),
                                               n0 = 250,
                                               n1 = 250,
                                               n01 = 250*0.5,
                                               n02 = 250-250*0.5,
                                               n11 = 250*0.5,
                                               n12 = 250-250*0.5,
                                               n22 = 250,
                                               trend = "inv_u",
                                               trend_param = "add",
                                               endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(1)
results_bin_inv_3_add_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_add_alpha_eq, endpoint = "binary")

write_csv(results_bin_inv_3_add_alpha_eq, "results/results_bin_inv_3_add_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_3_mult_alpha_eq <- data.frame(K = 1,
                                                p0 = 0.7,
                                                sigma = 1,
                                                OR1 = 0.4,
                                                OR2 = 1,
                                                lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                n0 = 250,
                                                n1 = 250,
                                                n01 = 250*0.5,
                                                n02 = 250-250*0.5,
                                                n11 = 250*0.5,
                                                n12 = 250-250*0.5,
                                                n22 = 250,
                                                trend = "inv_u",
                                                trend_param = "mult",
                                                endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(2)
results_bin_inv_3_mult_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_alpha_eq, endpoint = "binary")

write_csv(results_bin_inv_3_mult_alpha_eq, "results/results_bin_inv_3_mult_alpha_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - additive case - power

scenarios_bin_inv_3_add_pow_eq <- data.frame(K = 1,
                                             p0 = 0.7,
                                             sigma = 1,
                                             OR1 = 0.4,
                                             OR2 = 1.8,
                                             lambda0 = seq(-0.5, 0.5, length.out = 19),
                                             lambda1 = seq(-0.5, 0.5, length.out = 19),
                                             lambda2 = seq(-0.5, 0.5, length.out = 19),
                                             n0 = 250,
                                             n1 = 250,
                                             n01 = 250*0.5,
                                             n02 = 250-250*0.5,
                                             n11 = 250*0.5,
                                             n12 = 250-250*0.5,
                                             n22 = 250,
                                             trend = "inv_u",
                                             trend_param = "add",
                                             endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(3)
results_bin_inv_3_add_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_add_pow_eq, endpoint = "binary")

write_csv(results_bin_inv_3_add_pow_eq, "results/results_bin_inv_3_add_pow_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_3_mult_pow_eq <- data.frame(K = 1,
                                              p0 = 0.7,
                                              sigma = 1,
                                              OR1 = 0.4,
                                              OR2 = 1.8,
                                              lambda0 = seq(-0.5, 0.5, length.out = 19),
                                              lambda1 = seq(-0.5, 0.5, length.out = 19),
                                              lambda2 = seq(-0.5, 0.5, length.out = 19),
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "inv_u",
                                              trend_param = "mult",
                                              endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(4)
results_bin_inv_3_mult_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_pow_eq, endpoint = "binary")

write_csv(results_bin_inv_3_mult_pow_eq, "results/results_bin_inv_3_mult_pow_eq.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-3 trend - OR1>1 - equal time trends
############################################################################################################


# Binary scenario (equal time trends) (OR>1) - additive case - type one error

scenarios_bin_inv_3_add_alpha_OR1_eq <- data.frame(K = 1,
                                                   p0 = 0.7,
                                                   sigma = 1,
                                                   OR1 = 1.8,
                                                   OR2 = 1,
                                                   lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                   lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                   n0 = 250,
                                                   n1 = 250,
                                                   n01 = 250*0.5,
                                                   n02 = 250-250*0.5,
                                                   n11 = 250*0.5,
                                                   n12 = 250-250*0.5,
                                                   n22 = 250,
                                                   trend = "inv_u",
                                                   trend_param = "add",
                                                   endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(5)
results_bin_inv_3_add_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_add_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_3_add_alpha_OR1_eq, "results/results_bin_inv_3_add_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_3_mult_alpha_OR1_eq <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 1.8,
                                                    OR2 = 1,
                                                    lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(6)
results_bin_inv_3_mult_alpha_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_alpha_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_3_mult_alpha_OR1_eq, "results/results_bin_inv_3_mult_alpha_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - additive case - power 

scenarios_bin_inv_3_add_pow_OR1_eq <- data.frame(K = 1,
                                                 p0 = 0.7,
                                                 sigma = 1,
                                                 OR1 = 1.8,
                                                 OR2 = 1.8,
                                                 lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                 lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                 n0 = 250,
                                                 n1 = 250,
                                                 n01 = 250*0.5,
                                                 n02 = 250-250*0.5,
                                                 n11 = 250*0.5,
                                                 n12 = 250-250*0.5,
                                                 n22 = 250,
                                                 trend = "inv_u",
                                                 trend_param = "add",
                                                 endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(7)
results_bin_inv_3_add_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_add_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_3_add_pow_OR1_eq, "results/results_bin_inv_3_add_pow_OR1_eq.csv")
gc()

############################################################################################################

# Binary scenario (equal time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_3_mult_pow_OR1_eq <- data.frame(K = 1,
                                                  p0 = 0.7,
                                                  sigma = 1,
                                                  OR1 = 1.8,
                                                  OR2 = 1.8,
                                                  lambda0 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                  lambda2 = seq(-0.5, 0.5, length.out = 19),
                                                  n0 = 250,
                                                  n1 = 250,
                                                  n01 = 250*0.5,
                                                  n02 = 250-250*0.5,
                                                  n11 = 250*0.5,
                                                  n12 = 250-250*0.5,
                                                  n22 = 250,
                                                  trend = "inv_u",
                                                  trend_param = "mult",
                                                  endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(8)
results_bin_inv_3_mult_pow_OR1_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_pow_OR1_eq, endpoint = "binary")

write_csv(results_bin_inv_3_mult_pow_OR1_eq, "results/results_bin_inv_3_mult_pow_OR1_eq.csv")
gc()

############################################################################################################









############################################################################################################
# BINARY CASE - inverted U-3 trend - OR1<1 - different time trends, negative trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_3_mult_alpha_diff_neg <- data.frame(K = 1,
                                                      p0 = 0.7,
                                                      sigma = 1,
                                                      OR1 = 0.4,
                                                      OR2 = 1,
                                                      lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                      lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                      lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                      n0 = 250,
                                                      n1 = 250,
                                                      n01 = 250*0.5,
                                                      n02 = 250-250*0.5,
                                                      n11 = 250*0.5,
                                                      n12 = 250-250*0.5,
                                                      n22 = 250,
                                                      trend = "inv_u",
                                                      trend_param = "mult",
                                                      endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(9)
results_bin_inv_3_mult_alpha_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_alpha_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_3_mult_alpha_diff_neg, "results/results_bin_inv_3_mult_alpha_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_3_mult_pow_diff_neg <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 0.4,
                                                    OR2 = 1.8,
                                                    lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(10)
results_bin_inv_3_mult_pow_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_pow_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_3_mult_pow_diff_neg, "results/results_bin_inv_3_mult_pow_diff_neg.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-3 trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_3_mult_alpha_OR1_diff_neg <- data.frame(K = 1,
                                                          p0 = 0.7,
                                                          sigma = 1,
                                                          OR1 = 1.8,
                                                          OR2 = 1,
                                                          lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                          lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                          lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                          n0 = 250,
                                                          n1 = 250,
                                                          n01 = 250*0.5,
                                                          n02 = 250-250*0.5,
                                                          n11 = 250*0.5,
                                                          n12 = 250-250*0.5,
                                                          n22 = 250,
                                                          trend = "inv_u",
                                                          trend_param = "mult",
                                                          endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(11)
results_bin_inv_3_mult_alpha_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_alpha_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_3_mult_alpha_OR1_diff_neg, "results/results_bin_inv_3_mult_alpha_OR1_diff_neg.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_3_mult_pow_OR1_diff_neg <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1.8,
                                                        lambda0 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = -0.2282587, #-log(0.7/(1-0.7))-log((1/0.65)-1)
                                                        n0 = 250,
                                                        n1 = 250,
                                                        n01 = 250*0.5,
                                                        n02 = 250-250*0.5,
                                                        n11 = 250*0.5,
                                                        n12 = 250-250*0.5,
                                                        n22 = 250,
                                                        trend = "inv_u",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(12)
results_bin_inv_3_mult_pow_OR1_diff_neg <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_pow_OR1_diff_neg, endpoint = "binary")

write_csv(results_bin_inv_3_mult_pow_OR1_diff_neg, "results/results_bin_inv_3_mult_pow_OR1_diff_neg.csv")
gc()

############################################################################################################














############################################################################################################
# BINARY CASE - inverted U-3 trend - OR1<1 - different time trends, positive trends for C and T2
############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - type one error

scenarios_bin_inv_3_mult_alpha_diff_pos <- data.frame(K = 1,
                                                      p0 = 0.7,
                                                      sigma = 1,
                                                      OR1 = 0.4,
                                                      OR2 = 1,
                                                      lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                      lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                      lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                      n0 = 250,
                                                      n1 = 250,
                                                      n01 = 250*0.5,
                                                      n02 = 250-250*0.5,
                                                      n11 = 250*0.5,
                                                      n12 = 250-250*0.5,
                                                      n22 = 250,
                                                      trend = "inv_u",
                                                      trend_param = "mult",
                                                      endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)



set.seed(13)
results_bin_inv_3_mult_alpha_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_alpha_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_3_mult_alpha_diff_pos, "results/results_bin_inv_3_mult_alpha_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR1<1) - multiplicative case - power

scenarios_bin_inv_3_mult_pow_diff_pos <- data.frame(K = 1,
                                                    p0 = 0.7,
                                                    sigma = 1,
                                                    OR1 = 0.4,
                                                    OR2 = 1.8,
                                                    lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                    lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                    lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                    n0 = 250,
                                                    n1 = 250,
                                                    n01 = 250*0.5,
                                                    n02 = 250-250*0.5,
                                                    n11 = 250*0.5,
                                                    n12 = 250-250*0.5,
                                                    n22 = 250,
                                                    trend = "inv_u",
                                                    trend_param = "mult",
                                                    endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0", "H1"),
         timing = n01/n0)


set.seed(14)
results_bin_inv_3_mult_pow_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_pow_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_3_mult_pow_diff_pos, "results/results_bin_inv_3_mult_pow_diff_pos.csv")
gc()




############################################################################################################
# BINARY CASE - inverted U-3 trend - OR1>1 - different time trends
############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - type one error 

scenarios_bin_inv_3_mult_alpha_OR1_diff_pos <- data.frame(K = 1,
                                                          p0 = 0.7,
                                                          sigma = 1,
                                                          OR1 = 1.8,
                                                          OR2 = 1,
                                                          lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                          lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                          lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                          n0 = 250,
                                                          n1 = 250,
                                                          n01 = 250*0.5,
                                                          n02 = 250-250*0.5,
                                                          n11 = 250*0.5,
                                                          n12 = 250-250*0.5,
                                                          n22 = 250,
                                                          trend = "inv_u",
                                                          trend_param = "mult",
                                                          endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)



set.seed(15)
results_bin_inv_3_mult_alpha_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_alpha_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_3_mult_alpha_OR1_diff_pos, "results/results_bin_inv_3_mult_alpha_OR1_diff_pos.csv")
gc()

############################################################################################################

# Binary scenario (different time trends) (OR>1) - multiplicative case - power 

scenarios_bin_inv_3_mult_pow_OR1_diff_pos <- data.frame(K = 1,
                                                        p0 = 0.7,
                                                        sigma = 1,
                                                        OR1 = 1.8,
                                                        OR2 = 1.8,
                                                        lambda0 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        lambda1 = seq(-0.5, 0.5, length.out = 19),
                                                        lambda2 = 0.2513144, #-log(0.7/(1-0.7))-log((1/0.75)-1)
                                                        n0 = 250,
                                                        n1 = 250,
                                                        n01 = 250*0.5,
                                                        n02 = 250-250*0.5,
                                                        n11 = 250*0.5,
                                                        n12 = 250-250*0.5,
                                                        n22 = 250,
                                                        trend = "inv_u",
                                                        trend_param = "mult",
                                                        endpoint = "binary") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(OR2==1, "H0_1", "H1_1"),
         timing = n01/n0)


set.seed(16)
results_bin_inv_3_mult_pow_OR1_diff_pos <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_bin_inv_3_mult_pow_OR1_diff_pos, endpoint = "binary")

write_csv(results_bin_inv_3_mult_pow_OR1_diff_pos, "results/results_bin_inv_3_mult_pow_OR1_diff_pos.csv")
gc()

############################################################################################################









############################################################################################################
# CONTINUOUS CASE - equal time trends
############################################################################################################

# Continuous case (equal time trends) - linear trend - type one error


scenarios_cont_lin_alpha_eq <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0,
                                          lambda0 = seq(-0.15, 0.15, length.out = 19),
                                          lambda1 = seq(-0.15, 0.15, length.out = 19),
                                          lambda2 = seq(-0.15, 0.15, length.out = 19),
                                          n0 = 250,
                                          n1 = 250,
                                          n01 = 250*0.5,
                                          n02 = 250-250*0.5,
                                          n11 = 250*0.5,
                                          n12 = 250-250*0.5,
                                          n22 = 250,
                                          trend = "linear",
                                          endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)


set.seed(17)
results_cont_lin_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_alpha_eq, endpoint = "continuous")

write_csv(results_cont_lin_alpha_eq, "results/results_cont_lin_alpha_eq.csv")
gc()

############################################################################################################

# Continuous case (equal time trends) - stepwise trend - type one error


scenarios_cont_step_alpha_eq <- data.frame(K = 1,
                                           mu0 = 0,
                                           sigma = 1,
                                           delta1 = 0.25,
                                           delta2 = 0,
                                           lambda0 = seq(-0.15, 0.15, length.out = 19),
                                           lambda1 = seq(-0.15, 0.15, length.out = 19),
                                           lambda2 = seq(-0.15, 0.15, length.out = 19),
                                           n0 = 250,
                                           n1 = 250,
                                           n01 = 250*0.5,
                                           n02 = 250-250*0.5,
                                           n11 = 250*0.5,
                                           n12 = 250-250*0.5,
                                           n22 = 250,
                                           trend = "stepwise",
                                           endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(18)
results_cont_step_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_step_alpha_eq, endpoint = "continuous")

write_csv(results_cont_step_alpha_eq, "results/results_cont_step_alpha_eq.csv")
gc()

############################################################################################################

# Continuous case (equal time trends) - linear trend - power

scenarios_cont_lin_pow_eq <- data.frame(K = 1,
                                        mu0 = 0,
                                        sigma = 1,
                                        delta1 = 0.25,
                                        delta2 = 0.25,
                                        lambda0 = seq(-0.15, 0.15, length.out = 19),
                                        lambda1 = seq(-0.15, 0.15, length.out = 19),
                                        lambda2 = seq(-0.15, 0.15, length.out = 19),
                                        n0 = 250,
                                        n1 = 250,
                                        n01 = 250*0.5,
                                        n02 = 250-250*0.5,
                                        n11 = 250*0.5,
                                        n12 = 250-250*0.5,
                                        n22 = 250,
                                        trend = "linear",
                                        endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(19)
results_cont_lin_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_pow_eq, endpoint = "continuous")

write_csv(results_cont_lin_pow_eq, "results/results_cont_lin_pow_eq.csv")
gc()

############################################################################################################

# Continuous case (equal time trends) - stepwise trend - power


scenarios_cont_step_pow_eq <- data.frame(K = 1,
                                         mu0 = 0,
                                         sigma = 1,
                                         delta1 = 0.25,
                                         delta2 = 0.25,
                                         lambda0 = seq(-0.15, 0.15, length.out = 19),
                                         lambda1 = seq(-0.15, 0.15, length.out = 19),
                                         lambda2 = seq(-0.15, 0.15, length.out = 19),
                                         n0 = 250,
                                         n1 = 250,
                                         n01 = 250*0.5,
                                         n02 = 250-250*0.5,
                                         n11 = 250*0.5,
                                         n12 = 250-250*0.5,
                                         n22 = 250,
                                         trend = "stepwise",
                                         endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(20)
results_cont_step_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_step_pow_eq, endpoint = "continuous")

write_csv(results_cont_step_pow_eq, "results/results_cont_step_pow_eq.csv")
gc()

############################################################################################################








############################################################################################################
# CONTINUOUS CASE - different time trends
############################################################################################################

# Continuous case (different time trends) - linear trend - type one error


scenarios_cont_lin_alpha_diff <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0,
                                            lambda0 = 0.1,
                                            lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                            lambda2 = 0.1,
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "linear",
                                            endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)


set.seed(21)
results_cont_lin_alpha_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_alpha_diff, endpoint = "continuous")

write_csv(results_cont_lin_alpha_diff, "results/results_cont_lin_alpha_diff.csv")
gc()

############################################################################################################

# Continuous case (different time trends) - stepwise trend - type one error


scenarios_cont_step_alpha_diff <- data.frame(K = 1,
                                             mu0 = 0,
                                             sigma = 1,
                                             delta1 = 0.25,
                                             delta2 = 0,
                                             lambda0 = 0.1,
                                             lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                             lambda2 = 0.1,
                                             n0 = 250,
                                             n1 = 250,
                                             n01 = 250*0.5,
                                             n02 = 250-250*0.5,
                                             n11 = 250*0.5,
                                             n12 = 250-250*0.5,
                                             n22 = 250,
                                             trend = "stepwise",
                                             endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(22)
results_cont_step_alpha_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_step_alpha_diff, endpoint = "continuous")

write_csv(results_cont_step_alpha_diff, "results/results_cont_step_alpha_diff.csv")
gc()

############################################################################################################

# Continuous case (different time trends) - linear trend - power

scenarios_cont_lin_pow_diff <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0.25,
                                          lambda0 = 0.1,
                                          lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                          lambda2 = 0.1,
                                          n0 = 250,
                                          n1 = 250,
                                          n01 = 250*0.5,
                                          n02 = 250-250*0.5,
                                          n11 = 250*0.5,
                                          n12 = 250-250*0.5,
                                          n22 = 250,
                                          trend = "linear",
                                          endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(23)
results_cont_lin_pow_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_lin_pow_diff, endpoint = "continuous")

write_csv(results_cont_lin_pow_diff, "results/results_cont_lin_pow_diff.csv")
gc()

############################################################################################################

# Continuous case (different time trends) - stepwise trend - power


scenarios_cont_step_pow_diff <- data.frame(K = 1,
                                           mu0 = 0,
                                           sigma = 1,
                                           delta1 = 0.25,
                                           delta2 = 0.25,
                                           lambda0 = 0.1,
                                           lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                           lambda2 = 0.1,
                                           n0 = 250,
                                           n1 = 250,
                                           n01 = 250*0.5,
                                           n02 = 250-250*0.5,
                                           n11 = 250*0.5,
                                           n12 = 250-250*0.5,
                                           n22 = 250,
                                           trend = "stepwise",
                                           endpoint = "continuous") %>%
  mutate(timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(24)
results_cont_step_pow_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_step_pow_diff, endpoint = "continuous")

write_csv(results_cont_step_pow_diff, "results/results_cont_step_pow_diff.csv")
gc()

############################################################################################################








############################################################################################################
# CONTINUOUS CASE - equal time trends - inverted U trend
############################################################################################################

############################################################################################################

# Continuous case (equal time trends) - inverted-U trend, peak in the middle of period 1 - type one error


scenarios_cont_inv_1_alpha_eq <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0,
                                            lambda0 = seq(-0.15, 0.15, length.out = 19),
                                            lambda1 = seq(-0.15, 0.15, length.out = 19),
                                            lambda2 = seq(-0.15, 0.15, length.out = 19),
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(25)
results_cont_inv_1_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_1_alpha_eq, endpoint = "continuous")

write_csv(results_cont_inv_1_alpha_eq, "results/results_cont_inv_1_alpha_eq.csv")
gc()


############################################################################################################

# Continuous case (equal time trends) - inverted-U trend, peak in the middle of period 1  - power


scenarios_cont_inv_1_pow_eq <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0.25,
                                          lambda0 = seq(-0.15, 0.15, length.out = 19),
                                          lambda1 = seq(-0.15, 0.15, length.out = 19),
                                          lambda2 = seq(-0.15, 0.15, length.out = 19),
                                          n0 = 250,
                                          n1 = 250,
                                          n01 = 250*0.5,
                                          n02 = 250-250*0.5,
                                          n11 = 250*0.5,
                                          n12 = 250-250*0.5,
                                          n22 = 250,
                                          trend = "inv_u",
                                          endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(26)
results_cont_inv_1_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_1_pow_eq, endpoint = "continuous")

write_csv(results_cont_inv_1_pow_eq, "results/results_cont_inv_1_pow_eq.csv")
gc()



############################################################################################################

# Continuous case (equal time trends) - inverted-U trend, peak between periods 1 and 2 - type one error


scenarios_cont_inv_2_alpha_eq <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0,
                                            lambda0 = seq(-0.15, 0.15, length.out = 19),
                                            lambda1 = seq(-0.15, 0.15, length.out = 19),
                                            lambda2 = seq(-0.15, 0.15, length.out = 19),
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(27)
results_cont_inv_2_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_2_alpha_eq, endpoint = "continuous")

write_csv(results_cont_inv_2_alpha_eq, "results/results_cont_inv_2_alpha_eq.csv")
gc()

############################################################################################################

# Continuous case (equal time trends) - inverted-U trend, peak between periods 1 and 2  - power


scenarios_cont_inv_2_pow_eq <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0.25,
                                          lambda0 = seq(-0.15, 0.15, length.out = 19),
                                          lambda1 = seq(-0.15, 0.15, length.out = 19),
                                          lambda2 = seq(-0.15, 0.15, length.out = 19),
                                          n0 = 250,
                                          n1 = 250,
                                          n01 = 250*0.5,
                                          n02 = 250-250*0.5,
                                          n11 = 250*0.5,
                                          n12 = 250-250*0.5,
                                          n22 = 250,
                                          trend = "inv_u",
                                          endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(28)
results_cont_inv_2_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_2_pow_eq, endpoint = "continuous")

write_csv(results_cont_inv_2_pow_eq, "results/results_cont_inv_2_pow_eq.csv")
gc()


############################################################################################################

# Continuous case (equal time trends) - inverted-U trend, peak in the middle of period 2 - type one error


scenarios_cont_inv_3_alpha_eq <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0,
                                            lambda0 = seq(-0.15, 0.15, length.out = 19),
                                            lambda1 = seq(-0.15, 0.15, length.out = 19),
                                            lambda2 = seq(-0.15, 0.15, length.out = 19),
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(29)
results_cont_inv_3_alpha_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_3_alpha_eq, endpoint = "continuous")

write_csv(results_cont_inv_3_alpha_eq, "results/results_cont_inv_3_alpha_eq.csv")
gc()

############################################################################################################

# Continuous case (equal time trends) - inverted-U trend, peak in the middle of period 2  - power


scenarios_cont_inv_3_pow_eq <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0.25,
                                          lambda0 = seq(-0.15, 0.15, length.out = 19),
                                          lambda1 = seq(-0.15, 0.15, length.out = 19),
                                          lambda2 = seq(-0.15, 0.15, length.out = 19),
                                          n0 = 250,
                                          n1 = 250,
                                          n01 = 250*0.5,
                                          n02 = 250-250*0.5,
                                          n11 = 250*0.5,
                                          n12 = 250-250*0.5,
                                          n22 = 250,
                                          trend = "inv_u",
                                          endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(30)
results_cont_inv_3_pow_eq <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_3_pow_eq, endpoint = "continuous")

write_csv(results_cont_inv_3_pow_eq, "results/results_cont_inv_3_pow_eq.csv")
gc()







############################################################################################################
# CONTINUOUS CASE - different time trends - inverted U-1 trend
############################################################################################################

# Continuous case (different time trends) - inverted U-1 trend - type one error


scenarios_cont_inv_1_alpha_diff <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0,
                                            lambda0 = 0.1,
                                            lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                            lambda2 = 0.1,
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)


set.seed(21)
results_cont_inv_1_alpha_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_1_alpha_diff, endpoint = "continuous")

write_csv(results_cont_inv_1_alpha_diff, "results/results_cont_inv_1_alpha_diff.csv")
gc()

############################################################################################################

# Continuous case (different time trends) - inverted U-1 trend - power

scenarios_cont_inv_1_pow_diff <- data.frame(K = 1,
                                          mu0 = 0,
                                          sigma = 1,
                                          delta1 = 0.25,
                                          delta2 = 0.25,
                                          lambda0 = 0.1,
                                          lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                          lambda2 = 0.1,
                                          n0 = 250,
                                          n1 = 250,
                                          n01 = 250*0.5,
                                          n02 = 250-250*0.5,
                                          n11 = 250*0.5,
                                          n12 = 250-250*0.5,
                                          n22 = 250,
                                          trend = "inv_u",
                                          endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(23)
results_cont_inv_1_pow_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_1_pow_diff, endpoint = "continuous")

write_csv(results_cont_inv_1_pow_diff, "results/results_cont_inv_1_pow_diff.csv")
gc()

############################################################################################################





############################################################################################################
# CONTINUOUS CASE - different time trends - inverted U-2 trend
############################################################################################################

# Continuous case (different time trends) - inverted U-2 trend - type one error


scenarios_cont_inv_2_alpha_diff <- data.frame(K = 1,
                                              mu0 = 0,
                                              sigma = 1,
                                              delta1 = 0.25,
                                              delta2 = 0,
                                              lambda0 = 0.1,
                                              lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                              lambda2 = 0.1,
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "inv_u",
                                              endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)


set.seed(21)
results_cont_inv_2_alpha_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_2_alpha_diff, endpoint = "continuous")

write_csv(results_cont_inv_2_alpha_diff, "results/results_cont_inv_2_alpha_diff.csv")
gc()

############################################################################################################

# Continuous case (different time trends) - inverted U-2 trend - power

scenarios_cont_inv_2_pow_diff <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0.25,
                                            lambda0 = 0.1,
                                            lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                            lambda2 = 0.1,
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11),
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(23)
results_cont_inv_2_pow_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_2_pow_diff, endpoint = "continuous")

write_csv(results_cont_inv_2_pow_diff, "results/results_cont_inv_2_pow_diff.csv")
gc()

############################################################################################################





############################################################################################################
# CONTINUOUS CASE - different time trends - inverted U-3 trend
############################################################################################################

# Continuous case (different time trends) - inverted U-3 trend - type one error


scenarios_cont_inv_3_alpha_diff <- data.frame(K = 1,
                                              mu0 = 0,
                                              sigma = 1,
                                              delta1 = 0.25,
                                              delta2 = 0,
                                              lambda0 = 0.1,
                                              lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                              lambda2 = 0.1,
                                              n0 = 250,
                                              n1 = 250,
                                              n01 = 250*0.5,
                                              n02 = 250-250*0.5,
                                              n11 = 250*0.5,
                                              n12 = 250-250*0.5,
                                              n22 = 250,
                                              trend = "inv_u",
                                              endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)


set.seed(21)
results_cont_inv_3_alpha_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_3_alpha_diff, endpoint = "continuous")

write_csv(results_cont_inv_3_alpha_diff, "results/results_cont_inv_3_alpha_diff.csv")
gc()

############################################################################################################

# Continuous case (different time trends) - inverted U-3 trend - power

scenarios_cont_inv_3_pow_diff <- data.frame(K = 1,
                                            mu0 = 0,
                                            sigma = 1,
                                            delta1 = 0.25,
                                            delta2 = 0.25,
                                            lambda0 = 0.1,
                                            lambda1 = round(seq(0.1-0.15, 0.1+0.15, by=0.01), 4),
                                            lambda2 = 0.1,
                                            n0 = 250,
                                            n1 = 250,
                                            n01 = 250*0.5,
                                            n02 = 250-250*0.5,
                                            n11 = 250*0.5,
                                            n12 = 250-250*0.5,
                                            n22 = 250,
                                            trend = "inv_u",
                                            endpoint = "continuous") %>%
  mutate(N_peak = (n01+n11)+(n02+n12+n22)/2,
         timetrend = ifelse(lambda0==lambda1 & lambda1==lambda2, "EQ", "DIFF"),
         hypothesis = ifelse(delta2==0, "H0", "H1"),
         timing = n01/n0)



set.seed(23)
results_cont_inv_3_pow_diff <- allinone_simsce_par(nsim=nsim, scenarios=scenarios_cont_inv_3_pow_diff, endpoint = "continuous")

write_csv(results_cont_inv_3_pow_diff, "results/results_cont_inv_3_pow_diff.csv")
gc()

############################################################################################################