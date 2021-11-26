Model-based time trend adjustments in platform trials with
non-concurrent controls
================

This repository contains the main code and R functions to reproduce the
results of the simulation study in \`\`On model-based time trend
adjustments in platform trials with non-concurrent controls’’ (XXX).

We investigated frequentist, model-based approaches to adjust for time
trends in platform trials utilizing non-concurrent controls. We
investigated conditions under which the model-based approaches can
successfully adjust for time trends for the simple case of a two-period
trial with two experimental treatments and a shared control. We focus on
a platform trial with two treatment arms and a common control arm.
Assuming that the second treatment arm is added at a later time period,
we assess the robustness of recently proposed model-based approaches to
adjust for time trends when utilizing non-concurrent controls. The
following figure illustrates the design considered:

![Image Title](./Figs/trial_scheme.png)

Outline structure repository.

## Main simulation script and results:

  - **mainpaper\_script.R**: Script with all scenarios considered in the
    study. Simulates given scenarios and saves the results as .csv files
    in the folder *results*.

Running this script results in 100 separate .csv files, named
**results\_bin\_trendpattern\_parametrization\_hypothesis\_trendtype.csv**
or **results\_cont\_trendpattern\_hypothesis\_trendtype.csv** depending
on whether the simulated data were binary or continuous.

  - For binary endpoints, the given scenario is indicated in the file
    name as
    **results\_bin\_trendpattern\_parametrization\_hypothesis\_trendtype.csv**,
    where
    
      - for **trendpattern**, we considered the options:
        
          - **inv\_1** - inverted-U trend with peak in the middle of
            period 1
          - **inv\_2** - inverted-U trend with peak between period 1 and
            2
          - **inv\_3** - inverted-U trend with peak in the middle of
            period 2
          - **lin** - linear trend
          - **step** - stepwise trend
    
      - for **parametrization**, we considered:
        
          - **add** - additive parametrization
          - **mult** - multiplicative parametrization
    
      - for **hypothesis**, we considered:
        
          - **alpha** - simulations of the type I error under the null
            hypothesis
          - **pow** - simulations of the power under the alternative
            hypothesis
    
      - for **trendtype**, we considered:
        
          - **eq** - equal time trend in all arms
          - **diff\_pos** - positive time trend in control and arm 2,
            varying trend in arm 1
          - **diff\_neg** - negative time trend in control and arm 2,
            varying trend in arm 1
    
      - additionally, **OR1** refers to scenarios with OR1\>1, otherwise
        OR1\<1 is used

  - For continuous endpoints, the given scenario is indicated in the
    file name as follows:
    **results\_cont\_trendpattern\_hypothesis\_trendtype.csv**, where
    
      - for **trendpattern**, we considered:
        
          - **inv\_1** - inverted-U trend with peak in the middle of
            period 1
          - **inv\_2** - inverted-U trend with peak between period 1 and
            2
          - **inv\_3** - inverted-U trend with peak in the middle of
            period 2
          - **lin** - linear trend
          - **step** - stepwise trend
    
      - for **hypothesis**, we considered:
        
          - **alpha** - simulations of the type I error under the null
            hypothesis
          - **pow** - simulations of the power under the alternative
            hypothesis
    
      - for **trendtype**, we considered:
        
          - **eq** - equal time trend in all arms
          - **diff** - time trend of lambda=0.1 in control and arm 2,
            varying trend in arm 1

## R Functions:

  - **trend\_functions.R**: Contains functions `linear_trend()` and
    `sw_trend()` that are used to generate linear or stepwise trend.
  - **data\_sim.R**: The function `data_sim()` is used to simulate
    continuous or binary data for two-stage platform trial, based on
    given sample sizes, effect sizes and time trend specifications.
  - **z\_test.R**: Contains functions `z_test_pol()` and `z_test_sep()`
    that perform pooled or separate one-sided z-test with the given
    data.
  - **t\_test.R**: Contains functions `t_test_pol()` and `t_test_sep()`
    that perform pooled or separate one-sided t-test with the given
    data.
  - **linear\_model.R**: Contains functions `linear_model_a1()`,
    `linear_model_a1_int()`, `linear_model_a2()`,
    `linear_model_a2_int()`, `linear_model_b1()` and `linear_model_b2()`
    that are used to fit different linear regression models
    (ALLTC-Linear, ALLTCI-Linear, ALLTC-Step, ALLTCI-Step, TC-Linear and
    TC-Step) to given data.
  - **z\_prop\_test.R**: Contains function `z_prop_pol()` and
    `z_prop_sep()` that perform pooled or separate one-sided
    two-proportions z-test for with the given data.
  - **log\_model.R**: Contains functions `logistic_model_a1()`,
    `logistic_model_a1_int()`, `logistic_model_a2()`,
    `logistic_model_a2_int()`, `logistic_model_b1()`,
    `logistic_model_b2()`, `logistic_model_sep()` and
    `logistic_model_pol()` that are used to fit different logistic
    regression models (ALLTC-Linear, ALLTCI-Linear, ALLTC-Step,
    ALLTCI-Step, TC-Linear, TC-Step and separate and pooled logistic
    model) to given data.
  - **allinone\_model.R**: The function `allinone_model()` is programmed
    to call all the aforementioned functions and return the calculated
    metrics (T/F for reject \(H_{02}\) and bias) from all considered
    tests and models.
  - **allinone\_sim.R**: The function `allinone_simsce()` takes a data
    frame with different simulation scenarios as input and performs
    *nsim* replications for each of them to simulate the type I error
    rate/power, bias and rMSE.
  - **allinone\_sim\_par.R**: The function `allinone_simsce_par()` is a
    parallel version of `allinone_simsce()`

## Plots:

  - **Plots\_paper.Rmd**: generates all figures used in the paper and
    saves them in the folder *figures*.

<!-- ## Example?? -->

<!-- This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated. -->

<!-- ... Including Plots -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo=FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->
