Model-based time trend adjustments in platform trials with
non-concurrent controls
================

This repository contains the main code and R functions to reproduce the
results of the simulation study in **“On model-based time trend
adjustments in platform trials with non-concurrent controls”**
\[[Paper](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-022-01683-w)\]
\[[Preprint](https://arxiv.org/abs/2112.06574)\].

We investigated frequentist, model-based approaches to adjust for time
trends in platform trials utilizing non-concurrent controls. We
investigated conditions under which the model-based approaches can
successfully adjust for time trends in the simple case of a two-period
trial with two experimental treatments and a shared control. Assuming
that the second treatment arm is added at a later time period, we assess
the robustness of recently proposed model-based approaches to adjust for
time trends when utilizing non-concurrent controls. The following figure
illustrates the design considered:

<img src="./figures/trial_scheme.png" style="width:30.0%" />

Next, we explain the R functions we created to simulate the data and to
analyse it. After this, we describe the main simulation script and
results and briefly summarize the supplementary material and the figures
created.

# R Functions

The folder *functions* contains all R functions necessary to reproduce
the results and figures in the main paper and supplementary material:

## Functions for data generation

-   **trend_functions.R**: Contains functions `linear_trend()` and
    `sw_trend()` that are used to generate linear or stepwise trend in
    the main paper and the function `linear_trend2()` that is used to
    simulate linear trend that starts in the second period in the
    supplementary material.

-   **data_sim_block.R**: The function `data_sim_block()` simulates
    continuous or binary data for two-stage platform trial using block
    randomization and deterministic patient entry times, based on given
    sample sizes, effect sizes and time trend specifications. This
    function is used for data generation in the main paper.

-   **data_sim_brdt.R**: The function `data_sim_brdt()` simulates
    continuous or binary data for two-stage platform trial using block
    randomization and random patient entry times.

-   **data_sim.R**: The function `data_sim()` simulates continuous or
    binary data for two-stage platform trial using simple randomization
    per period and deterministic patient entry times.

-   **data_sim_rdt.R**: The function `data_sim_rdt()` simulates
    continuous or binary data for two-stage platform trial using simple
    randomization per period and random patient entry times.

## Functions for tests and models

-   **z_test.R**: Contains functions `z_test_pol()` and `z_test_sep()`
    that perform pooled or separate one-sided z-test with the given
    data.

-   **t_test.R**: Contains functions `t_test_pol()` and `t_test_sep()`
    that perform pooled or separate one-sided t-test with the given
    data.

-   **linear_model.R**: Contains functions `linear_model_a1()`,
    `linear_model_a1_int()`, `linear_model_a2()`,
    `linear_model_a2_int()`, `linear_model_b1()` and `linear_model_b2()`
    that are used to fit different linear regression models
    (ALLTC-Linear, ALLTCI-Linear, ALLTC-Step, ALLTCI-Step, TC-Linear and
    TC-Step) to given data.

-   **z_prop_test.R**: Contains function `z_prop_pol()` and
    `z_prop_sep()` that perform pooled or separate one-sided
    two-proportions z-test for with the given data.

-   **log_model.R**: Contains functions `logistic_model_a1()`,
    `logistic_model_a1_int()`, `logistic_model_a2()`,
    `logistic_model_a2_int()`, `logistic_model_b1()`,
    `logistic_model_b2()`, `logistic_model_sep()` and
    `logistic_model_pol()` that are used to fit different logistic
    regression models (ALLTC-Linear, ALLTCI-Linear, ALLTC-Step,
    ALLTCI-Step, TC-Linear, TC-Step and separate and pooled logistic
    model) to given data.

## Functions for simulations

-   **allinone_model.R**: The function `allinone_model()` is programmed
    to call all the aforementioned functions and return the calculated
    metrics (T/F for reject H02 and bias) from all considered tests and
    models.

-   **allinone_sim.R**: The function `allinone_simsce()` takes a data
    frame with different simulation scenarios as input and performs
    *nsim* replications for each of them to simulate the type I error
    rate/power, bias and rMSE.

-   **allinone_sim_par.R**: The function `allinone_simsce_par()` is a
    parallelized version of `allinone_simsce()`.

# Main paper

## Main simulation script and results

-   **mainpaper_script.R**: Script with all scenarios considered in the
    main study. Simulates given scenarios and saves the results as .csv
    files in the folder *results*.

Running this script results in 100 separate .csv files, named
**results_bin_trendpattern_parametrization_hypothesis_trendtype.csv** or
**results_cont_trendpattern_hypothesis_trendtype.csv** depending on
whether the simulated data were binary or continuous.

-   For binary endpoints, the given scenario is indicated in the file
    name as
    **results_bin_trendpattern_parametrization_hypothesis_trendtype.csv**,
    where

    -   for **trendpattern**, we considered the options:

        -   **inv_1** - inverted-U trend with peak in the middle of
            period 1
        -   **inv_2** - inverted-U trend with peak between period 1 and
            2
        -   **inv_3** - inverted-U trend with peak in the middle of
            period 2
        -   **lin** - linear trend
        -   **step** - stepwise trend

    -   for **parametrization**, we considered:

        -   **add** - additive parametrization
        -   **mult** - multiplicative parametrization

    -   for **hypothesis**, we considered:

        -   **alpha** - simulations of the type I error under the null
            hypothesis
        -   **pow** - simulations of the power under the alternative
            hypothesis

    -   for **trendtype**, we considered:

        -   **eq** - equal time trend in all arms
        -   **diff_pos** - positive time trend in control and arm 2,
            varying trend in arm 1
        -   **diff_neg** - negative time trend in control and arm 2,
            varying trend in arm 1

    -   additionally, **OR1** refers to scenarios with OR1\>1, otherwise
        OR1\<1 is used

-   For continuous endpoints, the given scenario is indicated in the
    file name as follows:
    **results_cont_trendpattern_hypothesis_trendtype.csv**, where

    -   for **trendpattern**, we considered:

        -   **inv_1** - inverted-U trend with peak in the middle of
            period 1
        -   **inv_2** - inverted-U trend with peak between period 1 and
            2
        -   **inv_3** - inverted-U trend with peak in the middle of
            period 2
        -   **lin** - linear trend
        -   **step** - stepwise trend

    -   for **hypothesis**, we considered:

        -   **alpha** - simulations of the type I error under the null
            hypothesis
        -   **pow** - simulations of the power under the alternative
            hypothesis

    -   for **trendtype**, we considered:

        -   **eq** - equal time trend in all arms
        -   **diff** - time trend of lambda=0.1 in control and arm 2,
            varying trend in arm 1

# Figures

-   **Plots_paper.Rmd**: generates all figures used in the main paper
    and section E of the supplementary material and saves them in the
    folder *figures*.

Figures presented in the main paper have the suffix **\_main**
(e.g. **cont_response_main.png**), while figures included in the
supplementary material are indicated by the suffix **\_supp**
(e.g. **cont_all_eq_bias_supp.png**).

# Supplementary material

The folder *Supp_mat* contains additional simulation results for
three-stage platform trials and additional results regarding different
randomization procedures discussed in the supplementary material.

Results that are reported in section F.1 of the supplementary material
are simulated in the **Supp_mat_randomization.Rmd** script and presented
in the corresponding HTML file.

For section F.2 of the supplementary material (platform trials with
three periods), there are analogous functions as for the main paper in
the subfolder *functions*, which were adapted for the three-stage
design. All considered three-stage scenarios can be simulated by running
the **additional_sim_script.R**, while the results are again saved in
the *results* subfolder. The files **Plots_additional_sim.Rmd** and the
corresponding HTML file are used for visualization of the additional
results. Figures presented in Section F.2 of the Supplementary material
are saved in the *figures* subfolder.

------------------------------------------------------------------------

**Funding**

[EU-PEARL](https://eu-pearl.eu/) (EU Patient-cEntric clinicAl tRial
pLatforms) project has received funding from the Innovative Medicines
Initiative (IMI) 2 Joint Undertaking (JU) under grant agreement No
853966. This Joint Undertaking receives support from the European
Union’s Horizon 2020 research and innovation programme and EFPIA
andChildren’s Tumor Foundation, Global Alliance for TB Drug Development
non-profit organisation, Spring works Therapeutics Inc. This publication
reflects the authors’ views. Neither IMI nor the European Union, EFPIA,
or any Associated Partners are responsible for any use that may be made
of the information contained herein.

<!-- ## Example?? -->
<!-- This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated. -->
<!-- ... Including Plots -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo=FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot. -->
