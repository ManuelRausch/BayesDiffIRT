<!-- Dieses RMarkdown Dokument dient dazu, ein Readme-Markdown-Dokument zu knitten.  Dieses Read-me wird man sehen, wenn man auf die Paket Development Page of Github geht.  -->

# 1 The BayesDiffIRT for R

The `BayesDiffIRT` package provides functions to sample posterior
distributions and posterior predictive distributions of item and subject
parameters of diffusion item response theory models for responses and
reaction times Kang, De Boeck, and Ratcliff (2022). `BayesDiffIRT` also
provides functions to visualize posterior distributions of Diffusion
item response theory model parameters and construct credible intervals.
Under the hood, the package relies on NUTS sampling with STAN (Carpenter
et al. 2017). Up to know, the following diffusion item response theory
models have been implemented:

- D-diffusion model (Tuerlinckx and Boeck 2005),
- Q-diffusion model (Van Der Maas et al. 2011),
- D-diffusion model with random variability (Kang, De Boeck, and
  Ratcliff 2022),
- Q-diffusion model with random variability (Kang, De Boeck, and
  Ratcliff 2022).

The two versions of the D-diffusion model are appropriate for survey
items where persons decide whether to accept or reject an item. The two
flavours of the Q-diffusion model were designed to model ability tests.

*Important*: BayesDiffIRT uses Stan through the cmdstanr interface.
Fitting models therefore requires a C++ toolchain and a separate
installation of CmdStan. Installing the BayesDiffIRT R package alone is
not sufficient.

# 2 Mathematical description of diffusion item response theory models

Diffusion item response theory combines item response theory with the
drift diffusion model of decision making. According to the drift
diffusion decision model (Stone 1960; Link and Heath 1975; Ratcliff et
al. 2016), the sensory system repeatedly generates momentary evidence
about which of two possible choice options is correct. This momentary
evidence is drawn from a Gaussian distribution and accumulated over
time. The newly acquired evidence is therefore continuously added to the
evidence collected up to that moment.The accumulation process is bounded
by an upper and a lower threshold, where each threshold represents one
of the two possible choice options. When the accumulated evidence
reaches one of the thresholds, a choice is made for the corresponding
option. The quality of information favouring one response option over
the other is reflected in the drift rate $`\delta`$, which quantifies
how quickly the accumulated evidence approaches the threshold associated
with the correct or preferred decision. The distance between the two
thresholds, $`\alpha`$, determines the amount of evidence required
before a decision is made; a larger distance means that decisions tend
to be slower because more evidence is required. The starting point
$`\beta`$ of the accumulation process reflects an a priori bias toward
one of the response options. To visualize the reaction time
distributions that follow from the drift diffusion model, interactive
tools such as the [diffusion model visualizer](https://osf.io/4en3b) are
openly available (Alexandrowicz 2020).

In diffusion item response theory models, two traditional parameters of
the drift diffusion model, boundary separation and drift rate, are
decomposed into person and item parameters (Tuerlinckx and Boeck 2005).
When person $`p`$ makes a decision about item $`i`$, the boundary
separation is given by

``` math

\alpha_{pi} = \frac{\gamma_p}{a_i},
```

where $`\gamma_p`$ represents person-specific response caution and
$`a_i`$ represents item-specific time pressure.

The D-diffusion and Q-diffusion models differ in how the drift rate is
decomposed. In the D-diffusion model (Tuerlinckx and Boeck 2005), which
is applicable to survey items, the drift rate is given by

``` math

\delta_{pi} = \theta_p - \nu_i.
```

According to the Q-diffusion model (Van Der Maas et al. 2011), which is
applicable to ability tests, the drift rate is given by

``` math

\delta_{pi} = \frac{\theta_p}{\nu_i}.
```

In both the D-diffusion and Q-diffusion models, the accumulation process
starts midway between the two response alternatives. Thus, there is
assumed to be no a priori bias toward either choice alternative.

Kang, De Boeck, and Ratcliff (2022) proposed extensions that include
random trial-to-trial variability in both the starting point $`\beta`$
and the drift rate $`\delta`$. In the Q- and D-diffusion models with
random variation, the relative starting point $`\beta_{pij}`$ for trial
$`j`$, item $`i`$, and person $`p`$ is sampled from a uniform
distribution,

``` math

\beta_{pij} \sim \mathcal{U}\left(0.5 - \frac{s_{\beta}}{2}, 0.5 + \frac{s_{\beta}}{2}\right).
```

The drift rate $`\delta_{pij}`$ for trial $`j`$, item $`i`$, and person
$`p`$ is sampled from a Gaussian distribution,

``` math

\delta_{pij} \sim \mathcal{N}\left(\delta_{pi}, s_{\delta}^2\right).
```

Random variability in starting points and drift rates accounts for the
conditional dependency between accuracy and reaction times (Kang, De
Boeck, and Ratcliff 2022), but note that sampling is considerably slower
for these models.

# 3 Installation

## 3.1 System requirements

The `BayesDiffIRT` package relies on STAN via cmdstanr. Stan models are
translated into C++ and compiled before they are run. Consequently, a
working C++ compiler and GNU Make are required. Please install the
appropriate toolchain for your operating system before installing
CmdStan:

- Windows: Install the version of
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) appropriate
  for your version of R.
- macOS: Install the Xcode Command Line Tools by entering the following
  command in the Terminal:

``` bash
xcode-select --install
```

- Debian/Ubuntu Linux: Install g++ and make:

``` bash
sudo apt update
sudo apt install g++ make
```

Further operating-system-specific information is available in the
[CmdStan installation
guide](https://mc-stan.org/docs/cmdstan-guide/installation.html).

## 3.2 Installing CmdStan

If `cmdstanr` is not yet installed on oyur system, install it from the
Stan R-universe repository:

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

First check whether the C++ toolchain is correctly configured:

``` r
cmdstanr::check_cmdstan_toolchain()
```

Once the toolchain check succeeds, install CmdStan:

``` r
cmdstanr::install_cmdstan()
```

Verify the installation with:

``` r
cmdstanr::cmdstan_version()
```

CmdStan is installed separately from the R package, normally in the
.cmdstan directory in the user’s home directory. The initial compilation
of CmdStan may take several minutes and require substantial memory.

## 3.3 Installing BayesDiffIRT

The development version of `BayesDiffIRT` is available from GitHub. The
easiest way to install it is using the `devtools` package:

<!-- without any dots, the code chunk will be shown, but not executed -->

    devtools::install_github("ManuelRausch/BayesDiffIRT")

# 4 Usage

## 4.1 Fitting drift diffusion item response theory models

The function `fitBayesDiffIRT` fits Bayesian diffusion item-response
theory models by sampling from the posterior distributions of item and
subject parameters using the No-U-Turn Sampler (NUTS) as implemented in
Stan (Carpenter et al. 2017). The data should be a dataframe with
columns identifying the subject, item, response, and response time.
Response times should be numeric and measured in seconds. For ability
tests, binary responses should be coded as 0 for incorrect responses and
1 for correct responses. For questionnaire items, binary responses
should be coded as 0 for rejected items and 1 for accepted items. Here,
we prepare a dataset contained in the diffIRT package as example.

``` r
library(tidyverse)
```

``` r
# Example for preparing the data set.

data(extraversion, package = "diffIRT")
Extra <- as.data.frame(extraversion)
names(Extra)[1:10]  <- paste0("Item", 1:10, "_resp")
names(Extra)[11:20] <- paste0("Item", 1:10, "_rt")
Extra$sbj <- 1:nrow(Extra)
  
Extra <- tidyr::pivot_longer(
  Extra,
  cols = tidyselect::matches("^Item\\d+_(resp|rt)$"),
  names_to = c("item", ".value"),
  names_pattern = "^(Item\\d+)_(resp|rt)$"
)

Extra$item <- factor(Extra$item)
Extra$item <- factor(Extra$item)
head(Extra)
```

    ## # A tibble: 6 × 4
    ##     sbj item   resp    rt
    ##   <int> <fct> <dbl> <dbl>
    ## 1     1 Item1     0 2.73 
    ## 2     1 Item2     1 0.915
    ## 3     1 Item3     1 3.48 
    ## 4     1 Item4     1 1.02 
    ## 5     1 Item5     1 1.11 
    ## 6     1 Item6     0 2.38

Priors distributions for parameter classes are specified using the
function `prior()`. The function creates objects of class
`BayesDiffIRTPrior`, which can be passed to the `priors` argument of
`fitBayesDiffIRT`. Priors are specified using Stan distribution syntax,
for example `normal(0, 1)` or `lognormal(0, 0.5)`. Priors on the
following parameters can be specified:

| Parameter | Description |
|----|----|
| `"omega_theta"` | Standard deviation of `theta`, the latent trait. |
| `"omega_gamma"` | Standard deviation of `gamma`, the item response tendency parameter. |
| `"nu"` | Item difficulty parameter. |
| `"a"` | Item boundary separation parameter. |
| `"tnd"` | Person-specific non-decision time. |
| `"s_delta"` | Standard deviation of Gaussian trial-to-trial variability in drift rate. Only relevant for models `"dRV"` and `"qRV"`. |
| `"s_beta"` | Range of uniform trial-to-trial variability in starting point. Only relevant for models `"dRV"` and `"qRV"`. |

If priors are supplied for only some parameter classes, the remaining
parameter classes are filled in with their model-specific defaults.

``` r
# Example for prior specification
library(BayesDiffIRT)
```

    ## BayesDiffIRT: Bayesian diffusion-IRT models with RTs.
    ## Backend: cmdstanr.

``` r
myPrior <-  list(prior(normal(0, 2.5), class = "omega_theta"),
               prior(normal(0, 0.5), class = "omega_gamma"),
               prior(lognormal(0, 0.75), class = "nu"),
               prior(lognormal(0, 0.5), class = "a"),
               prior(lognormal(-1.25, 0.3), class = "tnd"))
```

Finally, a string identifying the chosen model, the prepared data and a
list of priors should be passed to `fitBayesDiffIRT` to fit a Bayesian
diffusion item-response theory model. The following model names are
recognized:

- “d” for the D-diffusion model (for survey items, default),
- “dRV” for the D-diffusion model with random variability (for survey
  items),
- “q” for the Q-diffusion model (for ability tests),
- “qRV” for the Q-diffusion model with random variability (for ability
  tests).

``` r
samples <- 
  fitBayesDiffIRT(Extra,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "d")
```

    ## Running MCMC with 4 parallel chains...

    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 4 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 2 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 1 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 3 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 4 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 2 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 1 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 3 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 4 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 2 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 1 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 3 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 4 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 2 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 1 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 2 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 4 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 1 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 3 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 2 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 4 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 1 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 3 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 2 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 4 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 1 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 3 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 2 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 4 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 1 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 192.3 seconds.
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 193.3 seconds.
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 194.6 seconds.
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 198.6 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 194.7 seconds.
    ## Total execution time: 199.0 seconds.

## 4.2 Inspecting the results

The results of a fitted Bayesian diffusion item-response theory model
can be inspected using the `summary` method.

``` r
summary(samples)
```

    ## Summary of BayesDiffIRT model fit
    ## ---------------------------------
    ## Model: D-Diffusion model 
    ## 
    ## Data:
    ##   Observations: 1429 
    ##   Persons:      143 
    ##   Items:        10 
    ## 
    ## Call:
    ## fitBayesDiffIRT(data = Extra, rt = "rt", resp = "resp", sbj = "sbj", 
    ##     item = "item", model = "d")
    ## 
    ## Posterior summaries:
    ## 
    ## Hyperparameters:
    ## # A tibble: 2 × 9
    ##   variable     mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##   <chr>       <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ## 1 omega_theta  0.69   0.69  0.06  0.6   0.8      1    1196.    2002.
    ## 2 omega_gamma  0.2    0.2   0.03  0.16  0.25     1     927.    1560.
    ## 
    ## Item parameters:
    ## # A tibble: 20 × 9
    ##    variable  mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 nu[1]    -0.64  -0.64  0.11 -0.82 -0.47     1    1883.    2783.
    ##  2 nu[2]    -0.14  -0.14  0.11 -0.32  0.03     1    1928.    2852.
    ##  3 nu[3]    -1.23  -1.23  0.13 -1.44 -1.02     1    2038.    3015.
    ##  4 nu[4]    -1.7   -1.7   0.15 -1.94 -1.46     1    2378.    2650.
    ##  5 nu[5]    -0.21  -0.2   0.11 -0.38 -0.03     1    2011.    2663.
    ##  6 nu[6]    -1.3   -1.3   0.12 -1.5  -1.1      1    1991.    3005.
    ##  7 nu[7]    -1.69  -1.69  0.14 -1.93 -1.45     1    2621.    2885.
    ##  8 nu[8]    -1.91  -1.91  0.15 -2.15 -1.67     1    2472     2949.
    ##  9 nu[9]    -0.83  -0.83  0.1  -1    -0.66     1    1994.    2739.
    ## 10 nu[10]   -1.42  -1.42  0.14 -1.65 -1.2      1    2075.    2687.
    ## 11 a[1]      0.44   0.44  0.02  0.41  0.47     1    2248.    2846.
    ## 12 a[2]      0.5    0.49  0.02  0.46  0.53     1    2417.    2634.
    ## 13 a[3]      0.5    0.5   0.02  0.46  0.54     1    2194.    2912.
    ## 14 a[4]      0.51   0.51  0.03  0.47  0.55     1    1859.    2741.
    ## 15 a[5]      0.51   0.51  0.02  0.48  0.55     1    2455.    2364.
    ## 16 a[6]      0.43   0.43  0.02  0.4   0.47     1    2609.    2939.
    ## 17 a[7]      0.4    0.4   0.02  0.37  0.44     1    2379.    2787.
    ## 18 a[8]      0.42   0.42  0.02  0.39  0.46     1    2061.    2944.
    ## 19 a[9]      0.35   0.35  0.02  0.32  0.38     1    3313.    2720.
    ## 20 a[10]     0.55   0.55  0.03  0.51  0.6      1    2019.    2917.
    ## 
    ## Subject parameters:
    ## # A tibble: 429 × 9
    ##    variable  mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 tnd[1]    0.37   0.36  0.1   0.21  0.54     1    5128.    3013.
    ##  2 tnd[2]    0.4    0.41  0.09  0.25  0.54     1    4391.    3343.
    ##  3 tnd[3]    0.46   0.46  0.12  0.27  0.65     1    4489.    3092.
    ##  4 tnd[4]    0.27   0.27  0.06  0.18  0.36     1    6497.    3405.
    ##  5 tnd[5]    0.39   0.39  0.07  0.26  0.51     1    4797.    3017.
    ##  6 tnd[6]    0.32   0.32  0.06  0.21  0.41     1    4875.    3518.
    ##  7 tnd[7]    0.39   0.39  0.09  0.24  0.53     1    5252.    2968.
    ##  8 tnd[8]    0.52   0.52  0.15  0.27  0.77     1    4475.    3298.
    ##  9 tnd[9]    0.34   0.34  0.07  0.22  0.45     1    6166.    3389.
    ## 10 tnd[10]   0.43   0.42  0.13  0.23  0.66     1    4652.    2719.
    ## # ℹ 419 more rows
    ## 
    ## Other parameters:
    ## # A tibble: 429 × 9
    ##    variable     mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>       <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 z_theta[1]  -0.82  -0.82  0.36 -1.42 -0.23     1    4618.    2650.
    ##  2 z_theta[2]  -0.3   -0.31  0.56 -1.2   0.63     1    6430.    2810.
    ##  3 z_theta[3]   0.4    0.4   0.53 -0.45  1.26     1    5697.    2800.
    ##  4 z_theta[4]   0.49   0.48  0.46 -0.24  1.26     1    4973.    3109.
    ##  5 z_theta[5]   1.46   1.45  0.67  0.39  2.57     1    5935.    3102.
    ##  6 z_theta[6]  -0.07  -0.07  0.52 -0.91  0.79     1    5366.    2855.
    ##  7 z_theta[7]  -0.07  -0.08  0.5  -0.88  0.77     1    7256     3147.
    ##  8 z_theta[8]   0.33   0.31  0.54 -0.54  1.24     1    4945.    2614.
    ##  9 z_theta[9]   0.26   0.26  0.47 -0.5   1.04     1    5464.    3222.
    ## 10 z_theta[10]  0.38   0.37  0.44 -0.34  1.12     1    4243.    3271.
    ## # ℹ 419 more rows

The method `checkDiagnostics` provides common Stan diagnostics such as
number of divergeant transitions, R hat, effective sample size, maximum
treedepth hits, and E-BFMI.

``` r
checkDiagnostics(samples)
```

    ## Stan diagnostics
    ## ----------------
    ## Chains:               4
    ## Divergences:          0
    ## Max treedepth hits:   0
    ## R-hat warnings:       0
    ## Low ESS warnings:     0
    ## E-BFMI warnings:      0
    ## 
    ## Overall status:        OK

Marcov chains of selected parameters can be visualized using the `plot`
method:

``` r
plot(samples, parameter = "omega_theta", type = "trace")
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

There is also the possibility to plot posterior means with 50% and 95%
credible intervals as well as marginal posterior densities.

``` r
plot(samples, parameter = "theta", type = "interval")
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
plot(samples, parameter = "omega_theta",
     type = "density")
```

![](README_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

The function `fitBayesDiffIRT` returns a `BayesDiffIRTfit`-object. The
samples from a `BayesDiffIRTfit` can be extracted using the
extractSamples method.

``` r
samplesDf <- extractSamples(samples)
head(samplesDf)
```

    ## # A draws_df: 6 iterations, 1 chains, and 881 variables
    ##    lp__ z_theta[1] z_theta[2] z_theta[3] z_theta[4] z_theta[5] z_theta[6]
    ## 1 -1730      -0.38       0.46      0.537       1.02       1.83     -0.351
    ## 2 -1753      -0.87       0.20     -0.016       0.26       1.14      0.964
    ## 3 -1753      -1.33       0.65      0.306      -0.54       1.81     -0.659
    ## 4 -1731      -0.91      -0.71      0.121       1.01       0.97     -0.074
    ## 5 -1739      -0.93       0.20      0.690       0.53       1.94      0.052
    ## 6 -1761      -0.70      -0.14     -0.086       1.17       0.93     -0.263
    ##   z_theta[7]
    ## 1      0.713
    ## 2     -0.384
    ## 3     -0.004
    ## 4     -0.297
    ## 5     -0.016
    ## 6     -0.210
    ## # ... with 873 more variables
    ## # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

Should you really prefer working with point estimates, you can extract
them conveniently with the coef method:

``` r
pointEstim <- coef(samples, parameter = "theta")
pointEstim 
```

    ##     theta[1]     theta[2]     theta[3]     theta[4]     theta[5]     theta[6] 
    ## -0.566553822 -0.207144131  0.276347686  0.339556607  1.016869532 -0.051757098 
    ##     theta[7]     theta[8]     theta[9]    theta[10]    theta[11]    theta[12] 
    ## -0.049124830  0.225043448  0.177682710  0.262116813 -0.037225211  0.357596628 
    ##    theta[13]    theta[14]    theta[15]    theta[16]    theta[17]    theta[18] 
    ##  0.297683098 -0.038802429  0.506702750  0.323905827 -0.107285211 -0.052698875 
    ##    theta[19]    theta[20]    theta[21]    theta[22]    theta[23]    theta[24] 
    ## -0.251394224  0.674551131 -1.109402141  1.209010279  0.461525496  0.076621953 
    ##    theta[25]    theta[26]    theta[27]    theta[28]    theta[29]    theta[30] 
    ##  0.137701740 -0.920250066  1.125526110  0.236256192  0.522456391 -0.781602656 
    ##    theta[31]    theta[32]    theta[33]    theta[34]    theta[35]    theta[36] 
    ##  0.427789631 -0.863977841  0.717027610  0.048131979  0.081382202  0.154339075 
    ##    theta[37]    theta[38]    theta[39]    theta[40]    theta[41]    theta[42] 
    ## -0.133913959  0.274542779 -0.023923453  0.335863847 -1.082581325 -0.365929664 
    ##    theta[43]    theta[44]    theta[45]    theta[46]    theta[47]    theta[48] 
    ##  0.732519269  0.405392619  0.632389409 -0.025661872  0.033230806 -0.284628068 
    ##    theta[49]    theta[50]    theta[51]    theta[52]    theta[53]    theta[54] 
    ## -0.623016039  0.412333983  0.412316196 -0.210976060  1.176806917 -0.044019832 
    ##    theta[55]    theta[56]    theta[57]    theta[58]    theta[59]    theta[60] 
    ##  0.006202363  1.165195454  0.423766648  0.434375253 -1.157838265 -1.047370140 
    ##    theta[61]    theta[62]    theta[63]    theta[64]    theta[65]    theta[66] 
    ## -1.093013230 -0.075531164 -0.476655882  1.114525314 -0.125722589 -0.694534778 
    ##    theta[67]    theta[68]    theta[69]    theta[70]    theta[71]    theta[72] 
    ## -0.709283970  0.175478050  0.904396542  0.800220586 -0.393383129  0.494922719 
    ##    theta[73]    theta[74]    theta[75]    theta[76]    theta[77]    theta[78] 
    ## -0.391857260  0.995022431 -0.132605809 -1.681316677  0.265888587 -0.680846195 
    ##    theta[79]    theta[80]    theta[81]    theta[82]    theta[83]    theta[84] 
    ## -0.260345029 -0.187531201 -0.304048556 -0.509736767  0.407348854  0.126802226 
    ##    theta[85]    theta[86]    theta[87]    theta[88]    theta[89]    theta[90] 
    ##  0.360823200 -0.083953890  0.050689409 -0.394370843 -0.196775519  0.287670583 
    ##    theta[91]    theta[92]    theta[93]    theta[94]    theta[95]    theta[96] 
    ## -0.016298083 -1.355877315  0.320284680 -0.373635254 -0.212382886 -0.491910410 
    ##    theta[97]    theta[98]    theta[99]   theta[100]   theta[101]   theta[102] 
    ## -0.137917264 -0.342538949 -0.936902841 -0.643806293  0.136593266  0.910550152 
    ##   theta[103]   theta[104]   theta[105]   theta[106]   theta[107]   theta[108] 
    ## -0.200221599 -0.320735751  0.010879582  0.092861743  0.369232465 -0.147104553 
    ##   theta[109]   theta[110]   theta[111]   theta[112]   theta[113]   theta[114] 
    ##  0.236507181 -0.786256278 -0.190663344  0.541717764 -0.577828571  1.065456751 
    ##   theta[115]   theta[116]   theta[117]   theta[118]   theta[119]   theta[120] 
    ##  0.050503466  0.048421094 -0.383720220 -0.414369975  0.029222907  0.168790919 
    ##   theta[121]   theta[122]   theta[123]   theta[124]   theta[125]   theta[126] 
    ## -1.011600940  0.808630826  0.517489160 -1.383111954  0.084382017  0.778863460 
    ##   theta[127]   theta[128]   theta[129]   theta[130]   theta[131]   theta[132] 
    ##  0.452147943  0.059403985  0.016359253  0.625891467  0.159281637  0.482535762 
    ##   theta[133]   theta[134]   theta[135]   theta[136]   theta[137]   theta[138] 
    ##  0.959632221 -0.667602846  0.174106771  0.371190256 -0.859800730  1.230934718 
    ##   theta[139]   theta[140]   theta[141]   theta[142]   theta[143] 
    ## -0.330208235 -1.334808650  0.620353049 -0.181311390 -0.378247557

## 4.3 Posterior predcitive checks

The posterior predictive distributions can be visualized using
`ppCheck`. Set type = “response” to visualize the predicted probability
of a correct response / item acceptance as a function of item or person.

``` r
yrep <- posteriorPredict(samples, ndraws=20)
ppCheck(samples, type = "response", yrep=yrep)
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ppCheck(samples, type = "response", group = "item", yrep=yrep )
```

![](README_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
ppCheck(samples, type = "response", group = "person",
        index=1:10, yrep=yrep)
```

![](README_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

Set type = “rtQuantile” to compares observed and posterior-predictive
reaction-time quantiles:

``` r
ppCheck(samples, type = "rtQuantile")
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## 4.4 Plot item characteristics

According to drift diffusion item response theory, whether a person
solves a test item or accepts a survey item depends on two latent
variables, person ability $`\theta`$ and person response caution
$`\gamma`$. Thus, we can characterise the response to an item by surface
plot with the two latent variables on the x-axis and y-axis,
respectively, and the probability of solving / accepting as colours.
Warning: Creating response surfaces consumes some time, so please be
patient.

``` r
plotResponseSurface(samples, items = 1:10, facet.ncol = 5)
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

# 5 Contributing to the package

The package is under active development. Please feel free to [contact
us](mailto:manuel.rausch@aau.at) to suggest diffusion item response
theory models that we might have not yet implemented, or to volunteer
adding additional features.

# 6 Contact

For comments, bug reports, and feature suggestions please feel free to
either write to <manuel.rausch@aau.at> or [submit an
issue](https://github.com/ManuelRausch/BayesDiffIRT/issues).

# 7 References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-alexandrowicz_diffusion_2020" class="csl-entry">

Alexandrowicz, Rainer W. 2020. “The Diffusion Model Visualizer: An
Interactive Tool to Understand the Diffusion Model Parameters.”
*Psychological Research* 84 (4): 1157–65.
<https://doi.org/10.1007/s00426-018-1112-6>.

</div>

<div id="ref-carpenter_stan_2017" class="csl-entry">

Carpenter, Bob, Andrew Gelman, Matthew D. Hoffman, Daniel Lee, Ben
Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li,
and Allen Riddell. 2017. “*Stan* : A Probabilistic Programming
Language.” *Journal of Statistical Software* 76 (1).
<https://doi.org/10.18637/jss.v076.i01>.

</div>

<div id="ref-kang_modeling_2022" class="csl-entry">

Kang, Inhan, Paul De Boeck, and Roger Ratcliff. 2022. “Modeling
Conditional Dependence of Response Accuracy and Response Time with the
Diffusion Item Response Theory Model.” *Psychometrika* 87 (2): 725–48.
<https://doi.org/10.1007/s11336-021-09819-5>.

</div>

<div id="ref-link_sequential_1975" class="csl-entry">

Link, S. W., and R. A. Heath. 1975. “A Sequential Theory of
Psychological Discrimination.” *Psychometrika* 40 (1): 77–105.
<https://doi.org/10.1007/BF02291481>.

</div>

<div id="ref-molenaar_fitting_2015" class="csl-entry">

Molenaar, Dylan, Francis Tuerlinckx, and Han L. J. Van Der Maas. 2015.
“Fitting Diffusion Item Response Theory Models for Responses and
Response Times Using the *r* Package **diffIRT**.” *Journal of
Statistical Software* 66 (4). <https://doi.org/10.18637/jss.v066.i04>.

</div>

<div id="ref-Ratcliff2016" class="csl-entry">

Ratcliff, Roger, Philip L Smith, Scott D Brown, and Gail McKoon. 2016.
“Diffusion Decision Model : Current Issues and History.” *Trends in
Cognitive Sciences* 20 (4): 260–81.
<https://doi.org/10.1016/j.tics.2016.01.007>.

</div>

<div id="ref-stone_models_1960" class="csl-entry">

Stone, Mervyn. 1960. “Models for Choice-Reaction Time.” *Psychometrika*
25 (3): 251–60. <https://doi.org/10.1007/BF02289729>.

</div>

<div id="ref-tuerlinckx_two_2005" class="csl-entry">

Tuerlinckx, Francis, and Paul De Boeck. 2005. “Two Interpretations of
the Discrimination Parameter.” *Psychometrika* 70 (4): 629–50.
<https://doi.org/10.1007/s11336-000-0810-3>.

</div>

<div id="ref-van_der_maas_cognitive_2011" class="csl-entry">

Van Der Maas, Han L. J., Dylan Molenaar, Gunter Maris, Rogier A. Kievit,
and Denny Borsboom. 2011. “Cognitive Psychology Meets Psychometric
Theory: On the Relation Between Process Models for Decision Making and
Latent Variable Models for Individual Differences.” *Psychological
Review* 118 (2): 339–56. <https://doi.org/10.1037/a0022749>.

</div>

</div>
