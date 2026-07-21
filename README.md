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

    ## Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup)

    ## Chain 3 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 1 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 2 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 4 Iteration:  200 / 2000 [ 10%]  (Warmup) 
    ## Chain 3 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 1 Iteration:  400 / 2000 [ 20%]  (Warmup)

    ## Chain 2 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 4 Iteration:  400 / 2000 [ 20%]  (Warmup) 
    ## Chain 3 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 2 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 1 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 4 Iteration:  600 / 2000 [ 30%]  (Warmup) 
    ## Chain 3 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 1 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 2 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 4 Iteration:  800 / 2000 [ 40%]  (Warmup) 
    ## Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
    ## Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
    ## Chain 3 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 2 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 1 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 4 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
    ## Chain 3 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 2 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 1 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 4 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
    ## Chain 3 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 2 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 1 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 4 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
    ## Chain 3 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 2 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 1 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 4 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
    ## Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 3 finished in 209.3 seconds.
    ## Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 2 finished in 214.5 seconds.
    ## Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 4 finished in 214.9 seconds.
    ## Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
    ## Chain 1 finished in 210.8 seconds.
    ## 
    ## All 4 chains finished successfully.
    ## Mean chain execution time: 212.4 seconds.
    ## Total execution time: 217.2 seconds.

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
    ## 1 omega_theta  0.69   0.69  0.06  0.6   0.8      1    1109.    1589.
    ## 2 omega_gamma  0.2    0.2   0.03  0.15  0.25     1     725.    1464.
    ## 
    ## Item parameters:
    ## # A tibble: 20 × 9
    ##    variable  mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 nu[1]    -0.65  -0.65  0.11 -0.82 -0.48     1    1486.    2497.
    ##  2 nu[2]    -0.15  -0.15  0.11 -0.32  0.03     1    1605.    2279.
    ##  3 nu[3]    -1.23  -1.23  0.13 -1.44 -1.03     1    1783.    2943.
    ##  4 nu[4]    -1.7   -1.7   0.15 -1.94 -1.46     1    1705.    2632.
    ##  5 nu[5]    -0.21  -0.21  0.11 -0.39 -0.03     1    1633.    2463.
    ##  6 nu[6]    -1.3   -1.3   0.12 -1.51 -1.09     1    1733.    2594.
    ##  7 nu[7]    -1.69  -1.69  0.15 -1.93 -1.45     1    1999.    2621.
    ##  8 nu[8]    -1.91  -1.91  0.15 -2.16 -1.67     1    2038.    2553.
    ##  9 nu[9]    -0.83  -0.83  0.11 -1    -0.66     1    1306.    2561.
    ## 10 nu[10]   -1.42  -1.42  0.14 -1.65 -1.2      1    2113.    2612.
    ## 11 a[1]      0.44   0.44  0.02  0.41  0.47     1    1891.    2533.
    ## 12 a[2]      0.49   0.49  0.02  0.46  0.53     1    2143.    2732.
    ## 13 a[3]      0.5    0.49  0.02  0.46  0.54     1    1682.    3054.
    ## 14 a[4]      0.51   0.51  0.03  0.47  0.56     1    1508.    2357.
    ## 15 a[5]      0.51   0.51  0.02  0.48  0.55     1    1921.    2418.
    ## 16 a[6]      0.43   0.43  0.02  0.4   0.47     1    2117.    2981.
    ## 17 a[7]      0.4    0.4   0.02  0.37  0.45     1    2080.    2672.
    ## 18 a[8]      0.42   0.42  0.02  0.38  0.46     1    1905.    2425.
    ## 19 a[9]      0.35   0.35  0.02  0.32  0.38     1    2569.    2313.
    ## 20 a[10]     0.55   0.55  0.03  0.5   0.6      1    1779.    2337.
    ## 
    ## Subject parameters:
    ## # A tibble: 429 × 9
    ##    variable  mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 tnd[1]    0.37   0.36  0.1   0.22  0.54     1    4895.    3073.
    ##  2 tnd[2]    0.4    0.4   0.09  0.24  0.54     1    3672.    2885.
    ##  3 tnd[3]    0.46   0.47  0.12  0.27  0.65     1    3868.    3010.
    ##  4 tnd[4]    0.27   0.27  0.06  0.18  0.37     1    4644.    2926.
    ##  5 tnd[5]    0.39   0.39  0.08  0.26  0.51     1    4105.    3447.
    ##  6 tnd[6]    0.31   0.32  0.06  0.21  0.41     1    4428.    3376.
    ##  7 tnd[7]    0.39   0.39  0.09  0.24  0.53     1    4004.    2712.
    ##  8 tnd[8]    0.52   0.52  0.15  0.27  0.78     1    3288.    2695.
    ##  9 tnd[9]    0.34   0.34  0.07  0.21  0.46     1    4880.    2951.
    ## 10 tnd[10]   0.43   0.42  0.13  0.24  0.65     1    4587.    2806.
    ## # ℹ 419 more rows
    ## 
    ## Other parameters:
    ## # A tibble: 429 × 9
    ##    variable     mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>       <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 z_theta[1]  -0.82  -0.82  0.36 -1.4  -0.24     1    3458.    2830.
    ##  2 z_theta[2]  -0.3   -0.31  0.54 -1.17  0.62     1    5767.    2674.
    ##  3 z_theta[3]   0.4    0.39  0.54 -0.45  1.32     1    5866.    3116.
    ##  4 z_theta[4]   0.49   0.48  0.46 -0.26  1.26     1    4103.    3000.
    ##  5 z_theta[5]   1.47   1.48  0.67  0.4   2.59     1    5223.    3058.
    ##  6 z_theta[6]  -0.07  -0.07  0.53 -0.96  0.8      1    5785.    3105.
    ##  7 z_theta[7]  -0.07  -0.07  0.52 -0.91  0.8      1    5294.    2862.
    ##  8 z_theta[8]   0.33   0.32  0.56 -0.59  1.28     1    4870.    3072.
    ##  9 z_theta[9]   0.27   0.27  0.49 -0.53  1.07     1    5256.    2935.
    ## 10 z_theta[10]  0.37   0.36  0.46 -0.35  1.13     1    4306.    2974.
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
    ## 1 -1729      -0.81      -0.37       0.84      0.309       0.51      0.994
    ## 2 -1740      -1.05      -0.25       0.34      0.444       2.60     -1.189
    ## 3 -1733      -0.57      -0.49       0.31      0.035       0.76      0.702
    ## 4 -1711      -0.76       0.18       1.17      0.665       1.80     -0.443
    ## 5 -1713      -1.06      -1.03      -0.29      0.071       2.16     -0.249
    ## 6 -1706      -0.96      -1.03      -0.16      0.735       2.37     -0.026
    ##   z_theta[7]
    ## 1      -0.20
    ## 2       0.60
    ## 3      -0.50
    ## 4       0.61
    ## 5      -0.30
    ## 6      -0.71
    ## # ... with 873 more variables
    ## # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

Should you really prefer working with point estimates, you can extract
them conveniently with the coef method:

``` r
pointEstim <- coef(samples, parameter = "theta")
pointEstim 
```

    ##     theta[1]     theta[2]     theta[3]     theta[4]     theta[5]     theta[6] 
    ## -0.569520019 -0.209285584  0.277636153  0.336519997  1.023910705 -0.049829556 
    ##     theta[7]     theta[8]     theta[9]    theta[10]    theta[11]    theta[12] 
    ## -0.051630431  0.230456798  0.184561586  0.255965047 -0.035412350  0.347768629 
    ##    theta[13]    theta[14]    theta[15]    theta[16]    theta[17]    theta[18] 
    ##  0.291585559 -0.048192310  0.502675284  0.317207223 -0.112446303 -0.050642397 
    ##    theta[19]    theta[20]    theta[21]    theta[22]    theta[23]    theta[24] 
    ## -0.248443253  0.666722577 -1.109748113  1.226674588  0.454476878  0.089180505 
    ##    theta[25]    theta[26]    theta[27]    theta[28]    theta[29]    theta[30] 
    ##  0.129068258 -0.925067323  1.103761641  0.233593194  0.527437541 -0.780084583 
    ##    theta[31]    theta[32]    theta[33]    theta[34]    theta[35]    theta[36] 
    ##  0.428069266 -0.864653524  0.717329926  0.041588960  0.084318169  0.153425651 
    ##    theta[37]    theta[38]    theta[39]    theta[40]    theta[41]    theta[42] 
    ## -0.139378609  0.263625318 -0.037596780  0.336340478 -1.107169604 -0.371515273 
    ##    theta[43]    theta[44]    theta[45]    theta[46]    theta[47]    theta[48] 
    ##  0.736954869  0.403726323  0.638236034 -0.042707669  0.033432663 -0.283141711 
    ##    theta[49]    theta[50]    theta[51]    theta[52]    theta[53]    theta[54] 
    ## -0.625372419  0.411520380  0.407479429 -0.207486155  1.190459440 -0.046561946 
    ##    theta[55]    theta[56]    theta[57]    theta[58]    theta[59]    theta[60] 
    ##  0.005897704  1.159979454  0.406145197  0.423256696 -1.159252943 -1.048539675 
    ##    theta[61]    theta[62]    theta[63]    theta[64]    theta[65]    theta[66] 
    ## -1.100214038 -0.078644509 -0.482871591  1.120388319 -0.129883298 -0.695972741 
    ##    theta[67]    theta[68]    theta[69]    theta[70]    theta[71]    theta[72] 
    ## -0.711958675  0.163137250  0.892144985  0.792333848 -0.392710605  0.501243717 
    ##    theta[73]    theta[74]    theta[75]    theta[76]    theta[77]    theta[78] 
    ## -0.393733304  0.991091822 -0.138492906 -1.681032655  0.270416812 -0.684742756 
    ##    theta[79]    theta[80]    theta[81]    theta[82]    theta[83]    theta[84] 
    ## -0.266418176 -0.189490327 -0.313077470 -0.498725094  0.412683318  0.132127167 
    ##    theta[85]    theta[86]    theta[87]    theta[88]    theta[89]    theta[90] 
    ##  0.362474891 -0.088304637  0.042290239 -0.398056927 -0.194056580  0.283639296 
    ##    theta[91]    theta[92]    theta[93]    theta[94]    theta[95]    theta[96] 
    ## -0.025972026 -1.347483643  0.314232349 -0.374185652 -0.210352702 -0.494349631 
    ##    theta[97]    theta[98]    theta[99]   theta[100]   theta[101]   theta[102] 
    ## -0.136139282 -0.338163725 -0.941078150 -0.649550482  0.134793934  0.905775589 
    ##   theta[103]   theta[104]   theta[105]   theta[106]   theta[107]   theta[108] 
    ## -0.212478751 -0.322475715  0.008521914  0.088084597  0.379291563 -0.147629974 
    ##   theta[109]   theta[110]   theta[111]   theta[112]   theta[113]   theta[114] 
    ##  0.226490632 -0.788087726 -0.191515119  0.534028985 -0.578334130  1.058881166 
    ##   theta[115]   theta[116]   theta[117]   theta[118]   theta[119]   theta[120] 
    ##  0.041381601  0.042342203 -0.395949337 -0.422913056  0.029548343  0.171584846 
    ##   theta[121]   theta[122]   theta[123]   theta[124]   theta[125]   theta[126] 
    ## -1.015805697  0.802383888  0.523905988 -1.383991773  0.074167739  0.782436970 
    ##   theta[127]   theta[128]   theta[129]   theta[130]   theta[131]   theta[132] 
    ##  0.457817822  0.064473234  0.008655518  0.612048040  0.153635861  0.480023626 
    ##   theta[133]   theta[134]   theta[135]   theta[136]   theta[137]   theta[138] 
    ##  0.956621014 -0.669607382  0.173875118  0.374055617 -0.857826497  1.235726873 
    ##   theta[139]   theta[140]   theta[141]   theta[142]   theta[143] 
    ## -0.334005857 -1.332447561  0.618582328 -0.187328022 -0.394417040

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
