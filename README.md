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
    ## 1 omega_theta  0.69   0.69  0.06  0.6   0.8      1    1382.    2784.
    ## 2 omega_gamma  0.2    0.2   0.03  0.16  0.25     1     860.    1508.
    ## 
    ## Item parameters:
    ## # A tibble: 20 × 9
    ##    variable  mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 nu[1]    -0.65  -0.64  0.1  -0.82 -0.48     1    1531.    2138.
    ##  2 nu[2]    -0.15  -0.14  0.11 -0.32  0.03     1    1632.    2309.
    ##  3 nu[3]    -1.24  -1.23  0.13 -1.45 -1.03     1    1914.    2462.
    ##  4 nu[4]    -1.7   -1.7   0.14 -1.94 -1.46     1    2063.    2485.
    ##  5 nu[5]    -0.21  -0.21  0.11 -0.39 -0.03     1    1643.    2458.
    ##  6 nu[6]    -1.3   -1.3   0.12 -1.5  -1.11     1    1909.    2562.
    ##  7 nu[7]    -1.69  -1.69  0.14 -1.92 -1.46     1    2075.    2597.
    ##  8 nu[8]    -1.91  -1.92  0.14 -2.15 -1.68     1    2339.    2710.
    ##  9 nu[9]    -0.83  -0.83  0.1  -1    -0.66     1    1622.    1889.
    ## 10 nu[10]   -1.42  -1.42  0.14 -1.65 -1.2      1    2340.    2667.
    ## 11 a[1]      0.44   0.44  0.02  0.41  0.47     1    2583.    2448.
    ## 12 a[2]      0.5    0.5   0.02  0.46  0.53     1    2204.    2528.
    ## 13 a[3]      0.5    0.5   0.02  0.46  0.54     1    2063.    2975.
    ## 14 a[4]      0.51   0.51  0.03  0.47  0.56     1    1927.    2528.
    ## 15 a[5]      0.51   0.51  0.02  0.48  0.55     1    2127.    2075.
    ## 16 a[6]      0.43   0.43  0.02  0.4   0.47     1    2261.    2255.
    ## 17 a[7]      0.4    0.4   0.02  0.37  0.44     1    2490.    3025.
    ## 18 a[8]      0.42   0.42  0.02  0.38  0.46     1    2275.    3036.
    ## 19 a[9]      0.35   0.35  0.02  0.32  0.38     1    2835.    2670.
    ## 20 a[10]     0.55   0.55  0.03  0.5   0.6      1    1960.    3107.
    ## 
    ## Subject parameters:
    ## # A tibble: 429 × 9
    ##    variable  mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>    <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 tnd[1]    0.37   0.36  0.1   0.22  0.54     1    5419.    3416.
    ##  2 tnd[2]    0.4    0.4   0.09  0.25  0.54     1    3918.    3306.
    ##  3 tnd[3]    0.46   0.46  0.12  0.27  0.65     1    4303.    3422.
    ##  4 tnd[4]    0.27   0.27  0.06  0.18  0.36     1    5046.    2665.
    ##  5 tnd[5]    0.39   0.39  0.08  0.26  0.51     1    4206.    3282.
    ##  6 tnd[6]    0.32   0.32  0.06  0.21  0.41     1    4078.    2588.
    ##  7 tnd[7]    0.39   0.39  0.09  0.24  0.53     1    4374.    2987.
    ##  8 tnd[8]    0.51   0.51  0.15  0.27  0.77     1    3666.    3438.
    ##  9 tnd[9]    0.34   0.34  0.07  0.22  0.45     1    5381.    3265.
    ## 10 tnd[10]   0.42   0.41  0.13  0.24  0.64     1    4201.    2723.
    ## # ℹ 419 more rows
    ## 
    ## Other parameters:
    ## # A tibble: 429 × 9
    ##    variable     mean median    sd    q5   q95  rhat ess_bulk ess_tail
    ##    <chr>       <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>
    ##  1 z_theta[1]  -0.82  -0.82  0.36 -1.41 -0.23     1    3798.    2797.
    ##  2 z_theta[2]  -0.3   -0.31  0.53 -1.17  0.57     1    5354.    2841.
    ##  3 z_theta[3]   0.39   0.39  0.53 -0.47  1.27     1    5618.    2559.
    ##  4 z_theta[4]   0.49   0.48  0.46 -0.22  1.27     1    4994.    3233.
    ##  5 z_theta[5]   1.47   1.47  0.69  0.34  2.61     1    6581.    3105.
    ##  6 z_theta[6]  -0.07  -0.07  0.54 -0.97  0.81     1    6293.    3274.
    ##  7 z_theta[7]  -0.08  -0.07  0.51 -0.93  0.78     1    5113.    2857.
    ##  8 z_theta[8]   0.32   0.3   0.56 -0.57  1.27     1    5028.    3248.
    ##  9 z_theta[9]   0.25   0.27  0.47 -0.51  1.03     1    5371.    3100.
    ## 10 z_theta[10]  0.38   0.37  0.44 -0.33  1.11     1    4178.    2971.
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
    ## 1 -1744      -0.84      -1.50      -0.66       0.43       2.33       0.14
    ## 2 -1705      -1.46      -0.19       1.38       0.56       0.22      -0.50
    ## 3 -1707      -0.87      -1.41       1.07      -0.14       1.68       0.82
    ## 4 -1724      -0.71       0.90      -0.17       0.65       2.07      -1.21
    ## 5 -1759      -0.50      -1.34       0.80       0.33       0.49      -0.24
    ## 6 -1719      -0.49      -1.01       1.01       0.27       0.88      -0.50
    ##   z_theta[7]
    ## 1       0.87
    ## 2      -0.15
    ## 3       0.21
    ## 4       0.14
    ## 5       0.62
    ## 6      -0.55
    ## # ... with 873 more variables
    ## # ... hidden reserved variables {'.chain', '.iteration', '.draw'}

Should you really prefer working with point estimates, you can extract
them conveniently with the coef method:

``` r
pointEstim <- coef(samples)
pointEstim 
```

    ## NULL

## 4.3 Posterior predcitive checks

The posterior predictive distributions can be visualized using
`ppCheck`. Set type = “response” to visualize the predicted probability
of a correct response / item acceptance as a function of item or person.

``` r
ppCheck(samples, type = "response")
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ppCheck(samples, type = "response", group = "item")
```

![](README_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
ppCheck(samples, type = "response", group = "person",
        index=1:10)
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

``` r
 plotResponseSurface(samples, item = 1)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
 plotResponseSurface(samples, item = 2)
```

![](README_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
 plotResponseSurface(samples, item = 3)
```

![](README_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

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
