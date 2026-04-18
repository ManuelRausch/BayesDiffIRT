<!-- Dieses RMarkdown Dokument dient dazu, ein Readme-Markdown-Dokument zu knitten.  Dieses Read-me wird man sehen, wenn man auf die Paket Development Page of Github geht.  -->

# 1 The BayesDiffIRT for R

The `BayesDiffIRT` provides functions to sample posterior distributions
and posterior predictive distributions of item and subject parameters of
diffusion item response theory models for responses and reaction times
Tuerlinckx and Boeck (2005) `BayesDiffIRT` also provides functions to
visualize posterior distributions of Diffusion item response theory
model parameters and construct HDIs. Under the hood, the package relies
on NUTS sampling with STAN (**carpenter_stan_2017?**). Up to know, the
following diffusion item response theory models have been implemented:

- D-Diffusion model (Tuerlinckx and Boeck 2005),
- Q-Diffusion model (Van Der Maas et al. 2011).

# 2 Mathematical description of diffusion item response theory models

All models included in the `BayesDiffIRT` package are all based on the
drift diffusion model of decision making (**stone_models_1960?**; Link
and Heath 1975; Ratcliff et al. 2016).

To visualize the resulting reaction time distributions, interactive
tools such as the [diffusion model visualizer](https://osf.io/4en3b) are
openly available (Alexandrowicz 2020).

# 3 Contributing to the package

The package is under active development. Please feel free to [contact
us](mailto:manuel.rausch@aau.at) to suggest diffusion item response
theory models that we might have not yet implemented, or to volunteer
adding additional features.

# 4 Contact

For comments, bug reports, and feature suggestions please feel free to
either write to <manuel.rausch@aau.at> or [submit an
issue](https://github.com/ManuelRausch/BayesDiffIRT/issues).

# 5 References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-alexandrowicz_diffusion_2020" class="csl-entry">

Alexandrowicz, Rainer W. 2020. “The Diffusion Model Visualizer: An
Interactive Tool to Understand the Diffusion Model Parameters.”
*Psychological Research* 84 (4): 1157–65.
<https://doi.org/10.1007/s00426-018-1112-6>.

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
