---
title: 'BayesDiffIRT: An R Package for Bayesian estimation of Diffusion Item Response Theory Models for responses and response times'
tags:
- item response theory
- drift diffusion model
- Bayesian modelling
- R
- R package
date: "16 "
output:
  pdf_document: default
  word_document: default
  html_document:
    df_print: paged
authors:
- name: Manuel Rausch
  orcid: "0000-0002-5805-5544"
  affiliation: 1, 2
- name: Rainer W. Alexandrowicz
  orcid: "0000-0001-6928-4126"
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: "Alpe-Adria-Universität Klagenfurt, Institut für Psychologie, Abteilung für Methodenlehre, Klagenfurt, Austria"
  index: 1
- name: "Katholische Universität Eichstätt-Ingolstadt, Philosophisch-pädagogische
    Fakultät, Eichstätt, Germany"
  index: 2

---
  
# Summmary

The `BayesDiffIRT` package provides R functions to fits Bayesian drift diffusion
item-response theory models by sampling from the posterior distributions of item
and subject parameters using the No-U-Turn Sampler (NUTS) as implemented in Stan 
via `cmdstanr`[@carpenter_stan_2017]. Diffusion item response theory combines 
item response theory [@birnbaum_latent_1968, @rasch_probabilistic_1980] with the
drift diffusion  model of decision.making [@link_sequential_1975, 
@link_relative_1975, @Ratcliff2016, @alexandrowicz_diffusion_2020]. The drift 
diffusion model assumes that observers continuously accumulate evidence differentiating between
two binary decision options until the accumulated evidence hits one of the two decision 
boundaries, after which the corresponding choice option is selected. Drift diffusion
item-response theory models allow for the decomposition of two of the traditional parameters
of the drift diffusion decision model, boundary separation and drift rate, into person and item parameters.
The following diffusion item response theory models are currently implemented: 

* D-Diffusion model [@tuerlinckx_two_2005]
* Q-Diffusion model [@van_der_maas_cognitive_2011]. 

The packages allows for sampling for the posterior distribution, provies fit diagnostics, 
posterior predictive distributions, 

plot methods enable researchers to plot Marcov Chains, posterior distributions of selected parameters, 
as well as credible intervals. posterior predictive distributions. Diagnostic checks .
Other stuff that is IRT-realted?0

# Statement of need

The diffIRT package provides R functions to fit the Q-diffusion and the D-Diffusion
model using maximim likelihood estimation [@molenaar_2015]. 

# Statement of the field

# Software design

# Research impact statement

## Contact

For comments, bug reports, and feature suggestions please feel free to either write to 
[manuel.rausch@aau.at](mailto:manuel.rausch@aau.at) or [submit an issue](https://github.com/ManuelRausch/BayesDiffIRT/issues).

# Acknowledgements
    
# References
