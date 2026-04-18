<!-- Dieses RMarkdown Dokument dient dazu, ein Readme-Markdown-Dokument zu knitten.  Das Read-me wird dann auf Github angezeigt wird, wenn auf die Website des Pakets geht -->

# 1 The BayesDiffIRT for R

The `BayesDiffIRT` provides functions to sample posterior distributions
and posterior predictive distributions of item and subject parameters of
diffusion item response theory models for responses and reaction time
<!-- ToDo: Add references --> `BayesDiffIRT` also provides functions to
visualiize posterior distributions of Diffusion item response theory
model parameters and construct HDIs. Under the hood, the package relies
on NUTS sampling with STAN <!-- Citation missing -->. Up to know, the
following diffusion item response theory models have been implemented:

- D-Diffusion model <!--citation missing -->
- Q-Diffusion model

# 2 Mathematical description of diffusion item response theory models

All models included in the `BayesDiffIRT` package are all based on the
dirft diffusion model of decision making (Ratcliff 1978).

# 3 Contributing to the package

The package is under active development. We are planning to implement
new models of decision confidence when they are published. Please feel
free to [contact us](mailto:manuel.rausch@aau.at) to suggest new models
to implement in the package, or to volunteer adding additional models.

# 4 Contact

For comments, bug reports, and feature suggestions please feel free to
either write to <manuel.rausch@aau.at> or [submit an
issue](https://github.com/ManuelRausch/BayesDiffIRT/issues).

# 5 References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Ratcliff1978" class="csl-entry">

Ratcliff, Roger. 1978. “A Theory of Memory Retrieval.” *Psychological
Review* 85 (2): 59–108.

</div>

</div>
