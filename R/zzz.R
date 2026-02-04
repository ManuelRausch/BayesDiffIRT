#' BayesDiffIRT: Bayesian Sampling of Diffusion Item Response Theory Models for Responses and Response Times
#'
#' The BayesDiffIRT package provides functions for simulating, fitting,
#' and visualizing hierarchical confidence models.
#'
#' @docType package
#' @name BayesDiffIRT
#' @keywords internal
"_BayesDiffIRT"


.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "diffirt: Bayesian diffusion-IRT models with RTs.\n",
    "Backend: cmdstanr (run diffirt_setup() if needed)."
  )
}
