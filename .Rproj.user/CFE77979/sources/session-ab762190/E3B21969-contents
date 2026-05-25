#' @title Prior Definitions for **BayesDiffIRTfit** models
#' @description Define priors for specific parameters or classes of parameters.
#' @param dist  An expression defining a distribution in **Stan** language (see examples)
#' @param class `character` The parameter class. The following parameter classes are available:
#' * \code{omega_theta} standard deviation of \theta, the latent trait.
#' * \code{omega_gamma} standard deviation of \gamma, the response caution parameter.
#' * \code{nu} item difficulty,
#' * \code{a} item time pressure,
#' * \code{tnd} person-specific non-decision time.
#' @param  description coef currently not implemented.


#' @export
prior <- function(dist, class, coef = NULL) {
  if (!class %in% c("omega_theta",
                    "omega_gamma",
                    "nu", "a", "tnd")){
    stop(paste0(class, " not recognized. class should be on of the following: omega_theta,
                    omega_gamma, nu, a, tnd"))
  }
  if (!is.null(coef)) stop("Setting individual priors for subjects and items is not yet supported")
  structure(
    list(
      dist = substitute(dist),
      class = class
    ),
    class = "BayesDiffIRTPrior"
  )
}

normalizePriors <- function(priors) {
  if (is.null(priors)) {
    return(list())
  }
  if (inherits(priors, "BayesDiffIRTPrior")) {
    return(list(priors))
  }
  if (is.list(priors) &&
      all(vapply(priors, inherits, logical(1), "BayesDiffIRTPrior"))) {
    return(priors)
  }
  stop("priors must be NULL, a BayesDiffIRTPrior, or a list of BayesDiffIRTPrior objects",
       call. = FALSE)
}

completePriors <- function(priors, model = "d") {
  priors <- normalizePriors(priors)
  defaults <- defaultPriors(model = model)

  userClasses <-
    vapply(priors, function(x) x$class, character(1))
  defaultClasses <-
    vapply(defaults, function(x) x$class, character(1)) # extract the potentially available classes from the default because you don't know if you want to add a model later on with different defaults

  missing <- !(defaultClasses %in% userClasses)

  c(priors, defaults[missing])
}

defaultPriors <- function(model = "d") {
  if (model == "d"){
    # to do: Think hard about weakly informative values
    return(list(
      prior(normal(0, 2.5), "omega_theta"),
      prior(normal(0, 0.5), "omega_gamma"),
      prior(normal(0, 2.5), "nu"),
      prior(lognormal(0, .25), "a"),
      prior(lognormal(-1.25, 0.3), "tnd")))
  }
  if (model == "q"){
    # to do: Think hard about weakly informative values
    return(list(
      prior(normal(0, 1), "omega_theta"),
      prior(normal(0, 1), "omega_gamma"),
      prior(normal(.01, 100), "nu"),
      prior(uniform(0.1, 0.5), "a"),
      prior(lognormal(0, 1), "tnd")))
  }
  stop("Model not recognized", call. = FALSE)
}

