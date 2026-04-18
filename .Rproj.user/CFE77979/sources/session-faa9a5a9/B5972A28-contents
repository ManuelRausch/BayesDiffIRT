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

normalize_priors <- function(prior) {
  if (is.null(prior)) {
    return(list())
  }
  if (inherits(prior, "BayesDiffIRTPrior")) {
    return(list(prior))
  }
  if (is.list(prior) && all(vapply(prior, inherits, logical(1), "BayesDiffIRTPrior"))) {
    return(prior)
  }
  stop("prior must be NULL, a BayesDiffIRTPrior, or a list of BayesDiffIRTPrior objects",
       call. = FALSE)
}

complete_priors <- function(prior, model = "d") {
  prior <- normalize_priors(prior)
  defaults <- default_priors(model = model)

  user_classes <- vapply(prior, function(x) x$class, character(1))
  default_classes <- vapply(defaults, function(x) x$class, character(1)) # extract the potentially available classes from the default because you don't know if you want to add a model later on with different defaults

  missing <- !(default_classes %in% user_classes)

  c(prior, defaults[missing])
}

default_priors <- function(model = "d") {
  if (model == "d"){
    # to do: Think hard about weakly informative values
    return(list(
      prior(normal(0, 1), "omega_theta"),
      prior(normal(0, 1), "omega_gamma"),
      prior(normal(0, 2), "nu"),
      prior(lognormal(0, 0.5), "a"),
      prior(lognormal(-1, 0.5), "tnd")))
  }
  if (model == "q"){
    # to do: Think hard about weakly informative values
    return(list(
      prior(normal(0, 1), "omega_theta"),
      prior(normal(0, 1), "omega_gamma"),
      prior(normal(0, 2), "nu"),
      prior(lognormal(0, 0.5), "a"),
      prior(lognormal(-1, 0.5), "tnd")))
  }
  stop("Model not recognized", call. = FALSE)
}
