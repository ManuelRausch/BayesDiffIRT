#' @export
extractSamples <- function(object, ...) {
  UseMethod("extractSamples")
}

#' @export
extractSamples.BayesDiffIRTfit <- function(object) {
  posterior::as_draws_df(object$fit$draws())
}
