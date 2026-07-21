#' Extract point estimates from a BayesDiffIRT model
#'
#' Extract posterior point estimates for the model parameters.
#'
#' @param object A fitted object of class \code{"BayesDiffIRTfit"}.
#' @param estimate Point estimate to return. Either \code{"mean"} or
#'   \code{"median"}.
#' @param parameters Character vector specifying the parameters to extract.
#'   By default, all substantive model parameters are returned.
#' @param ... Additional arguments passed to \code{summary()}.
#'
#' @return A named numeric vector containing posterior point estimates.
#'
#' @method coef BayesDiffIRTfit
#' @export

coef.BayesDiffIRTfit <- function(
    object,
    estimate = "mean",
    parameters = NULL,
    ...) {

  if (!inherits(object, "BayesDiffIRTfit")) {
    stop("`object` must inherit from class \"BayesDiffIRTfit\".",
         call. = FALSE)
  }

  modelSummary <- summary(object, ...)
  variables <- modelSummary$variables

  requiredColumns <- c("variable", estimate)

  # Remove the index from names such as theta[1].
  baseName <- sub("\\[.*$", "", variables$variable)

  if (is.null(parameters)) {
    parameters <-
      c("theta","gamma","nu","a","tnd","omega_theta","omega_gamma")
  }

  keep <- baseName %in% parameters

  result <- variables[[estimate]][keep]
  names(result) <- variables$variable[keep]

  result
}
