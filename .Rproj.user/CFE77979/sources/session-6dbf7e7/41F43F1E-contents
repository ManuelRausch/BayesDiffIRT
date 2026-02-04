#' @export
summary.DiffIRTfit <- function(object, ...) {
  validate_diffirt_fit(object)
  # Return a summary object (often a list) with its own class
  out <- list(
    model = object$model,
    diagnostics = object$diagnostics,
    # e.g. posterior::summarise_draws(as_draws(object), ...)
  )
  class(out) <- "summary.DiffIRTfit"
  out
}

#' @export
print.summary.DiffIRTfit <- function(x, ...) {
  cat("<summary.diffirt_fit>\n")
  cat("Model:", x$model$name, "\n")
  # print diagnostics, parameter summaries etc.
  invisible(x)
}
