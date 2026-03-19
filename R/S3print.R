#' @export
print.BayesDiffIRTfit <- function(x, ...) {
  cat("<DiffIRTfit>\n")
  cat("Model:", x$model$name, "\n")
  invisible(x)
}
