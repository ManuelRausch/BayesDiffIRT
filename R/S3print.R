#' @export
print.DiffIRT <- function(x, ...) {
  cat("<DiffIRT>\n")
  cat("Model:", x$model$name, "\n")
  invisible(x)
}
