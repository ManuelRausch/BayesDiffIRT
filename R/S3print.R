#' @export
print.BayesDiffIRTfit <- function(x, ...) {

  cat("BayesDiffIRT model fit\n")
  cat("-----------------------\n")

  # Model
  cat("Model:", switch(x$model,
                       "d" = "D-Diffusion model",
                       "q" = "Q-diffusion model",
                       "dRV" = "D-diffusion model with random variation",
                       "qRV" = "Q-diffusion model with random variation",
                       x$model), "\n")

  # Data
  cat("\nData:\n")
  cat("  Observations:", x$stan_data$nObs, "\n")
  cat("  Persons:     ", x$stan_data$nPerson, "\n")
  cat("  Items:       ", x$stan_data$nItem, "\n")

  # Sampling info
  fit <- x$fit
  meta <- fit$metadata()

  cat("\nSampling:\n")
  cat("  Chains:      ", meta$chains, "\n")
  cat("  Iter (warmup):", meta$iter_warmup, "\n")
  cat("  Iter (sample):", meta$iter_sampling, "\n")

  # Call
  cat("\nCall:\n")
  print(x$call)

  cat("\nUse summary() for posterior summaries and diagnostics() for detailed diagnostics.\n")

  invisible(x)
}
