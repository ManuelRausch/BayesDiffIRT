#' @export
summary.BayesDiffIRTfit <- function(object, ...) {

  fitSummary <- object$fit$summary()
  fitSummary <-
    fitSummary[fitSummary$variable != "lp__", ,
                drop = FALSE]

  out <- list(
    call = object$call,
    model = object$model,
    dataInfo = list(
      nObs = object$stanData$nObs,
      nPerson = object$stanData$nPerson,
      nItem = object$stanData$nItem
    ),
    samplingInfo = object$fit$metadata(),
    variables = fitSummary
  )

  class(out) <- "summary.BayesDiffIRTfit"
  out
}

#' @export
print.summary.BayesDiffIRTfit <-
  function(x, digits = 2, nSubjects = NULL,
           nItems = NULL, nHypers = NULL, ...) {

    cat("Summary of BayesDiffIRT model fit\n")
    cat("---------------------------------\n")

    cat("Model:", switch(x$model,
                         "d" = "D-diffusion model",
                         "q" = "Q-diffusion model",
                         x$model), "\n")

    cat("\nData:\n")
    cat("  Observations:", x$dataInfo$nObs, "\n")
    cat("  Persons:     ", x$dataInfo$nPerson, "\n")
    cat("  Items:       ", x$dataInfo$nItem, "\n")

    cat("\nCall:\n")
    print(x$call)

    vars <- x$variables

    keep <- intersect(
      c("variable", "mean", "median", "sd",
        "q5", "q95", "q2.5", "q97.5",
        "rhat", "ess_bulk", "ess_tail"),
      names(vars)
    )

    vars <- vars[, keep, drop = FALSE]
    vars <- roundDf(vars, digits = digits)

    # classify parameters
    isHyper <- vars$variable %in% c("omega_theta", "omega_gamma")
    isItem <- grepl("^(nu|a)\\[", vars$variable)
    isSubject <- grepl("^(theta|gamma|tnd)\\[", vars$variable)

    varsHyper <- vars[isHyper, , drop = FALSE]
    varsItem <- vars[isItem, , drop = FALSE]
    varsSubject <- vars[isSubject, , drop = FALSE]
    #vars_other <- vars[!(is_hyper | is_item | is_subject),, drop = FALSE]

    headN <- function(df, n) {
      if (nrow(df) == 0) {
        return(df)
      }
      if (is.null(n) || n >= nrow(df)) {
        return(df)
      }
      df[seq_len(n), , drop = FALSE]
    }

    printSection <- function(title, df, n) {
      if (nrow(df) == 0) {
        return(invisible(NULL))
      }

      cat("\n", title, ":\n", sep = "")
      print(headN(df, n), row.names = FALSE)

      if (!is.null(n) && nrow(df) > n) {
        cat("... (", nrow(df) - n, " more)\n", sep = "")
      }

      invisible(NULL)
    }

    cat("\nPosterior summaries:\n")

    printSection("Hyperparameters", varsHyper, nHypers)
    printSection("Item parameters", varsItem, nItems)
    printSection("Subject parameters", varsSubject, nSubjects)

    invisible(x)
  }


roundDf <- function(df, digits = 2) {
  out <- df
  isNum <- vapply(out, is.numeric, logical(1))
  out[isNum] <- lapply(out[isNum], round, digits = digits)
  out
}
