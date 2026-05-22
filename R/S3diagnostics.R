#' Run diagnostic checks for BayesDiffIRTfit objects.
#'
#' @param object A fitted BayesDiffIRTfit object.
#' @param ... Further arguments passed to methods.
#'
#' @export
checkDiagnostics <- function(object, ...) {
  UseMethod("checkDiagnostics")
}

#' @export
checkDiagnostics.BayesDiffIRTfit <- function(object) {
  md <- object$fit$metadata()
  if (!is.null(md$max_treedepth) && !is.na(md$max_treedepth)){
    maxTreedepth <- md$max_treedepth
  } else {
    maxTreedepth = 15
  }

  checkStanDiagnostics(object$fit, maxTreedepth)
}


checkStanDiagnostics <-
  function(fit,
           max.treedepth,
           rhat.threshold = 1.01,
           ess.per.chain = 100,
           ebfmi.threshold = 0.3) {
    summ <- fit$summary()
    diag <- fit$sampler_diagnostics(format = "df")

    nChains <- length(unique(diag$.chain))
    essThreshold <- ess.per.chain * nChains

    # Divergences by chain

    chains <- sort(unique(diag$.chain))

    divergences <- do.call(
      rbind,
      lapply(chains, function(ch) {
        x <- diag[diag$.chain == ch, , drop = FALSE]

        data.frame(
          chain = ch,
          nDivergent = sum(x$divergent__, na.rm = TRUE),
          nTransitions = nrow(x),
          propDivergent = mean(x$divergent__, na.rm = TRUE)
        )
      })
    )

    rownames(divergences) <- NULL

    # R-hat problems

    rhatBad <- summ[
      !is.na(summ$rhat) & summ$rhat > rhat.threshold,
      c("variable", "rhat"),
      drop = FALSE
    ]

    if (nrow(rhatBad) > 0) {
      rhatBad <- rhatBad[order(rhatBad$rhat, decreasing = TRUE), , drop = FALSE]
      rownames(rhatBad) <- NULL
    }

    # effetive sample size
    essBad <- summ[
      (!is.na(summ$ess_bulk) & summ$ess_bulk < essThreshold) |
        (!is.na(summ$ess_tail) & summ$ess_tail < essThreshold),
      c("variable", "ess_bulk", "ess_tail"),
      drop = FALSE
    ]

    if (nrow(essBad) > 0) {
      essOrder <- pmin(essBad$ess_bulk, essBad$ess_tail, na.rm = TRUE)
      essBad <- essBad[order(essOrder), , drop = FALSE]
      rownames(essBad) <- NULL
    }

    # Treedepth

    treedepth <- do.call(
      rbind,
      lapply(chains, function(ch) {
        x <- diag[diag$.chain == ch, , drop = FALSE]

        data.frame(
          chain = ch,
          maxTreedepthObserved = max(x$treedepth__, na.rm = TRUE),
          nMaxTreedepth = sum(x$treedepth__ >= max.treedepth, na.rm = TRUE),
          propMaxTreedepth = mean(x$treedepth__ >= max.treedepth, na.rm = TRUE)
        )
      })
    )
    rownames(treedepth) <- NULL

    # E-BFMI

    ebfmi <- do.call(
      rbind,
      lapply(chains, function(ch) {
        x <- diag[diag$.chain == ch, , drop = FALSE]

        if (".iteration" %in% names(x)) {
          x <- x[order(x$.iteration), , drop = FALSE]
        }

        energy <- x$energy__

        val <- mean(diff(energy)^2, na.rm = TRUE) / stats::var(energy, na.rm = TRUE)

        data.frame(
          chain = ch,
          ebfmi = val,
          warning = is.na(val) || val < ebfmi.threshold
        )
      })
    )

    rownames(ebfmi) <- NULL

    # Summary

    ok <- sum(divergences$nDivergent, na.rm = TRUE) == 0 &&
      nrow(rhatBad) == 0 &&
      nrow(essBad) == 0 &&
      !any(ebfmi$warning, na.rm = TRUE) &&
      sum(treedepth$nMaxTreedepth, na.rm = TRUE) == 0

    out <- list(
      ok = ok,
      nChains = nChains,
      thresholds = list(
        rhat = rhat.threshold,
        ess = essThreshold,
        essPerChain = ess.per.chain,
        ebfmi = ebfmi.threshold,
        maxTreedepth = max.treedepth
      ),
      divergences = divergences,
      rhat = rhatBad,
      ess = essBad,
      treedepth = treedepth,
      ebfmi = ebfmi
    )

    class(out) <- "diagnostics.BayesDiffIRTfit"
    out
  }

#' @export
print.diagnostics.BayesDiffIRTfit <- function(x, ...) {
  cat("Stan diagnostics\n")
  cat("----------------\n")

  cat("Chains:               ", x$nChains, "\n", sep = "")
  cat("Divergences:          ", sum(x$divergences$nDivergent), "\n", sep = "")
  cat("Max treedepth hits:   ", sum(x$treedepth$nMaxTreedepth), "\n", sep = "")
  cat("R-hat warnings:       ", nrow(x$rhat), "\n", sep = "")
  cat("Low ESS warnings:     ", nrow(x$ess), "\n", sep = "")
  cat("E-BFMI warnings:      ", sum(x$ebfmi$warning), "\n", sep = "")

  cat("\nOverall status:        ")
  if (isTRUE(x$ok)) {
    cat("OK\n")
  } else {
    cat("potential problems detected\n")
  }

  invisible(x)
}

