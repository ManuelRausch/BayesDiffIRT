#' Check MCMC diagnostics
#'
#' Check convergence and sampling diagnostics for fitted BayesDiffIRT models.
#'
#' `checkDiagnostics()` is an S3 generic for extracting and summarizing MCMC
#' diagnostics from fitted model objects. The method for objects of class
#' `BayesDiffIRTfit` checks common Stan diagnostics, including divergent
#' transitions, high \eqn{\hat R}, low effective sample size, maximum treedepth
#' hits, and low E-BFMI.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return For objects of class `BayesDiffIRTfit`, an object of class
#'   `diagnostics.BayesDiffIRTfit`. The returned object is a list with elements:
#'   \describe{
#'     \item{\code{ok}}{Logical. \code{TRUE} if no potential diagnostic problems
#'       were detected.}
#'     \item{\code{nChains}}{Number of MCMC chains.}
#'     \item{\code{thresholds}}{Diagnostic thresholds used for the checks.}
#'     \item{\code{divergences}}{Divergent transitions by chain.}
#'     \item{\code{rhat}}{Parameters with \eqn{\hat R} values above the
#'       threshold.}
#'     \item{\code{ess}}{Parameters with low bulk or tail effective sample size.}
#'     \item{\code{treedepth}}{Maximum treedepth information by chain.}
#'     \item{\code{ebfmi}}{E-BFMI values and warnings by chain.}
#'   }
#'
#' @details
#' For `BayesDiffIRTfit` objects, the method extracts sampler diagnostics and
#' posterior summaries from the underlying Stan fit. The maximum treedepth is
#' read from the Stan fit metadata when available. If it cannot be recovered, a
#' default value of `15` is used.
#'
#' The default diagnostic thresholds are:
#' \describe{
#'   \item{\eqn{\hat R}}{Values larger than `1.01` are flagged.}
#'   \item{Effective sample size}{Bulk or tail effective sample size smaller
#'     than `100` per chain is flagged.}
#'   \item{E-BFMI}{Values smaller than `0.3` are flagged.}
#'   \item{Treedepth}{Transitions reaching the maximum treedepth are flagged.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fitBayesDiffIRT(data)
#'
#' diagnostics <- checkDiagnostics(fit)
#' diagnostics
#'
#' diagnostics$divergences
#' diagnostics$rhat
#' diagnostics$ess
#' diagnostics$treedepth
#' diagnostics$ebfmi
#' }
#'
#' @seealso
#' \code{\link{fitBayesDiffIRT}}
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

