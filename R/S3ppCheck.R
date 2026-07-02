#' Posterior predictive checks
#'
#' @param object A fitted model object or posterior predictive object.
#' @param ... Additional arguments passed to methods.
#'
#' @export
ppCheck <- function(object, ...) {
  UseMethod("ppCheck")
}

#' @rdname ppCheck
#' @method ppCheck BayesDiffIRTfit
#' @export
ppCheck.BayesDiffIRTfit <-
  function(object,
           type = c( "response"),
           ndraws = 10,
           group = c("none", "item", "person"),
           seed = NULL,
           ...)
  {
    type <- match.arg(type)
    group <- match.arg(group)

    yrep <- posteriorPredict(object, ndraws = ndraws, seed = seed)

    if (type == "response") {
      return(ppCheckResponse(object, ndraws = ndraws, group = group, seed = seed))
    }

    stop("PPC type not implemented.", call. = FALSE)
  }


ppCheckResponse <- function(object,
                            ndraws = 10,
                            group = c("item", "person"),
                            seed = NULL) {

  group <- match.arg(group)

  yrep <- posteriorPredict(object, ndraws = ndraws, seed = seed)

  if (is.null(object$stanData$resp)) {
    stop("object$stanData$resp is required for response PPCs.", call. = FALSE)
  }

  obs <- data.frame(
    obs  = seq_len(object$stanData$nObs),
    sbj  = object$stanData$person,
    item = object$stanData$item,
    resp = object$stanData$resp
  )

  if (group == "item") {

    obsStat <- aggregate(resp ~ item, obs, mean)
    names(obsStat)[2] <- "observed"

    predStat <- aggregate(resp ~ draw + item,
                          data = yrep, mean)

    predSummary <- aggregate(resp ~ item, data = predStat,
                             FUN = function(x) {
                               c(mean = mean(x),
                                 lo50 = stats::quantile(x, 0.25),
                                 hi50 = stats::quantile(x, 0.75),
                                 lo95 = stats::quantile(x, 0.025),
                                 hi95 = stats::quantile(x, 0.975))
                             })

    predSummary <- data.frame(item = predSummary$item,
                              predSummary$resp,row.names = NULL)

    plotDf <- merge(obsStat, predSummary, by = "item")

    p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = item)) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = `lo95.2.5.`, ymax = `hi95.97.5.`),
        linewidth = 0.4) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = `lo50.25.`, ymax = `hi50.75.`),
        linewidth = 1.1) +
      ggplot2::geom_point(ggplot2::aes(y = mean),
                          shape = 21,size = 2,fill = "white") +
      ggplot2::geom_point(
        ggplot2::aes(y = observed), color = "purple",
        shape = 21,, size = 2.2) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = group,
        y = "Pr(response = 1)") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))
  }

  if (group == "person") {

    obsStat <- aggregate(resp ~ sbj, obs, mean)
    names(obsStat)[2] <- "observed"

    predStat <- aggregate(resp ~ draw + sbj,
                          data = yrep, mean)

    predSummary <-
      aggregate(resp ~ sbj, data = predStat,
                FUN = function(x) {
                  c(mean = mean(x),
                    lo50 = stats::quantile(x, 0.25),
                    hi50 = stats::quantile(x, 0.75),
                    lo95 = stats::quantile(x, 0.025),
                    hi95 = stats::quantile(x, 0.975))
                })

    predSummary <- data.frame(
      sbj = predSummary$sbj,
      predSummary$resp,row.names = NULL)

    plotDf <- merge(obsStat, predSummary, by = "sbj")

    p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = sbj)) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = `lo95.2.5.`, ymax = `hi95.97.5.`),
        linewidth = 0.4) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = `lo50.25.`, ymax = `hi50.75.`),
        linewidth = 1.1) +
      ggplot2::geom_point(ggplot2::aes(y = mean),
                          shape = 21,size = 1.2,fill = "white") +
      ggplot2::geom_point(
        ggplot2::aes(y = observed), color = "darkgreen",
        shape = 21,, size = 1.2) +
      #ggplot2::coord_flip() +
      ggplot2::labs(
        x = group,
        y = "Pr(response = 1)") +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))
  }

  stop("Unknown group option.", call. = FALSE)
}
