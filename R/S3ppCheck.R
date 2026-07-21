#' Posterior predictive checks
#'
#' Perform posterior predictive checks for fitted Bayesian diffusion item
#' response models.
#'
#' Posterior predictive checks compare summaries of the observed data with the
#' same summaries computed from data simulated from the posterior predictive
#' distribution.
#'
#' @param object A fitted  `BayesDiffIRTfit` object
#' @param ... Additional arguments passed to methods.
#'
#' @return A plot object.
#'
#' @export
ppCheck <- function(object, ...) {
  UseMethod("ppCheck")
}

#' @rdname ppCheck
#' @aliases ppCheck.BayesDiffIRTfit

#' @param type Character string specifying the type of posterior predictive
#'   check. Currently supported options are \code{"response"} and
#'   \code{"rtQuantile"}. \code{"response"}  compares observed and
#'   posterior-predictive response proportions. \code{"rtQuantile"}
#'   compares observed and posterior-predictive reaction-time quantiles
#'   separately by response category.
#' @param group Character string specifying the grouping level.
#'   Options are \code{"none"}, \code{"item"}, and \code{"person"}. If
#'   \code{"none"}, the check is computed for the full data set. If
#'   \code{"item"}, summaries are computed separately by item. If
#'   \code{"person"}, summaries are computed separately by person.
#' @param probs Numeric vector of probabilities used for reaction-time quantile
#'   checks. Values must lie between 0 and 1. Only used if type is \code{"rtQuantile"}.
#' @param minN Integer. Minimum number of observations required to compute an
#'   observed or posterior-predictive reaction-time quantile within a response
#'   category or group. Only used if type is \code{"rtQuantile"}.
#' @param yrep Optional BayesDiffIRTpp object created with  \code{\link{posteriorPredict}},
#' @param ndraws Optional Integer. Number of posterior predictive draws to generate. Only used
#'   if yrep is NULL.
#' @param seed Optional integer seed passed to \code{\link{posteriorPredict}}.
#' @param index Optional integer vector selecting specific items or persons when
#'   \code{group = "item"} or \code{group = "person"}. Ignored when
#'   \code{group = "none"}.
#'
#' @details
#' For \code{type = "response"}, the function compares the observed proportion
#' of \code{resp = 1} with the posterior predictive distribution of the same
#' quantity. With \code{group = "none"}, the plot shows the posterior predictive
#' distribution of the overall response proportion and marks the observed value.
#' With \code{group = "item"} or \code{group = "person"}, the plot shows
#' posterior predictive 50 percent and 95 percent intervals for each item or
#' person together with the observed response proportion.
#'
#' For \code{type = "rtQuantile"}, the function compares observed
#' reaction-time quantiles with posterior predictive intervals for the same
#' quantiles. Quantiles are computed separately for \code{resp = 0} and
#' \code{resp = 1}. With \code{group = "item"} or \code{group = "person"}, the checks are
#' computed separately for each item or person.

#' @return returns a \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{posteriorPredict}}, \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#' \dontrun{
#' # Overall response proportion check
#' ppCheck(fit, type = "response")
#'
#' # Response proportion check by item
#' ppCheck(fit, type = "response", group = "item")
#'
#' # Restrict item-level checks to selected items
#' ppCheck(fit, type = "response", group = "item", index = 1:5)
#'
#' # Reaction-time quantile check by response category
#' ppCheck(fit, type = "rtQuantile", probs = c(.1, .3, .5, .7, .9))
#'
#' # Reaction-time quantile check by item
#' ppCheck(fit, type = "rtQuantile", group = "item", index = 1:4)
#' }
#'

#' @method ppCheck BayesDiffIRTfit
#' @export


ppCheck.BayesDiffIRTfit <-
  function(object,
           type = c( "response", "rtQuantile"),
           group = c("none", "item", "person"),
           probs = c(.1, .3, .5, .7, .9),
           minN = 10,
           yrep = NULL,
           ndraws = 10,
           seed = NULL,
           index = NULL,
           ...) {
    type <- match.arg(type)
    group <- match.arg(group)

    if (!inherits(object, "BayesDiffIRTfit")) {
      stop("`object` must inherit from class \"BayesDiffIRTfit\".",
           call. = FALSE)
    }

    if (is.null(yrep)){
      yrep <- posteriorPredict(object, ndraws = ndraws, seed = seed)
    }

    if (!inherits(yrep, "BayesDiffIRTpp")) {
      stop("yrep must be a BayesDiffIRTpp object, created with posteriorPredict",
           call. = FALSE)
    }

    if (type == "response") {
      return(ppCheckResponse(object, yrep = yrep,group = group, seed = seed, index=index))
    }
    if (type ==  "rtQuantile"){
      return(ppCheckRtQuantile(
        object = object, yrep = yrep,group = group,
        probs = probs,minN = minN,seed = seed, index=index))
    }
    stop("PPC type not implemented.", call. = FALSE)
  }


ppCheckResponse <- function(object,
                            yrep = yrep,
                            group = c("none", "item", "person"),
                            index = NULL,
                            seed = NULL) {

  group <- match.arg(group)

  #yrep <- posteriorPredict(object, ndraws = ndraws, seed = seed)

  if (is.null(object$stanData$resp)) {
    stop("object$stanData$resp is required for response PPCs.", call. = FALSE)
  }

  obs <- data.frame(
    obs  = seq_len(object$stanData$nObs),
    sbj  = object$stanData$person,
    item = object$stanData$item,
    resp = object$stanData$resp
  )

  if (group == "none"){
    obsStat <- mean(object$stanData$resp)
    predStat <- stats::aggregate(stats::as.formula("resp ~ draw"),
                                 data = yrep, mean)

    p <- ggplot2::ggplot(predStat, ggplot2::aes(x=resp)) +
      ggplot2::geom_histogram(binwidth=.01,
                              colour="black",
                              fill="seagreen") +
      ggplot2::labs(x = "Pr(response = 1)",y = "Number of samples") +
      ggplot2::geom_vline(ggplot2::aes(xintercept=obsStat), size=2) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))
  }

  if (group == "item") {

    if (!is.null(index)) {
      obs <- obs[obs$item %in% index, , drop = FALSE]
      yrep <- yrep[yrep$item %in% index, , drop = FALSE]
    }

    obsStat <- stats::aggregate(
      obs["resp"],by = list(item = obs$item),FUN = mean)
    names(obsStat)[2] <- "observed"

    predStat <- stats::aggregate(
      yrep["resp"],by = list(draw = yrep$draw, item = yrep$item), FUN = mean)

    predSummary <- stats::aggregate(
      stats::as.formula("resp ~ item"),
      data = predStat,
      FUN = function(x) {
        c(mean = mean(x),
          lo50 = unname(stats::quantile(x, 0.25)),
          hi50 = unname(stats::quantile(x, 0.75)),
          lo95 = unname(stats::quantile(x, 0.025)),
          hi95 = unname(stats::quantile(x, 0.975))
        )})

    predSummary <- data.frame(
      item = predSummary$item,
      predSummary$resp,
      row.names = NULL)

    plotDf <- merge(obsStat, predSummary, by = "item")

    plotDf$item <- factor(plotDf$item,
                          levels = sort(unique(plotDf$item)))

    p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = .data[["item"]])) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data[["lo95"]], ymax = .data[["hi95"]]),
        linewidth = 0.4
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data[["lo50"]], ymax = .data[["hi50"]]),
        linewidth = 1.1
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data[["mean"]]),
        shape = 21,
        size = 2,
        fill = "green4"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data[["observed"]]),
        size = 2.8
      ) +
      ggplot2::coord_flip(ylim = c(0, 1)) +
      ggplot2::labs(
        x = "Item",
        y = "Pr(response = 1)"
      ) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))
  }

  if (group == "person") {

    if (!is.null(index)) {
      obs <- obs[obs$sbj %in% index, , drop = FALSE]
      yrep <- yrep[yrep$sbj %in% index, , drop = FALSE]
    }

    obsStat <- stats::aggregate(obs["resp"],
                                by = list(sbj = obs$sbj),FUN = mean)
    names(obsStat)[2] <- "observed"

    predStat <- stats::aggregate(
      yrep["resp"],by = list(draw = yrep$draw, sbj = yrep$sbj),FUN = mean)

    predSummary <- stats::aggregate(
      stats::as.formula("resp ~ sbj"),
      data = predStat,
      FUN = function(x) {
        c(
          mean = mean(x),
          lo50 = unname(stats::quantile(x, 0.25)),
          hi50 = unname(stats::quantile(x, 0.75)),
          lo95 = unname(stats::quantile(x, 0.025)),
          hi95 = unname(stats::quantile(x, 0.975))
        )
      }
    )

    predSummary <- data.frame(
      sbj = predSummary$sbj,
      predSummary$resp,
      row.names = NULL
    )

    plotDf <- merge(obsStat, predSummary, by = "sbj")

    plotDf$sbj <- factor(plotDf$sbj, levels = sort(unique(plotDf$sbj)))

    p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = .data[["sbj"]])) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data[["lo95"]], ymax = .data[["hi95"]]),
        linewidth = 0.4
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data[["lo50"]], ymax = .data[["hi50"]]),
        linewidth = 1.1
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data[["mean"]]),
        shape = 21,
        size = 1.2,
        fill = "white"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data[["observed"]]),
        fill = "lightgoldenrod",
        shape = 21,
        size = 2
      ) +
      ggplot2::labs(
        x = "Person",
        y = "Pr(response = 1)"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))

  }

  stop("Unknown group option.", call. = FALSE)

}

ppCheckRtQuantile <- function(object,
                              yrep = yrep,
                              group = "none",
                              probs = c(.1, .3, .5, .7, .9),
                              minN = 5,
                              seed = NULL,
                              index=NULL) {

  if (is.null(object$stanData$rt)) {
    stop("object$stanData$rt is required for RT quantile PPCs.", call. = FALSE)
  }

  if (is.null(object$stanData$resp)) {
    stop("object$stanData$resp is required for RT quantile PPCs.", call. = FALSE)
  }

  if (is.null(object$stanData$nObs)) {
    stop("object$stanData$nObs is required for RT quantile PPCs.", call. = FALSE)
  }

  if (length(object$stanData$rt) != object$stanData$nObs ||
      length(object$stanData$resp) != object$stanData$nObs) {
    stop("rt and resp must both have length object$stanData$nObs.", call. = FALSE)
  }

  if (!all(probs > 0 & probs < 1)) {
    stop("probs must contain values between 0 and 1.", call. = FALSE)
  }

  #yrep <- posteriorPredict(object, ndraws = ndraws, seed = seed)

  obs <- data.frame(
    obs  = seq_len(object$stanData$nObs),
    sbj  = object$stanData$person,
    item = object$stanData$item,
    rt   = object$stanData$rt,
    resp = object$stanData$resp
  )

  obs <- obs[!is.na(obs$rt) & !is.na(obs$resp), , drop = FALSE]
  yrep <- yrep[!is.na(yrep$rt) & !is.na(yrep$resp), , drop = FALSE]

  if (group == "none") {

    obsStat <- computeObservedRtQuantiles(
      data = obs,probs = probs, minN = minN)

    predStat <- computePredictedRtQuantiles(
      data = yrep,probs = probs,minN = minN)

    plotDf <- merge(obsStat,predStat,
                    by = c("resp", "prob"),all = FALSE)

    plotDf$resp <- factor(plotDf$resp)

    p <- ggplot2::ggplot(
      plotDf, ggplot2::aes(x = .data[["prob"]])) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["lo95"]], ymax = .data[["hi95"]]),
        alpha = 0.20, fill= "green4") +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["lo50"]], ymax = .data[["hi50"]]),
        alpha = 0.35, fill="green4") +
      ggplot2::geom_line(
        ggplot2::aes(y = .data[["mean"]])) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data[["observed"]]),
        size = 2) +
      ggplot2::facet_wrap(
        stats::as.formula("~resp"),
        labeller = ggplot2::labeller(
          resp = c(`0` = "resp = 0", `1` = "resp = 1"))) +
      ggplot2::scale_x_continuous(breaks = probs, limits = c(0,1)) +
      ggplot2::labs(
        x = "RT quantile",
        y = "Reaction time [s]") +
      ggplot2::theme_classic()


    base::print(p)
    return(invisible(p))
  }

  if (group %in% c("item", "person")) {

    groupVar <- if (group == "item") "item" else "sbj"
    if (!is.null(index)) {
      obs <- obs[obs[,groupVar] %in% index, , drop = FALSE]
      yrep <- yrep[yrep[,groupVar] %in% index, , drop = FALSE]
    }

    obs2 <- obs
    obs2$group <- obs2[[groupVar]]

    yrep2 <- yrep
    yrep2$group <- yrep2[[groupVar]]

    obsStat <- computeObservedRtQuantiles(
      data = obs2, probs = probs, minN = minN, group = TRUE)

    predStat <- computePredictedRtQuantiles(
      data = yrep2, probs = probs, minN = minN, group = TRUE )

    plotDf <- merge(obsStat, predStat,
                    by = c("group", "resp", "prob"), all = FALSE)

    plotDf$resp <- factor(plotDf$resp)
    plotDf$group <- factor(plotDf$group)

    p <- ggplot2::ggplot(
      plotDf, ggplot2::aes(x = .data[["prob"]])) +
      ggplot2::facet_grid(
        stats::as.formula("resp ~ group"),
        labeller = ggplot2::labeller(
          resp = c(`0` = "resp = 0", `1` = "resp = 1"),
          group = function(x) paste0(groupVar, " ", x))) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["lo95"]], ymax = .data[["hi95"]]),
        alpha = 0.20, fill= "green4") +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[["lo50"]], ymax = .data[["hi50"]]),
        alpha = 0.35, fill="green4") +
      ggplot2::geom_line(
        ggplot2::aes(y = .data[["mean"]])) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data[["observed"]]),
        size = 2) +
      ggplot2::scale_x_continuous(breaks = probs, limits = c(0,1)) +
      ggplot2::labs(
        x = "RT quantile",
        y = "Reaction time [s]") +
      ggplot2::theme_classic()



    base::print(p)
    return(invisible(p))
  }

  stop("Unknown group option.", call. = FALSE)
}

computeObservedRtQuantiles <- function(data,
                                       probs,
                                       minN = 5,
                                       group = FALSE) {

  splitVars <- if (isTRUE(group)) {
    c("group", "resp")
  } else {
    "resp"
  }

  splitFac <- interaction(data[splitVars], drop = TRUE, sep = "___")
  splitData <- split(data, splitFac)

  out <- lapply(splitData, function(d) {

    if (nrow(d) < minN) {
      return(NULL)
    }

    qs <- stats::quantile(d$rt, probs = probs, na.rm = TRUE, names = FALSE)

    res <- data.frame(
      resp = d$resp[1],
      prob = probs,
      observed = as.numeric(qs))

    if (isTRUE(group)) {
      res$group <- d$group[1]
      res <- res[, c("group", "resp", "prob", "observed")]
    }

    res
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

computePredictedRtQuantiles <- function(data,
                                        probs,
                                        minN = 5,
                                        group = FALSE) {

  splitVars <- if (isTRUE(group)) {
    c("draw", "group", "resp")
  } else {
    c("draw", "resp")
  }

  splitFac <- interaction(data[splitVars], drop = TRUE, sep = "___")
  splitData <- split(data, splitFac)

  perDraw <- lapply(splitData, function(d) {

    if (nrow(d) < minN) {
      return(NULL)
    }

    qs <- stats::quantile(d$rt, probs = probs, na.rm = TRUE, names = FALSE)

    res <- data.frame(
      draw = d$draw[1],
      resp = d$resp[1],
      prob = probs,
      value = as.numeric(qs))

    if (isTRUE(group)) {
      res$group <- d$group[1]
      res <- res[, c("draw", "group", "resp", "prob", "value")]
    }

    res
  })

  perDraw <- do.call(rbind, perDraw)
  rownames(perDraw) <- NULL

  summaryVars <- if (isTRUE(group)) {
    c("group", "resp", "prob")
  } else {
    c("resp", "prob")
  }
  splitFac <- interaction(perDraw[summaryVars], drop = TRUE, sep = "___")
  splitSummary <- split(perDraw, splitFac)

  out <- lapply(splitSummary, function(d) {
    vals <- d$value
    res <- data.frame(
      resp = d$resp[1],
      prob = d$prob[1],
      mean = mean(vals, na.rm = TRUE),
      lo50 = unname(stats::quantile(vals, 0.25, na.rm = TRUE)),
      hi50 = unname(stats::quantile(vals, 0.75, na.rm = TRUE)),
      lo95 = unname(stats::quantile(vals, 0.025, na.rm = TRUE)),
      hi95 = unname(stats::quantile(vals, 0.975, na.rm = TRUE))
    )
    if (isTRUE(group)) {
      res$group <- d$group[1]
      res <- res[, c("group", "resp", "prob",
                     "mean", "lo50", "hi50", "lo95", "hi95")]
    }
    res
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}
