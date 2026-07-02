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
                            group = c("none", "item", "person"),
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

  if (group == "none") {

    obsStat <- data.frame(
      response = c(0L, 1L),
      observed = c(mean(obs$resp == 0L), mean(obs$resp == 1L))
    )

    predStat <- aggregate(
      resp ~ draw,
      data = yrep,
      FUN = function(x) mean(x == 1L)
    )

    predLong <- rbind(
      data.frame(
        response = 1L,
        draw = predStat$draw,
        value = predStat$resp
      ),
      data.frame(
        response = 0L,
        draw = predStat$draw,
        value = 1 - predStat$resp
      )
    )

    predSummary <- aggregate(
      value ~ response,
      data = predLong,
      FUN = function(x) {
        c(
          mean = mean(x),
          lo50 = stats::quantile(x, 0.25),
          hi50 = stats::quantile(x, 0.75),
          lo95 = stats::quantile(x, 0.025),
          hi95 = stats::quantile(x, 0.975)
        )
      }
    )

    predSummary <- data.frame(
      response = predSummary$response,
      do.call(rbind, predSummary$value),
      row.names = NULL
    )

    plotDf <- merge(obsStat, predSummary, by = "response")
    plotDf$response <- factor(plotDf$response, levels = c(0, 1))

    p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = response)) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo95, ymax = hi95),
        linewidth = 0.5
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo50, ymax = hi50),
        linewidth = 1.4
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = mean),
        shape = 21,
        size = 2.5,
        fill = "white"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = observed),
        size = 2.8
      ) +
      ggplot2::labs(
        x = "Response category",
        y = "Response proportion",
        subtitle =
          "Observed proportions and posterior predictive intervals"
      ) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))
  }

  if (group %in% c("item", "person")) {

    groupVar <- if (group == "item") "item" else "sbj"

    obsStat <- aggregate(
      resp ~ .,
      data = obs[, c(groupVar, "resp")],
      FUN = function(x) mean(x == 1L)
    )

    names(obsStat) <- c("group", "observed")

    predStat <- aggregate(
      resp ~ draw + group,
      data = transform(yrep, group = yrep[[groupVar]]),
      FUN = function(x) mean(x == 1L)
    )

    names(predStat)[names(predStat) == "resp"] <- "value"

    predSummary <- aggregate(
      value ~ group,
      data = predStat,
      FUN = function(x) {
        c(
          mean = mean(x),
          lo50 = stats::quantile(x, 0.25),
          hi50 = stats::quantile(x, 0.75),
          lo95 = stats::quantile(x, 0.025),
          hi95 = stats::quantile(x, 0.975)
        )
      }
    )

    predSummary <- data.frame(
      group = predSummary$group,
      do.call(rbind, predSummary$value),
      row.names = NULL
    )

    plotDf <- merge(obsStat, predSummary, by = "group")
    plotDf$group <- factor(plotDf$group, levels = plotDf$group)

    p <- ggplot2::ggplot(plotDf, ggplot2::aes(x = group)) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo95, ymax = hi95),
        linewidth = 0.4
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = lo50, ymax = hi50),
        linewidth = 1.1
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = mean),
        shape = 21,
        size = 2,
        fill = "white"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = observed),
        size = 2.2
      ) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = group,
        y = "Pr(response = 1)",
        subtitle = paste0(
          "Observed response proportions and posterior predictive intervals by ",
          group
        )
      ) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_classic()

    base::print(p)
    return(invisible(p))
  }

  stop("Unknown group option.", call. = FALSE)
}
