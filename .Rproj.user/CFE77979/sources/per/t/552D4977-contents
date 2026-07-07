#' @title
#' Plot response-probability surfaces

#' @description
#' `plotResponseSurface()` plots the posterior mean response probability for a
#' selected item over a two-dimensional grid of person parameters. The x-axis
#' represents the latent trait parameter theta, and the y-axis represents the
#' response caution parameter gamma. For each grid point, the response probability
#' is averaged over a random subset of posterior draws of the selected item's parameters.
#' The resulting surface is shown as a raster plot with contour lines.

#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A ggplot2 object.

#' @importFrom rtdists pdiffusion
#' @import ggplot2

#' @export
plotResponseSurface <- function(object, ...) {
  UseMethod("plotResponseSurface")
}

#' @rdname plotResponseSurface
#' @aliases plotResponseSurface.BayesDiffIRTfit
#' @method plotResponseSurface BayesDiffIRTfit

#' @param item Integer. Indexes the item for which the response-probability
#' surface should be plotted.
#' @param theta.range. Numeric vector of length two. Range of theta
#' values shown on the x-axis. If NULL, the range is set to the 1st and 99th
#' percentiles of the posterior draws of thetas form all items.
#' @param gamma.range. Numeric vector of length two. Range of gamma
#' values shown on the y-axis. If NULL, the range is set to the 1st and 99th
#' percentiles of the posterior draws of gamma from all items,
#' @param grid.size Integer. Number of grid points used for each axis. The total
#' number of evaluated grid points is grid.size * grid.size.
#' Defaults to 60.
#' @param contours Numeric vector. Response-probability values at which contour
#' lines should be drawn. Defaults to c(.1, .3, .5, .7, .9).
#' @param ndraws Integer. Number of posterior draws used to average the response
#' probability at each grid point. Defaults to 200.

#' @details
#' The method computes the predicted probability of a correct response (for ability
#' tests or of item acceptance (for survey questions) for a selected item across a grid of theta and gamma values.
#' The response probabilities are calculated using code{\link[rtdists:Diffusion]{rdiffusion}}\insertCite{Singmann2026}{BayesDiffIRT}..
#' @return a ggplot object

#' @seealso
#' \code{\link{fitBayesDiffIRT}},
#' \code{\link[rtdists:Diffusion]{rdiffusion}}

#' @importFrom rtdists rdiffusion
#' @importFrom posterior as_draws_df
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}

#' @examples
#' \dontrun{
#' data("bayesDiffIRTexample")
#' fit <- fitBayesDiffIRT(data=bayesDiffIRTexample, model = "d")
#'
#' plotResponseSurface(fit, item = 1)
#'
#' plotResponseSurface(
#' fit,
#' item = 2,
#' theta.range = c(-3, 3),
#' gamma.range = c(0.5, 3),
#' grid.size = 80,
#' contours = c(.25, .5, .75),
#' ndraws = 500
#' )
#' }
#
#' @author
#' Manuel Rausch, \email{manuel.rausch@@aau.at}

#' @export
plotResponseSurface.BayesDiffIRTfit <- function(
    object, item,
    theta.range = NULL,
    gamma.range = NULL,
    grid.size = 60,
    contours = c(.1, .3, .5, .7, .9),
    ndraws = 200) {

  if (!inherits(object, "BayesDiffIRTfit")) {
    stop("object must be a BayesDiffIRTfit object.", call. = FALSE)
  }

  if (length(item) != 1L || is.na(item) ||
      item < 1L || item > object$stanData$nItem) {
    stop("item must be a single valid item index.", call. = FALSE)
  }


  model <- object$model

  if(is.null(theta.range)){
    theta.range <-
      quantile(x = posterior::as_draws_matrix(object$fit$draws("theta")),
               probs = c(.01, .99))
  }
  if(is.null(gamma.range)){
    gamma.range <-
      quantile(x = posterior::as_draws_matrix(object$fit$draws("gamma")),
               probs = c(.01, .99))
  }

  # construct a grid of person parameters
  theta.grid <- seq(theta.range[1], theta.range[2],
                    length.out = grid.size)
  gamma.grid <- seq(gamma.range[1], gamma.range[2],
                    length.out = grid.size)

  surface <- expand.grid(
    theta = theta.grid,
    gamma = gamma.grid)

  # extract item parameter for each posterior draw

  pars <- c(paste0("nu[", item, "]"), paste0("a[", item, "]"))

  if (model %in% c("dRV", "qRV")) {
    pars <- c(pars, "s_delta", "s_beta")
  }

  drawsAll <- posterior::as_draws_df(object$fit$draws(pars))
  if (ndraws > nrow(drawsAll)){
    stop("ndraws should not be larger than the number of available posterior draws.")
  }
  drawsIdx <- sample(seq_len(nrow(drawsAll)), size = ndraws, replace = FALSE)
  draws <- as.data.frame(drawsAll[drawsIdx, , drop = FALSE])


  nu <- draws[[paste0("nu[", item, "]")]]
  a  <- draws[[paste0("a[", item, "]")]]

  if (model %in% c("dRV", "qRV")) {
    s_delta <- draws[["s_delta"]]
    s_beta  <- draws[["s_beta"]]
  } else {
    s_delta <- rep(0, ndraws)
    s_beta  <- rep(0, ndraws)
  }

  probSum <- numeric(nrow(surface))

  for (i in seq_len(ndraws)) {
    alpha <- surface$gamma / a[i]

    if (model %in% c("d", "dRV")) drift <- surface$theta - nu[i]
    if (model %in% c("q", "qRV")) drift <- surface$theta / nu[i]

    probSum <- probSum + rtdists::pdiffusion(
      rt = rep(Inf, nrow(surface)),
      response = rep("upper", nrow(surface)),
      a = alpha, v = drift, t0 = 0, z = alpha / 2,
      d = 0, sz = s_beta[i] * alpha,sv = s_delta[i])
  }
  surface$prob <- probSum / ndraws

  gg <- ggplot2::ggplot(
    surface,ggplot2::aes(x = theta,y = gamma)) +
    ggplot2::geom_raster(ggplot2::aes(fill = prob)) +
    ggplot2::geom_contour(
      ggplot2::aes(z = prob),
      breaks = contours,color = "white") +
    ggplot2::scale_fill_viridis_c(
      breaks = contours,
      limits = c(0, 1),
      name = "Pr(response = 1)")+
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      x = expression(person ~ ability ~ theta),
      y = expression(response ~ caution ~ gamma)) +
    ggplot2::theme_classic()

  base::print(gg)
  invisible(gg)

}
