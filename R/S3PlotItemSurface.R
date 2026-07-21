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

#' @param items Integer or integer vector. Indexes the items for which the response-probability
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
#' @importFrom posterior as_draws_df as_draws_matrix
#' @importFrom Rdpack reprompt
#' @importFrom utils txtProgressBar setTxtProgressBar
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
    object,
    items,
    theta.range = NULL,
    gamma.range = NULL,
    grid.size = 60,
    contours = c(.1, .3, .5, .7, .9),
    ndraws = 200,
    facet.ncol = NULL) {

  if (!inherits(object, "BayesDiffIRTfit")) {
    stop("object must be a BayesDiffIRTfit object.", call. = FALSE)
  }

  # Validate item indices
  if (!is.numeric(items) ||
      length(items) < 1L ||
      anyNA(items) ||
      any(!is.finite(items)) ||
      any(items != floor(items)) ||
      any(items < 1L | items > object$stanData$nItem)) {
    stop(
      "item must contain valid integer item indices.",
      call. = FALSE
    )
  }

  # Remove duplicate indices while retaining their order
  items <- unique(as.integer(items))

  if (length(grid.size) != 1L ||
      is.na(grid.size) ||
      grid.size < 2L ||
      grid.size != floor(grid.size)) {
    stop("grid.size must be a single integer greater than 1.",
         call. = FALSE)
  }

  if (length(ndraws) != 1L || is.na(ndraws) ||
      ndraws < 1L || ndraws != floor(ndraws)) {
    stop("ndraws must be a positive integer.", call. = FALSE)
  }

  if (!is.null(facet.ncol) &&
      (length(facet.ncol) != 1L ||
       is.na(facet.ncol) || facet.ncol < 1L ||
       facet.ncol != floor(facet.ncol))) {
    stop("facet.ncol must be NULL or a positive integer.",
         call. = FALSE)
  }

  model <- object$model

  if (!(model %in% c("d", "dRV", "q", "qRV"))) {
    stop("Unsupported model type.", call. = FALSE)
  }

  if (is.null(theta.range)) {
    theta.range <- stats::quantile(
      posterior::as_draws_matrix(object$fit$draws("theta")),
      probs = c(.01, .99))
  }

  if (is.null(gamma.range)) {
    gamma.range <- stats::quantile(
      posterior::as_draws_matrix(object$fit$draws("gamma")),
      probs = c(.01, .99))
  }

  checkRange <- function(x, name) {
    if (!is.numeric(x) ||length(x) != 2L ||
        anyNA(x) || any(!is.finite(x)) ||  x[1] >= x[2]) {
      stop(
        paste0(name, " must contain two finite, increasing values."),
        call. = FALSE)}
  }
  checkRange(theta.range, "theta.range")
  checkRange(gamma.range, "gamma.range")

  # Construct the common grid of person parameters
  theta.grid <- seq(theta.range[1],theta.range[2],length.out = grid.size)
  gamma.grid <- seq(gamma.range[1],gamma.range[2],length.out = grid.size)

  surface.grid <- expand.grid(
    theta = theta.grid,gamma = gamma.grid,KEEP.OUT.ATTRS = FALSE)

  n.grid <- nrow(surface.grid)

  # Extract item parameters for all requested items
  nu.pars <- paste0("nu[", items, "]")
  a.pars  <- paste0("a[", items, "]")

  pars <- c(nu.pars, a.pars)

  if (model %in% c("dRV", "qRV")) {
    pars <- c(pars, "s_delta", "s_beta")
  }

  draws.all <- posterior::as_draws_matrix(object$fit$draws(pars))

  if (ndraws > nrow(draws.all)) {
    stop("ndraws must not be larger than the number of available posterior draws.",
      call. = FALSE)
  }

  # Use the same posterior draws for every item
  draws.idx <-
    sample.int(n = nrow(draws.all),
               size = ndraws,
               replace = FALSE)

  draws <- draws.all[draws.idx, , drop = FALSE]

  if (model %in% c("dRV", "qRV")) {
    s_delta <- draws[, "s_delta"]
    s_beta  <- draws[, "s_beta"]
  } else {
    s_delta <- numeric(ndraws)
    s_beta  <- numeric(ndraws)
  }

  item.labels <- paste("Item", items)
  surface.list <- vector("list", length(items))

  # Values that do not change across draws or items
  rt.inf <- rep(Inf, n.grid)
  response.upper <- rep("upper", n.grid)

  n.steps <- ndraws * length(items)
  step <- 0L

  pb <- utils::txtProgressBar(min = 0, max = n.steps,style = 3)

  for (j in seq_along(items)) {

    current.item <- items[j]

    nu <- draws[, paste0("nu[", current.item, "]")]
    a  <- draws[, paste0("a[", current.item, "]")]

    prob.sum <- numeric(n.grid)

    for (i in seq_len(ndraws)) {

      alpha <- surface.grid$gamma / a[i]

      if (model %in% c("d", "dRV")) {
        drift <- surface.grid$theta - nu[i]
      } else {
        drift <- surface.grid$theta / nu[i]
      }

      prob.sum <- prob.sum + rtdists::pdiffusion(
        rt = rt.inf,response = response.upper,
        a = alpha, v = drift, t0 = 0,
        z = alpha / 2, d = 0,
        sz = s_beta[i] * alpha,sv = s_delta[i])

      step <- step + 1L
      utils::setTxtProgressBar(pb, step)
    }

    surface.item <- surface.grid
    surface.item$prob <- prob.sum / ndraws
    surface.item$item <- factor(paste("Item", current.item),
      levels = item.labels)
    surface.list[[j]] <- surface.item
  }

  close(pb)

  surface <- do.call(rbind, surface.list)
  rownames(surface) <- NULL

  gg <- ggplot2::ggplot(
    surface,
    ggplot2::aes(x = theta, y = gamma)) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = prob)) +
    ggplot2::geom_contour(
      ggplot2::aes(z = prob),
      breaks = contours,
      color = "white") +
    ggplot2::scale_fill_viridis_c(
      breaks = contours,
      limits = c(0, 1),
      name = "Pr(response = 1)") +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      x = expression(person ~ ability ~ theta),
      y = expression(response ~ caution ~ gamma)) +
    ggplot2::theme_classic()

  if (length(items) > 1L) {
    gg <- gg +
      ggplot2::facet_wrap(
        ggplot2::vars(item),
        ncol = facet.ncol)
  }

  base::print(gg)
  invisible(gg)
}
