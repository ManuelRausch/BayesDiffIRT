#' @title
#' Draw from the posterior predictive distribution

#' @description
#' Generate posterior predictive response and reaction-time data from a fitted
#' BayesDiffIRT model.
#' `posteriorPredict()` simulates replicated data sets from the posterior
#' predictive distribution of a fitted \code{BayesDiffIRTfit} object. For each
#' selected posterior draw, the function simulates one replicated response and
#' reaction time for each observed row in the original data set. The original
#' person-item observation structure is preserved. Thus, missing-by-design items,
#' repeated items, and unbalanced designs are retained in the posterior
#' predictive samples.
#' The function extracts posterior draws of the
#' model parameters from the fitted \code{BayesDiffIRTfit} object and uses these draws to simulate
#' replicated data with \code{\link[rtdists:Diffusion]{rdiffusion}}\insertCite{Singmann2026}{BayesDiffIRT}.

#' @param object A fitted model object of class \code{BayesDiffIRTfit}, typically
#'   returned by \code{\link{fitBayesDiffIRT}}.
#' @param ndraws Integer. Number of posterior draws to use for posterior
#'   predictive simulation. Defaults to \code{10}. The value must not exceed the
#'   number of available posterior draws in \code{object}.
#' @param seed Optional integer. Random seed used before drawing posterior
#'   samples and simulating replicated data. Defaults to \code{NULL}.

#' @return A \code{BayesDiffIRTpp} object, which is a  \code{data.frame} with \code{ndraws * object$stanData$nObs} rows.
#'   Each row corresponds to one replicated observation from one posterior draw.
#'   The returned data frame contains the following columns:
#'   \describe{
#'     \item{\code{draw}}{Index of the posterior predictive draw, from
#'       \code{1} to \code{ndraws}.}
#'     \item{\code{obs}}{Index of the original observation row.}
#'     \item{\code{sbj}}{Person index for the observation.}
#'     \item{\code{item}}{Item index for the observation.}
#'     \item{\code{rt}}{Simulated reaction time.}
#'     \item{\code{resp}}{Simulated binary response. Responses are coded as
#'       \code{1} for the upper boundary and \code{0} for the lower boundary.}
#'   }

#' @details
#' The replicated data use the same observation-level design as the original data
#' passed to the model. That is, for each posterior draw, the function simulates
#' one replicated response and reaction time for each original observation
#' \code{n = 1, ..., object$stanData$nObs}, using the corresponding person and
#' item indices stored in the \code{BayesDiffIRTfit} object.

#' @seealso
#' \code{\link{fitBayesDiffIRT}},
#' \code{\link[rtdists:Diffusion]{rdiffusion}}

#' @examples
#' data(bayesDiffIRTexample)
#' \dontrun{
#' fit <- fitBayesDiffIRT(data=bayesDiffIRTexample, model = "d")
#'
#' # Simulate 100 posterior predictive data sets
#' yrep <- posteriorPredict(fit, ndraws = 100, seed = 123)
#'
#' head(yrep)
#' }

#' @importFrom rtdists rdiffusion
#' @importFrom posterior as_draws_df
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}

#' @author
#' Manuel Rausch, \email{manuel.rausch@@aau.at}

#' @export
posteriorPredict <- function(object, ndraws = 10, seed = NULL) {

  if (!inherits(object, "BayesDiffIRTfit")) {
    stop("posteriorPredict requires a fitted BayesDiffIRTfit object.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # 1) extract info from BayesDiffIRTfit

  nSbjs <- object$stanData$nPerson
  nItems <-  object$stanData$nItem
  model <- object$model

  parameter <- c("theta", "gamma", "nu", "a", "tnd")

  if (model %in% c("dRV", "qRV")){
    parameter <- c(parameter, "s_delta", "s_beta")
  }

  drawsAll <- posterior::as_draws_df(object$fit$draws(parameter))

  if (ndraws > nrow(drawsAll)) {
    stop("ndraws cannot exceed the number of available posterior draws.")
  }

  draws <- sample(seq_len(nrow(drawsAll)), size = ndraws, replace = FALSE)
  myDraws <- as.data.frame(drawsAll[draws, , drop = FALSE])

  Out <- data.frame(
    draw = rep(seq_len(ndraws), each = object$stanData$nObs),
    obs  = rep(seq_len(object$stanData$nObs), times = ndraws),
    sbj  = rep(object$stanData$person, times = ndraws),
    item = rep(object$stanData$item, times = ndraws))

  thetas <- as.matrix(myDraws[paste0("theta[", seq_len(nSbjs), "]")])
  gammas <- as.matrix(myDraws[paste0("gamma[", seq_len(nSbjs), "]")])
  nus    <- as.matrix(myDraws[paste0("nu[", seq_len(nItems), "]")])
  as     <- as.matrix(myDraws[paste0("a[", seq_len(nItems), "]")])
  tnds   <- as.matrix(myDraws[paste0("tnd[", seq_len(nSbjs), "]")])

  idxSbj  <- cbind(Out$draw, Out$sbj)
  idxItem <- cbind(Out$draw, Out$item)

  gammas <- gammas[idxSbj]
  thetas <- thetas[idxSbj]
  tnds   <- tnds[idxSbj]
  nus    <- nus[idxItem]
  as     <- as[idxItem]

  # 2) Select a model specific simulaton function
  simFun <-  switch(
    model,
    "d" = simulateResponseD,
    "q"= simulateResponseQ,
    "dRV" = simulateResponseDRV,
    "qRV" = simulateResponseQRV,
    stop('Error! Model should be "d", "q", "dRV" or "qRV"'))

  # 3) Simulating data

  if (model %in% c("dRV", "qRV")){
    s_deltas <- myDraws[["s_delta"]][Out$draw]
    s_betas  <- myDraws[["s_beta"]][Out$draw]

    Sim <- mapply(
      simFun, gamma = gammas,
      theta = thetas, tnd = tnds,
      nu = nus, a = as,
      s_delta = s_deltas,
      s_beta = s_betas)
  } else {
    Sim <- mapply(
      simFun, gamma = gammas,
      theta = thetas, tnd = tnds,
      nu = nus, a = as)
  }

  Out$rt   <- as.numeric(Sim[1, ])
  Out$resp <- as.integer(Sim[2, ])
  Out
  class(Out) <- c(class(Out), "BayesDiffIRTpp")
  Out
}





simulateResponseD <- function(gamma, theta, tnd, nu, a){
  x <- rtdists::rdiffusion(n=1, a = gamma/a,
                           v = theta - nu,
                           t0 = tnd,
                           z = gamma/a * .5,
                           d=0, sz = 0, sv = 0, s = 1)
  resp <- as.character(x$response)
  resp <- as.integer(resp == "upper")
  c(x$rt, resp)
}

simulateResponseQ <- function(gamma, theta, tnd, nu, a){
  x <- rtdists::rdiffusion(n=1, a = gamma/a,
                           v = theta / nu,
                           t0 = tnd,
                           z = gamma/a * .5,
                           d=0, sz = 0, sv = 0, s = 1)
  resp <- as.character(x$response)
  resp <- as.integer(resp == "upper")
  c(x$rt, resp)
}

simulateResponseDRV <-
  function(gamma, theta, tnd, nu, a, s_delta, s_beta){
    x <- rtdists::rdiffusion(
      n=1, a = gamma/a, v = theta - nu,
      t0 = tnd, z = gamma/a * .5, # rtdists::rdiffusion is parametrized with the **absolute starting points!**
      d=0, sz = s_beta * gamma / a, # the stan module uses relative starting points
      sv = s_delta, s = 1)
    resp <- as.character(x$response)
    resp <- as.integer(resp == "upper")
    c(x$rt, resp)
  }

simulateResponseQRV <-
  function(gamma, theta, tnd, nu, a, s_delta, s_beta){
    x <- rtdists::rdiffusion(
      n=1, a = gamma/a, v = theta / nu,
      t0 = tnd, z = gamma/a * .5, # rtdists::rdiffusion is parametrized with the **absolute starting points!**
      d=0, sz = s_beta * gamma / a, # the stan module uses relative starting points
      sv = s_delta, s = 1)
    resp <- as.character(x$response)
    resp <- as.integer(resp == "upper")
    c(x$rt, resp)
  }
