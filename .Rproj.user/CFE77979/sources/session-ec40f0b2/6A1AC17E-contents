#' @title Fit Bayesian diffusion item-response theory models
#'
#' @description
#' `fitBayesDiffIRT()` fits Bayesian diffusion item-response theory models by
#' sampling from the posterior distributions of item and subject parameters using
#' the No-U-Turn Sampler (NUTS) as implemented in Stan \insertCite{carpenter_stan_2017}{BayesDiffIRT} via \pkg{cmdstanr}.
#'
#' @param data A `data.frame` with one row per response. The data must contain
#' columns identifying the subject, item, response, and response time.
#'  * Response times should be numeric and measured in seconds.
#'  * For ability tests, binary responses should be coded as 0 for incorrect responses and 1 for correct responses.
#'  * For questionnaire items, binary responses should be coded as 0 for rejected items and 1 for accepted items.
#' @param rt Character string. Name of the column in `data` containing response
#' times. Defaults to "rt".
#' @param resp `character`. Name of the column in data containing
#' the binary response. Defaults to "resp".
#' @param sbj `character` Name of the column in data identifying subjects.
#' Defaults to "sbj".
#' @param item `character`. Name of the column in data identifying items.
#' Defaults to "item".
#' @param model `character`. Name of the diffusion item-response theory model to fit.
#' Currently implemented models:
#'  * "d" for the D-diffusion model (for survey items, default),
#'  * "dRV" for the D-diffusion model with random variability (for survey items),
#'  * "q" for the Q-diffusion model (for ability tests),
#'  * "qRV" for the Q-diffusion model with random variability (for ability tests).
#' @param priors Prior specification. Either NULL, a
#' BayesDiffIRTPrior object, or a list of BayesDiffIRTPrior objects created
#' with [BayesDiffIRT::prior()]. See [BayesDiffIRT::prior()] for details on specifying priors.
#' If NULL, a set of default priors is used,
#' @param seed Optional integer seed passed to Stan for reproducible sampling.
#' @param nChains Integer. Number of Markov chains. Defaults to 4.
#' @param nCores Integer. Number of chains to run in parallel. Defaults to
#' nChains.
#' @param nWarmup Integer. Number of warmup iterations per chain. Defaults to
#' 1000.
#' @param nSamples Integer. Number of post-warmup sampling iterations per chain.
#' Defaults to 1000.
#' @param adaptDelta Numeric. Target average proposal acceptance probability
#' passed to Stan. Higher values may reduce divergent transitions but can
#' increase computation time. Defaults to 0.95.
#' @param maxTreeDepth Integer. Maximum tree depth for the NUTS sampler.
#' Increasing this value may help when transitions hit the maximum tree depth,
#' but can increase computation time. Defaults to 12.
#' @param init Initial values passed to Stan. See
#' [cmdstanr::model-method-sample] for supported formats. Defaults to NULL.
#' @param refresh Integer. Number of iterations between progress messages printed
#' by Stan. Use 0 to suppress sampling progress output. Defaults to 200.
#' @param diagnostic.warnings Logical. If TRUE, warnings are issued when common
#' Stan diagnostics indicate potential sampling problems. Defaults to TRUE.
#' @param na.rm Logical. If TRUE, rows with missing values in required columns
#' are removed during data validation. If FALSE, missing values cause an
#' error. Defaults to TRUE.
#' @param ... Additional arguments passed to `cmdstanr::CmdStanModel$sample()`.

#' @return
#' An object of class BayesDiffIRTfit, which is a list containing:
#'
#' \describe{
#' \item{fit}{The fitted cmdstanr object.}
#' \item{stanData}{The data list passed to Stan.}
#' \item{model}{The fitted model name.}
#' \item{call}{The matched function call.}
#' \item{diag}{A list of Stan diagnostic summaries.}
#' }
#'
#' Methods for BayesDiffIRTfit objects include print(), summary(), plot(), and
#' checkDiagnostics().
#'
#' @author
#' Manuel Rausch, \email{manuel.rausch@@aau.at}
#'
#' @details
#' `fitBayesDiffIRT` samples from the posterior distributions of item and subject
#' parameters of responses and reaction times for
#' diffusion item response theory models
#' \insertCite{van_der_maas_cognitive_2011,tuerlinckx_two_2005,kang_modeling_2022}{BayesDiffIRT}.
#' ## Description of models
#' Diffusion item response theory combines item response theory with the drift diffusion
#' model of decision making. According to the drift diffusion decision model
#' \insertCite{Ratcliff2016}{BayesDiffIRT}, in each moment, the sensory system generates new evidence about which of the two possible choice options
#' is the correct one. This momentary evidence is drawn from a Gaussian
#' distribution and accumulated over time, i.e. the newly
#' acquired evidence is constantly added to the evidence collected up to that moment.
#' This accumulation process is bounded by an upper and lower threshold, where each threshold represents
#' one of the two possible choice options.  When the accumulated evidence reaches
#' one of the thresholds, a choice is made for the corresponding choice option.
#' The quality of information favouring one response option over the other is reflected
#' in the drift rate \eqn{\delta}, which quantifies s how quickly the accumulated evidence approaches
#' the threshold associated with the correct or preferred decision.
#' The distance between the two thresholds \eqn{\alpha} determines the amount of
#' evidence required before a decision is made; a larger distance means that decisions tend to
#' be made later because more evidence is required.
#' The starting point \eqn{\beta} of the accumulation process reflects a priori bias toward one of the response options.
#' In diffusion item response theory models,
#' two of the traditional parameters from the drift diffusion decision model,
#' boundary separation and drift rate, are decomposed into person and item parameters.
#' When person p makes a decision about item i, the boundary separation is given
#' by \eqn{\alpha_{pi} = \gamma_p / a_i}, where \eqn{\gamma_p} represents the
#' person-specific response caution and \eqn{a_i} item-specific time pressure.
#' The D-diffusion and the Q-diffusion model are characterised by the way
#' the drift rate is decomposed. According to the D-diffusion model
#' \insertCite{tuerlinckx_two_2005}{BayesDiffIRT}, which is applicable for survey items,
#' the drift rate is given by \eqn{\delta_{pi} = \theta_p - \nu_i}.
#' According to the Q-diffusion model \insertCite{van_der_maas_cognitive_2011}{BayesDiffIRT},
#' which is applicable for ability tests, the drift is given by
#' \eqn{\delta_{pi} = \theta_p / \nu_i}. In both the D-diffusion and the Q-diffusion model,
#' the accumulation always starts midway between the two response alternatives, i.e.,
#' there is no a prior bias for any of the choice alternatives.
#' \insertCite{kang_modeling_2022}{BayesDiffIRT} proposed to include random
#' trial-to-trial variability in both starting value \eqn{\beta} and drift rate \eqn{\delta}.
#' In the Q- and D-diffusion model with random variation, the starting point
#' \eqn{\beta_{pij}} for trial j given item i and person p is sampled from a uniform
#' distribution \eqn{\beta_{pij} \sim \mathcal{U}(0.5 - s_{\beta}/2, 0.5 + s_{\beta}/2)}.
#' The drift rate \eqn{\delta_{pij}} in trial j given item i and person p is sampled
#' from a Gaussian distribution \eqn{\delta_{pij} \sim \mathcal{N}(\delta_{ij}, s_{\delta}^2)}.
#' Random variability in starting values and drift rates accounts for the conditional
#' dependency of accuracy and reaction times \insertCite{kang_modeling_2022}{BayesDiffIRT},
#' but note that sampling will be considerably slower.

#' @examples
#' data("bayesDiffIRTexample")
#' \dontrun{
#' fit <- fitBayesDiffIRT(
#'   data = bayesDiffIRTexample,
#'   rt = "rt",
#'   response = "response",
#'   person = "person",
#'   item = "item",
#'   model = "d",
#'   nChains = 4,
#'   nWarmup = 1000,
#'   nSampling = 1000
#' )
#'
#' summary(fit)
#' plot(samples, parameter = "gamma", type = "interval",index=1:10)
#' plot(samples, parameter = "a", type = "interval",index=1:10)
#'
#' diagnostics(fit)
#'
#'
#' }

#' @seealso
#' \code{\link{prior}},
#' \code{\link{plot}},
#' \code{\link{checkDiagnostics}},
#' \code{\link[cmdstanr]{cmdstan_model}},
#' \code{\link[cmdstanr]{CmdStanModel}}
#'
#' @import cmdstanr
#' @importFrom Rdpack reprompt

#' @references
#' \insertAllCited{}

#' @export
fitBayesDiffIRT <- function(
    data, rt = "rt", resp = "resp", sbj = "sbj", item = "item",
    model = "d",
    priors = NULL,
    seed = NULL,
    nChains = 4,
    nCores = nChains,
    nWarmup = 1000,
    nSamples = 1000,
    adaptDelta = 0.95,
    maxTreeDepth = 12,
    init = NULL,
    refresh = 200,
    diagnostic.warnings = TRUE,
    na.rm = TRUE,
    ...
){
  call <- match.call()

  # 0) Check stan

  # the fail should happen when the function is used,
  # not when the package is loaded.
  checkCmdstan()

  # 1) validate data
  data <- validateData(data, rt, resp,
                        sbj, item, na.rm = na.rm)

  # 2) Complete priors

  priors <- completePriors(priors, model = model)

  # 3) transform the data into stan-friendly format & append priors

  stanData <-
    makeStanData(data, rt = rt, resp = resp,
                   sbj = sbj, item = item,
                   priors = priors)

  # 4) Compile + sample (cmdstanr)

  stanFile <- switch(
    model,
    "q" = "qdiffusion.stan",
    "d" = "ddiffusion.stan",
    "qRV" = "Q-DIRT-RV.stan",
    "dRV" =  "D-DIRT-RV.stan",
    stop("Error. Unknown model.\nPlease select one out of the followng: c('d', 'q')"), .call=FALSE)

  fullNameStanFile <-
    system.file("stan", stanFile,
                package = "BayesDiffIRT")
  if (fullNameStanFile == "") {
    stop("Stan file not found. Package installation may be corrupted.")
  }

  mod <- cmdstan_model(fullNameStanFile,
                       compile = TRUE, quiet = TRUE)
  fit <- mod$sample(data = stanData,
                    seed = seed,
                    chains = nChains,
                    parallel_chains = nCores,
                    iter_warmup = nWarmup,
                    iter_sampling = nSamples,
                    adapt_delta = adaptDelta,
                    max_treedepth = maxTreeDepth,
                    init = init,
                    refresh = refresh, ...)

  # 5) Postprocess:
  diag <- checkStanDiagnostics(fit, maxTreeDepth)
  if(diagnostic.warnings) warnStanDiagnostics(diag)

  # 6) Construct return object
  newBayesDiffIRTfit(
    fit = fit,
    stanData = stanData,
    model = model,
    call = call,
    diag = diag # divergences, rhats, etc.
  )
}

newBayesDiffIRTfit <- function(
    fit, stanData, model, call, diag) {
  structure(list(
    fit = fit,
    stanData = stanData,
    model = model,
    call = call,
    diag = diag),
    class = "BayesDiffIRTfit"
  )
}


warnStanDiagnostics <- function(x) {
  problems <- character()
  nDiv <- sum(x$divergences$nDivergent)

  if (nDiv > 0) {
    problems <-
      c(problems, paste0(nDiv, " divergent transition(s)"))
  }
  if (nrow(x$rhat) > 0) {
    problems <-
      c(problems, paste0(nrow(x$rhat), " parameter(s) with R-hat > 1.01"))
  }
  if (nrow(x$ess) > 0) {
    problems <-
      c(problems, paste0(nrow(x$ess), " parameter(s) with low ESS"))
  }
  if (any(x$ebfmi$warning, na.rm = TRUE)) {
    problems <-
      c(problems, "low E-BFMI in at least one chain")
  }
  if (sum(x$treedepth$n_max_treedepth) > 0) {
    problems <-
      c(problems, "at least one transition hit max treedepth")
  }
  if (length(problems) > 0) {
    warning(
      "Some Stan diagnostics indicate potential sampling problems:\n- ",
      paste(problems, collapse = "\n- "),
      "\nUse checkDiagnostics() for details.",
      call. = FALSE
    )
  }
  invisible(x)
}



