#' @title Fit Bayesian diffusion item-response theory models
#'
#' @description
#' `fitBayesDiffIRT()` fits Bayesian diffusion item-response theory models by
#' sampling from the posterior distributions of item and subject parameters using
#' the No-U-Turn Sampler (NUTS) as implemented in Stan via `cmdstanr`.
#'
#' @param data A `data.frame` with one row per response. The data must contain
#' columns identifying the subject, item, response, and response time.
#'  * Response times should be numeric and measured in seconds.
#'  * For ability tests, binary responses should be coded as 0 for incorrect responses and 1 for correct responses.
#'  * For questionnaire items, binary responses should be coded as 0 for rejected items and 1 for accepted items.
#' @param rt Character string. Name of the column in `data` containing response
#' times. Defaults to "rt".
#' @param resp `character`. Name of the column in data containing
#' the binary response. Defaults to `resp`.
#' @param sbj `character` Name of the column in data identifying subjects.
#' Defaults to "sbj".
#' @param item `character`. Name of the column in data identifying items.
#' Defaults to "item".
#' @param model `character` . Name of the diffusion item-response theory model to fit.
#' Currently implemented models are "d" for the D-diffusion model (for survey items) and
#' "q" for the Q-diffusion model (for ability tests). Defaults to "d".
#' @param priors Prior specification. Either NULL, a
#' BayesDiffIRTPrior object, or a list of BayesDiffIRTPrior objects created
#' with [BayesDiffIRT::prior()]. ee [BayesDiffIRT::prior()] for details on specifying priors.
#' If NULL, a set of default priors is used,
#' @param seed Optional integer seed passed to Stan for reproducible sampling.
#' @param n.chains Integer. Number of Markov chains. Defaults to 4.
#' @param n.cores Integer. Number of chains to run in parallel. Defaults to
#' n.chains.
#' @param n.warmup Integer. Number of warmup iterations per chain. Defaults to
#' 1000.
#' @param n.samples Integer. Number of post-warmup sampling iterations per chain.
#' Defaults to 1000.
#' @param adapt.delta Numeric. Target average proposal acceptance probability
#' passed to Stan. Higher values may reduce divergent transitions but can
#' increase computation time. Defaults to 0.95.
#' @param max.treedepth Integer. Maximum tree depth for the NUTS sampler.
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
#' Methods for BayesDiffIRTfit objects include print(), summary(),
#' plot(), and checkDiagnostics().
#'
#' @author
#' Manuel Rausch, \email{manuel.rausch@@aau.at}
#'
#' @details
#' `BayesDiffIRT` samples from the posterior distributions of item and subject
#' parameters of responses and reaction times for
#' diffusion item response theory models \insertCite{van_der_maas_cognitive_2011,molenaar_fitting_2015,tuerlinckx_two_2005}{BayesDiffIRT}.
#' ## Description of models
#' Diffusion item response theory combines item response theory with the drift diffusion
#' model of decision making.
#' According to the drift diffusion decision model, in each moment, the sensory
#' system generates new evidence about which of the two possible choice options
#' is the correct one. This momentary evidence is drawn from a Gaussian
#' distribution and accumulated over time, i.e. the newly
#' acquired evidence is constantly added to the evidence collected up to that moment.
#' This accumulation process is bounded by an upper and lower threshold, where each threshold represents
#' one of the two possible choice options.  When the accumulated evidence reaches
#' one of the thresholds, a choice is made for the corresponding choice option.
#' The quality of information favouring one response option over the other is reflected
#' in the drift rate  \eqn{delta}, which quantifies s how quickly the accumulated evidence approaches
#' the threshold associated with the correct or preferred decision.
#' The distance between the two thresholds \eqn{alpha} determines the amount of evidence required
#' before a decision is made; a larger distance means that decisions tend to
#' be made later because more evidence is required.
#' The starting point \eqn{beta} of the accumulation process could be used to describe an a priori bias toward one of the response options.
#' However, in the current set of implemented that the accumulation always starts midway between the two response alternatives, i.e.,
#' there is no a prior bias for any of the choice alternatives .
#' In diffusion item response theory models, two of the traditional parameters
#' from the drift diffusion decision model, boundary separation and drift rate,
#'  are decomposed into person and item parameters. When person p makes a decision
#'  about item i, the boundary separation is given by \eqn{alpha_pi = gamma_p/a_i}, where
#'  \eqn{gamma_p} represents the person-specific response caution and \eqn{a_i}
#'  item-specific time pressure. The D-Diffusion and the Q-Diffusion model differ
#'   in the way the drift decomposed. According to the D-Diffusion model, which is applicable for
#'   survey item, the drift is given by \eqn{delta_ji = theta_j - nu_i}. According
#'   to the Q-Diffusion model, the drift is given by \eqn{delta_ji = theta_j / nu_i}.
#'

#' @examples
#' print("Coming!")

#' @import cmdstanr
#' @importFrom Rdpack reprompt

#' @references \insertAllCited{}

#' @export
fitBayesDiffIRT <- function(
    data, rt = "rt", resp = "resp", sbj = "sbj", item = "item",
    model = "d",
    priors = NULL,
    seed = NULL,
    n.chains = 4,
    n.cores = n.chains,
    n.warmup = 1000,
    n.samples = 1000,
    adapt.delta = 0.95,
    max.treedepth = 12,
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
                    chains = n.chains,
                    parallel_chains = n.cores,
                    iter_warmup = n.warmup,
                    iter_sampling = n.samples,
                    adapt_delta = adapt.delta,
                    max_treedepth = max.treedepth,
                    init = init,
                    refresh = refresh, ...)

  # 5) Postprocess:
  diag <- checkStanDiagnostics(fit, max.treedepth)
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



