#' @title Fit Bayesian diffusion item-response theory models
#'
#' @description `fitBayesDiffIRT` samples item and subject parameters
#' from the posterior distributions for diffusion item response theory
#' models. For this purpose, NUTS sampling as implemented in STAN is used. See Details for the mathematical specification of the implemented models and
#' their parameters.
#'
#' @param data  a `data.frame` with different row indicating responses to different items.
#' It must include at least three columns:
#'  * One column must provide response times. They should be numeric and measured in seconds.
#'  * One column must provide integer responses.
#'    - For ability tests, the value should be 0 for incorrect responses
#'    or correct responses.
#'    - For questionnaires, the value should be 0 for rejected items and 1 for accepted items.
#' @param rt `character` the name of the column in data that provides reaction times.
#' @param response `character` the name of the column in data that provides reaction times.
#' @param id `character` the name of the column in data in which responses are found.
#' @param model `character` the name of the model of which the parameters should be sampled. So far,
#' The following models have been implemented:  'D', 'Q'
#' @param priors BayesDiffIRTPrior object, or a list of BayesDiffIRTPrior objects, created with prior().
#' @param seed A seed for the random number generators to be passed to Stan.
#' @param n.chains Number of Markov chains.
#' @param n.cores Number of chains to run in parallel.
#' @param n.warmup Number of warmup iterations per chain.
#' @param n.samples Number of post-warmup sampling iterations per chain.
#' @param ... Additional arguments passed to `cmdstanr::CmdStanModel$sample()`.
#'   This is intended for advanced users. Common sampling controls such as
#'   `chains`, `parallel_chains`, `iter_warmup`, `iter_sampling`,
#'   `adapt_delta`, `max_treedepth`, and `init` should be supplied via the
#'   corresponding arguments of `fitBayesDiffIRT()`.


#' @return An object of class `BayesDiffIRTfit`. Available methods include `print`, `summary`, and `plot`.

#' @author
#' Manuel Rausch, \email{manuel.rausch@aau.at}
#'
#' @details
#' `BayesDiffIRT` samples from the posterior distributions of item and subject parameters
#' and from posterior predictive distributions of responses and reaction times for
#' diffusion item response theory models
#' \insertCite{van_der_maas_cognitive_2011,molenaar_fitting_2015,tuerlinckx_two_2005}{BayesDiffIRT}.
#' ## Mathematical description of models
#' According to the drift diffusion decision model, in each moment, the sensory
#' system generates new sensory evidence about which of the two possible choice options
#' is the correct one. This momentary evidence is drawn from Gaussian
#' distributions and accumulated over time, i.e. the newly
#' acquired evidence is constantly added to the evidence collected up
#'  to that moment. This accumulation process is bounded by an upper

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



