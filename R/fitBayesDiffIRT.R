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
#' @param chains `integer` the number of Markov chains to be run.
#' @param parallel_chains `integer` The maximum number of MCMC chains to run in parallel.
#' @param iter_warmup `integer` number of warmup (aka burnin) iterations.
#' @param iter_sampling `integer The number of post-warmup iterations to run per chain.
#' @return An object of class `BayesDiffIRTfit`. Available methods include `print`, `summary`, and `plot`.
#'
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
fitBayesDiffIRT <- function(data, rt = "rt", resp = "resp", sbj = "sbj",
                            item = "item", model = "d", priors = NULL,
                            seed = 42, chains = 3, parallel_chains = 1,
                            iter_warmup = 1000,
                            iter_sampling = 2000,
                            na.rm = TRUE){
  call <- match.call()

  # 0) Check stan

    # the fail should happen when the function is used,
    # not when the package is loaded.
  check_cmdstan()

  # 1) validate data
  data <- validate_data(data, rt, resp, sbj, item, na.rm = na.rm)

  # 2) Complete priors

  priors <- complete_priors(priors, model=model)

  # 3) transform the data into stan-friendly format & append priors

  stan_data <- make_stan_data(rt = data[,rt], resp = data[,resp],
                              sbj = data[,sbj], item = data[,item],
                              model = model, priors = priors)

  # 4) Compile + sample (cmdstanr)

  stanfile = switch(model,
                    "q" = "qdiffusion.stan",
                    "d" = "ddiffusion.stan",
                    stop("Error. Unknown model.\nPlease select one out of the followng: c('d', 'q')"), .call=FALSE)

  fullNameStanFile <-
    system.file("stan", stanfile,
                package = "BayesDiffIRT")
  if (fullNameStanFile == "") {
    stop("Stan file not found. Package installation may be corrupted.")
  }

  mod <- cmdstan_model(fullNameStanFile)
  fit <- mod$sample(data = stan_data, seed = seed, chains = chains,
                    parallel_chains = parallel_chains,
                    iter_warmup = iter_warmup,
                    iter_sampling = iter_sampling)

  # 5) Postprocess:
  # diagnostics <- collect_diagnostics(fit)
  # don't know if model diagnostics can are created automatically.

  # 6) Construct return object
  new_BayesDiffIRTfit(
    fit = fit,
    stan_data = stan_data,
    call = call#,
    #  diagnostics = diagnostics # divergences, rhats, etc.
  )
}

new_BayesDiffIRTfit <- function(fit, stan_data, model, call, diagnostics = NULL) {
  structure(
    list(
      fit = fit,
      stan_data = stan_data,
      model = model,
      call = call
    # , diagnostics = diagnostics - die müssen noch hinzugefügt werden.
    ),
    class = "BayesDiffIRTfit"
  )
}

print.BayesDiffIRTfit <- function(object){
  cat("Call:\n")
  print(object$call)
  summary(object)
}
