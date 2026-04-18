#' @title Fit Bayesian diffusion item-response theory models
#'
#' @description `fitBayesDiffIR` samples item and subject parameters
#' from the posterior distributions for diffusion item response theory
#' models. For this purpose, NUTS sampling is used. See Details for the mathematical specification of the implemented models and
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
#' @param id `character` the name of the column in data in which responses are foundedn.
#' @param model `character` the name of the model of which the parameters should be sampled. So far,
#' the following models have been implemented:  Models implemented so far: 'D', 'Q'
#' @param prior BayesDiffIRTPrior object, or a list of BayesDiffIRTPrior objects, created with prior().
#' @param seed A seed for the random number generators to pass to Stan.
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
#' diffusion item response theory models (REFERENCE MISSING).
#' ## Mathematical description of models
#' According to the drift diffusion decision model, in each moment, the sensory
#' system generates new sensory evidence about which of the two possible choice options
#' is the correct one. This momentary evidence is drawn from Gaussian
#' distributions and accumulated over time, i.e. the newly
#' acquired evidence is constantly added to the evidence collected up
#'  to that moment. This accumulation process is bounded by an upper
#'


#'
#' @examples
#' # 1) Prepare dataset
#' data(extraversion, package="diffIRT")
#' x=extraversion[,1:10]
#' x <- cbind(1:nrow(x), x)
#' colnames(x) <- c("sbj", paste0("Item", 1:10))
#' x <- tidyr::pivot_longer(as.data.frame(x), cols=Item1:Item10,
#'                          names_to="item",
#'                          values_to = "resp")
#'
#' rt=extraversion[,11:20]
#' rt <-  cbind(1:nrow(rt), rt)
#' colnames(rt) <- c("sbj", paste0("Item", 1:10))
#' rt <- tidyr::pivot_longer(as.data.frame(rt), cols=Item1:Item10,
#'                           names_to="item",
#'                           values_to = "rt")
#' Extra <- merge(x, rt)
#'
#' # 2) sample from the posterior
#' samples <-
#'   fitBayesDiffIRT(Extra,
#'                   rt = "rt", resp = "resp", sbj = "sbj",
#'                   item = "item", model = "d",
#'                   prior = NULL,
#'                   seed = 42, chains = 1,
#'                   parallel_chains = 1,
#'                   iter_warmup =  1000,
#'                   iter_sampling = 1000)
#' summary(samples)

#' @import cmdstanr

#' @references Alexandrowicz, R. W. (2020). The diffusion model visualizer: An interactive tool to understand the diffusion model parameters. Psychological Research, 84(4), 1157-1165. https://doi.org/10.1007/s00426-018-1112-6
#' @references Molenaar, D., Tuerlinckx, F., & Maas, H. L. J. V. D. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. Journal of Statistical Software, 66(4). https://doi.org/10.18637/jss.v066.i04
#' @references Van Der Maas, H. L. J., Molenaar, D., Maris, G., Kievit, R. A., & Borsboom, D. (2011). Cognitive psychology meets psychometric theory: On the relation between process models for decision making and latent variable models for individual differences. Psychological Review, 118(2), 339-356. https://doi.org/10.1037/a0022749


#' @export
fitBayesDiffIRT <- function(data, rt = "rt", resp = "resp", sbj = "sbj",
                            item = "item", model = "d", prior = NULL,
                            seed = 42, chains = 3, parallel_chains = 1,
                            iter_warmup =  1000,
                            iter_sampling = 2000){
  call <- match.call()

  # 0) Check stan
  # According to CRAN guidelines, the fail should happen when the function is used, not when the package is loaded.
  check_cmdstan()

  # 1) Preprocessing.

  # Complete priors

  priors <- complete_priors(priors, model=model)

  if (!rt %in% colnames(data)) stop(paste0("No column ", rt, "found in data" ),call. = FALSE)
  if (!resp %in% colnames(data)) stop(paste0("No column ", resp, "found in data" ),call. = FALSE)
  if (!sbj %in% colnames(data)) stop(paste0("No column ", sbj, "found in data" ),call. = FALSE)
  if (!item %in% colnames(data)) stop(paste0("No column ", item, "found in data" ),call. = FALSE)

  # transform the data into stan-friedly format & append priors

  stan_data <- make_stan_data(rt = data[,rt], resp = data[,resp],
                              sbj = data[,sbj], item = data[,item],
                              model = model, priors = priors)

  # 2) Compile + sample (cmdstanr)
  stanfile = switch(model, "q" = "qdiffusion.stan",
                    "d" = "ddiffusion.stan",
                    stop("Error. Unknown model. Please select one out of the followng: c('d', 'q')"), .call=FALSE)
  fullNameStanFile <-
    system.file("stan", stanfile,
                package = "BayesDiffIRT")
  if (fullNameStanFile == "") {
    stop("Stan file not found. Package installation may be corrupted.")
  }

  mod <- cmdstan_model(stanfile)
  fit <- mod$sample(data = stan_data, seed = seed, chains = chains,
                    parallel_chains = parallel_chains,
                    iter_warmup = iter_warmup,
                    iter_sampling = iter_sampling)

  # 3) Postprocess
  # diagnostics <- collect_diagnostics(fit)
  # don't know if model diagnostics can are created automatically.

  # 4) Construct return object
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
      call = call,
      diagnostics = diagnostics
    ),
    class = "BayesDiffIRTfit"
  )
}

#' keyword internal
print.BayesDiffIRTfit <- function(object){
  cat("Call:\n")
  print(object$call)
  summary(object)
}
