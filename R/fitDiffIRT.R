#' @title Fit Bayesian diffusion item-response theory models

#' @description  fitDiffIRT` function fits
#' @param data  a `data.frame` with different row indicating respones to different items
#' @param rt
#' @param response
#' @param id
#' @param model
#' @return An object of class `BayesDiffIRTfit`. Available methods include print, summary,
#' #' @author
#' Manuel Rausch, \email{manuel.rausch@aau.at}
#' @details
#' @references
#' Alexandrowicz, R. W. (2020). The diffusion model visualizer: An interactive tool to understand the diffusion model parameters. Psychological Research, 84(4), 1157-1165. https://doi.org/10.1007/s00426-018-1112-6
#' Molenaar, D., Tuerlinckx, F., & Maas, H. L. J. V. D. (2015). Fitting Diffusion Item Response Theory Models for Responses and Response Times Using the R Package diffIRT. Journal of Statistical Software, 66(4). https://doi.org/10.18637/jss.v066.i04
#' Van Der Maas, H. L. J., Molenaar, D., Maris, G., Kievit, R. A., & Borsboom, D. (2011). Cognitive psychology meets psychometric theory: On the relation between process models for decision making and latent variable models for individual differences. Psychological Review, 118(2), 339-356. https://doi.org/10.1037/a0022749

#' @export
fitDiffIRT <- function(data, rt = "rt", resp = "resp", sbj = "sbj",
                       item = "item", model = "Q", priors = NULL,
                       seed = 42, chains = 3, parallel_chains = 1,
                       iter_warmup =  1000,
                       iter_sampling = 2000
                       ){
  call <- match.call()

  # 1) Preprocessing. Prios should be included int othe data to be passed to STAN
  priors <- complete_priors(priors, model)
  stan_data <- make_stan_data(data, rt, resp, sbj, item, priors)

  # 2) Compile + sample (cmdstanr)
  stanfile = switch(model, "q" = "qdiffusion.stan",
                    "p" = "pdiffusion.stan",
                    stop("Error. Unknown model. Please select one out of the followng: c('d', 'q')"))
  mod <- cmdstan_model(stanfile)
  fit <- mod$sample(data = stan_data, seed = seed, chains = chains,
                    parallel_chains = parallel_chains,
                    iter_warmup = iter_warmup,
                    iter_sampling = iter_sampling)

  # 3) Postprocess
  # diagnostics <- collect_diagnostics(fit)   # don't know if model diagnostics can are created automatically.

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
#' @examples

