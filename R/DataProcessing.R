make_stan_data <- function(rt, resp, sbj, item, model, priors){

  stan_data <- list(
    nObs = length(rt),
    nPerson = length(unique(sbj)),
    nItem = length(unique(item)),
    person = as.integer(sbj),
    item = as.integer(item),
    rt = as.numeric(rt),
    resp = as.integer(resp)
  )
  priors <- complete_priors(priors, model = model)

  c(stan_data,  priors_to_stan_data(priors))
}

priors_to_stan_data <- function(priors){
  out <- list()
  for (pr in priors) {
    parsed <- parse_dist(pr$dist)
    cls <- pr$class
    out[[paste0(cls, "_prior_family")]] <- prior_family_code(parsed$family)
    out[[paste0(cls, "_prior_par1")]] <- parsed$par1
    out[[paste0(cls, "_prior_par2")]] <- parsed$par2
  }
  out
}

parse_dist <- function(dist_call) {
  fam <- as.character(dist_call[[1]])
  args <- as.list(dist_call[-1])
  num <- function(x) {
    if (is.numeric(x)) return(as.numeric(x))
    eval(x, envir = baseenv())
  }
  if (fam == "normal") {
    if (length(args) != 2) {
      stop("normal() must have two arguments: mean and sd", call. = FALSE)
    }
    return(list(
      family = "normal",
      par1 = num(args[[1]]),
      par2 = num(args[[2]])
    ))
  }
  if (fam == "lognormal") {
    if (length(args) != 2) {
      stop("lognormal() must have two arguments: meanlog and sdlog", call. = FALSE)
    }
    return(list(
      family = "lognormal",
      par1 = num(args[[1]]),
      par2 = num(args[[2]])
    ))
  }
  if (fam == "gamma") {
    if (length(args) != 2) {
      stop("gamma() must have two arguments: shape and rate", call. = FALSE)
    }
    return(list(
      family = "gamma",
      par1 = num(args[[1]]),
      par2 = num(args[[2]])
    ))
  }
  stop("Unsupported prior family: ", fam, call. = FALSE)
}

prior_family_code <- function(family) {
  switch(family,
         normal = 1L,
         lognormal = 2L,
         gamma = 3L,
         stop("Unsupported prior family: ", family, call. = FALSE)
  )
}

