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
  c(stan_data,  priors_to_stan_data(priors))
}

priors_to_stan_data <- function(priors){
  out <- list()

  for (pr in priors) {

    parsed <- parse_dist(pr$dist)
    cls <- pr$class
    validate_prior_family(cls, parsed$family)

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

allowed_prior_families <- function(class) {
  switch(class,
         omega_theta = c("lognormal", "halfnormal", "uniform"),
         omega_gamma = c("lognormal", "halfnormal", "uniform"),
         nu          = c("normal", "uniform"),
         a           = c("lognormal", "halfnormal", "uniform"),
         tnd         = c("lognormal", "halfnormal", "uniform"),
         stop("Unknown prior class: ", class, call. = FALSE)
  )
}

validate_prior_family <- function(class, family) {
  allowed <- allowed_prior_families(class)
  if (!family %in% allowed) {
    stop(
      "Prior family '", family, "' is not supported for class '", class,
      "'. Allowed families are: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
}

validate_data <- function(data, rt, resp, sbj, item, na.rm = TRUE){
  # 1. Check if  columns exist
  vars <- c(rt, resp, sbj, item)

  missing_vars <- vars[!vars %in% names(data)]
  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing in 'data': ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  # 2. Remove useless columns
  df <- data[, vars, drop = FALSE]

  # 3. Check for NA values
  na_rows <- !stats::complete.cases(df)

  if (any(na_rows)) {
    if (na.rm) {
      df <- df[!na_rows, , drop = FALSE]
    } else {
      stop(
        "Missing values detected in the data. ",
        "Consider seting  na.rm = TRUE to remove incomplete cases.",
        call. = FALSE
      )
    }
  }

  # 4. Check data Type types

    # RT must be numeric and > 0
  if (!is.numeric(df[[rt]])) {
    stop("'", rt, "' must be numeric.", call. = FALSE)
  }
  if (any(df[[rt]] <= 0)) {
    stop("'", rt, "' must contain only positive values.", call. = FALSE)
  }

    # Response must be 0/1 or logical
  if (is.logical(df[[resp]])) {
    df[[resp]] <- as.integer(df[[resp]])
  }

  if (!is.numeric(df[[resp]])) {
    stop("'", resp, "' must be numeric (0/1) or logical.", call. = FALSE)
  }

  if (!all(df[[resp]] %in% c(0, 1))) {
    stop("'", resp, "' must contain only 0 and 1.", call. = FALSE)
  }

    # sbj and item: allow anything coercible to factor
  df[[sbj]]  <- as.factor(df[[sbj]])
  df[[item]] <- as.factor(df[[item]])

  # 5. Return cleaned data
  return(df)
}


