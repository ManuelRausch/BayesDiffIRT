make_stan_data <- function(data, rt, resp, sbj, item, priors, eps = 1e-3) {

  person_levels <- unique(data[[sbj]])
  item_levels   <- unique(data[[item]])

  person_id <- match(data[[sbj]], person_levels)
  item_id   <- match(data[[item]], item_levels)

  tau_upper <- vapply(seq_along(person_levels), function(p) {
    min(data[[rt]][person_id == p]) - eps
  }, numeric(1))

  if (any(tau_upper <= 0)) {
    stop("Some tau_upper values are non-positive.", call. = FALSE)
  }

  # Critical consistency check
  if (!all(tau_upper[person_id] < data[[rt]])) {
    bad <- which(!(tau_upper[person_id] < data[[rt]]))
    stop(
      "Subject indexing mismatch between 'person' and 'tau_upper'. ",
      "First failing row: ", bad[1],
      call. = FALSE
    )
  }

  c(
    list(
      nObs = as.integer(nrow(data)),
      nPerson = as.integer(length(person_levels)),
      nItem = as.integer(length(item_levels)),
      person = as.integer(person_id),
      item = as.integer(item_id),
      rt = as.numeric(data[[rt]]),
      resp = as.integer(data[[resp]]),
      tau_upper = as.numeric(tau_upper)
    ),
    priors_to_stan_data(priors)
  )
}


priors_to_stan_data <- function(priors){
  out <- list()

  for (pr in priors) {

    parsed <- parse_dist(pr$dist)
    cls <- pr$class
    validate_prior_family(cls, parsed$family)

    out[[paste0(cls, "_prior_family")]] <-
      prior_family_code(parsed$family, class = cls)
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
  if (fam == "uniform") {
    if (length(args) != 2) {
      stop("uniform() must have two arguments: min and max", call. = FALSE)
    }
    return(list(
      family = "uniform",
      par1 = num(args[[1]]),
      par2 = num(args[[2]])
    ))
  }

  stop("Unsupported prior family: ", fam, call. = FALSE)
}

prior_family_code <- function(family, class) {

  if (class %in% c("omega_theta", "omega_gamma", "tnd", "a")) {
    return(
      switch(family,
             lognormal = 1L,
             normal    = 2L,
             uniform   = 3L,
             stop(paste0("Unsupported prior family: ",
                  family, " for ", class), call. = FALSE)))
  }
  if (class %in% c("nu")) {
    return(
      switch(family,
             normal  = 1L,
             uniform = 2L,
             stop(paste0("Unsupported prior family: ",
                         family, " for ", class), call. = FALSE)))
    }
  stop("Unsupported prior class: ", class, call. = FALSE)
}

allowed_prior_families <- function(class) {
  switch(class,
         omega_theta = c("lognormal", "normal", "uniform"),
         omega_gamma = c("lognormal", "normal", "uniform"),
         nu          = c("normal", "uniform"),
         a           = c("lognormal", "normal", "uniform"),
         tnd         = c("lognormal", "normal", "uniform"),
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
  if (any(df[[rt]] <= 1e-3)) {
    stop("'", rt, "' must contain only positive values larger than 1e-3.", call. = FALSE)
  }

  if (any(df[[rt]] < 0.1)) {
    warning(
      "Some response times are below 0.1 seconds. ",
      "This may cause instability for person-specific non-decision times.",
      " Consider filtering fast response times. "
    )
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


