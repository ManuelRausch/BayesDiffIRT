
check_cmdstan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "The 'cmdstanr' package is not installed.\n",
      "Install it with install.packages('cmdstanr').",
      call. = FALSE
    )
  }

  cs <- cmdstanr::cmdstan_version(quiet = TRUE)
  if (!isTRUE(cs$ok)) {
    stop(
      "CmdStan is not installed or not configured.\n",
      "Run cmdstanr::install_cmdstan().",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
