# All files with the purpose to validate some stuff

checkCmdstan <- function() {

  # Check cmdstanr package
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "The 'cmdstanr' package is not installed.\n",
      "Install it with install.packages('cmdstanr').",
      call. = FALSE
    )
  }

  # Check CmdStan installation
  cs <- try(cmdstanr::cmdstan_version(), silent = TRUE)

  if (inherits(cs, "try-error") || is.null(cs)) {
    stop(
      "CmdStan is not installed or not configured.\n",
      "Run cmdstanr::install_cmdstan() to install it.\n",
      "See https://mc-stan.org/cmdstanr/ for details.",
      call. = FALSE
    )
  }

  invisible(TRUE)

}
