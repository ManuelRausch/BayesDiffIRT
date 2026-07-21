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
  cmdstanVersion <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(e) NULL
  )

  if (is.null(cmdstanVersion)) {
    stop(
      "CmdStan could not be found. Install it using ",
      "cmdstanr::install_cmdstan().",
      call. = FALSE
    )
  }

  tryCatch(
    cmdstanr::check_cmdstan_toolchain(quiet = TRUE),
    error = function(e) {
      stop(
        "A working C++ toolchain is required to compile the Stan models.\n",
        conditionMessage(e),
        "\nSee https://mc-stan.org/docs/cmdstan-guide/installation.html#cpp-toolchain.",
        call. = FALSE
      )
    }
  )

  invisible(TRUE)

}
