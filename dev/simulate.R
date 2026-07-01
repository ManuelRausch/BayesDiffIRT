#' @title title
#' #' @param model `character`. Name of the diffusion item-response theory model to fit.
#' Currently implemented models:
#'  * "d" for the D-diffusion model (for survey items, default),
#'  * "dRV" for the D-diffusion model with random variability (for survey items),
#'  * "q" for the Q-diffusion model (for ability tests),
#'  * "qRV" for the Q-diffusion model with random variability (for ability tests).
#'
#' @export

simDiffIRT <- function(model = "d",
                       theta, gamma, nu, a, tnd, s_beta = 0, s_delta= 0){
  nItems <- length(theta)
  nSbjs <- length(gamma)

  out <- data.frame(rt = NA, resp=NA,sbj=1:nSbj, item=1:nItem)
  f <- switch(model, "d" =
                simDiffIRTd,
              "q"=simDiffIRTq,
              "dRV" = simDiffIRTdRV,
              "qRV" = simDiffIRTqRV)
  res <- maply(f, sbj=1:nSbj, item = 1:nItem)
}

simDiffIRTd <-
