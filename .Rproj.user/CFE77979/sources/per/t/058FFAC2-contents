#' @title Sample from the Posterior Predictive diffusion item-response theory models

#' @importFrom rtdists rdiffusion


posteriorPredict <- function(object, ndraws){

}

# mit rtdists::rdiffusion simulieren


nSbj <- 200
nItems <- 10

ndt <- rlnorm(nSbj, -1.25, 0.3) # non-decision time
a <- rlnorm(nItems, 0, .25)      # item time pressure
theta <- rnorm(nSbj, 0, 2) # person drift variability
gamma <- rlnorm(nSbj, 0, 0.5) # person drift variability
nu <- rnorm(nItems, 0, 2) # item drifts

quantile(rlnorm(nSbj, 0, 0.5) / runif(nItems, 0.5, 2) )
quantile(rnorm(nSbj, 0, 2) - rnorm(nItems, 0, 2) )

SimData <- expand.grid(item = 1:nItems, sbj = 1:nSbj,
                       resp = NA, rt = NA)

library(RWiener)
for (i in 1:nrow(SimData)){
  print(paste0("alpha = ",
               round(gamma[SimData$sbj[i]] / a[SimData$item[i]],2),
               ", ndt = ",
               round(ndt[SimData$sbj[i]],2),
               ", delta = ",
               round(theta[SimData$sbj[i]]  - nu[SimData$item[i]],2)))

  try(sim <- rwiener(
    1,
    alpha = gamma[SimData$sbj[i]] / a[SimData$item[i]],
    tau = ndt[SimData$sbj[i]],
    beta = 0.5,
    delta = theta[SimData$sbj[i]] - nu[SimData$item[i]]))
  SimData$resp[i] <- ifelse(sim$resp == "upper", "yes"=1, "no"=0)
  SimData$rt[i] <- sim$q
}
SimData$resp
