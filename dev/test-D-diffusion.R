### Test the D-Diffusion model
# 0) Simulate some data
# 1) diffIRT fits for comparison
# 2) sample from the posterior using the simulated data set
# 3) Check diagnostics
# 4) Plot traces
# 5 Plot posteriors
# 6) Plot parameter recovery
# 7) fit the data set from the diffIRT package


# 0) Simulate some data

rm(list=ls())

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

# 1) For comparision: diffIRT package

rt <- tidyr::pivot_wider(
  SimData, id_cols = sbj,
  values_from = rt,
  names_from = item)

x <- tidyr::pivot_wider(
  SimData, id_cols = sbj,
  values_from = resp,
  names_from = item)


fit_diffIRT <-
  diffIRT::diffIRT(rt[,-1], x[,-1], model="D", se=F)
summary(fit_diffIRT )

library(ggplot2)
library(patchwork)
ggplot(
  data.frame(true = a, estimated = coef(fit_diffIRT)$item[,1]),
  aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  ggplot(
    data.frame(true = nu, estimated = coef(fit_diffIRT)$item[,2]),
    aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic()


# 1) test thr data preprocessing functions

data <- BayesDiffIRT:::validateData(
  data = SimData,rt = "rt", resp = "resp",
  sbj = "sbj", item = "item",
  na.rm = TRUE)
head(data, n = 25)

priors <- BayesDiffIRT:::completePriors(
  priors = list(BayesDiffIRT:::prior(lognormal(0, 1), class = "tnd")),
  model= "d")

stanData <- BayesDiffIRT:::makeStanData(
  data, "rt", "resp", "sbj","item",
  priors = priors)

minRT <- SimData |>
  dplyr::group_by(sbj) |>
  dplyr::summarise(minRt = min(rt),.groups = "drop")
stanData$tauUpper - minRT$minRt
summary(stanData$tauUpper)

# 2) sample from the posterior

library(BayesDiffIRT)
samples <-
  fitBayesDiffIRT(SimData,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "d",
                  nChains = 3,  nCores = 3,
                  nWarmup =  10^3)

summary(samples)

# 3) Check diagnostics
checkDiagnostics(samples)

# 4) Plot traces

plot(samples, parameter = "theta", type = "trace")
plot(samples, parameter = "gamma", type = "trace")
plot(samples, parameter = "nu", type = "trace")
plot(samples, parameter = "a", type = "trace")
plot(samples, parameter = "tnd", type = "trace")

plot(samples, parameter = "omega_gamma",
     type = "trace")
plot(samples, parameter = "omega_theta",
     type = "trace")

# 5 Plot posteriors

plot(samples, parameter = "theta", type = "interval")
plot(samples, parameter = "gamma", type = "interval")
plot(samples, parameter = "nu", type = "interval")
plot(samples, parameter = "a", type = "interval")
plot(samples, parameter = "tnd", type = "interval")

plot(samples, parameter = "omega_gamma",
     type = "density")
plot(samples, parameter = "omega_theta",
     type = "density")
plot(samples, parameter = "a",
     type = "density")


# 6) Plot parameter recovery

library(ggplot2)
library(patchwork)
ggplot(
  data.frame(true = a,
             estimated =
               apply(samples$fit$draws("a"), 3, mean)),
  aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  labs(x = "true parameter", y = "posterior mean",
       title = "item time pressure") +
  ggplot(
    data.frame(true = nu,
               estimated = apply(samples$fit$draws("nu"), 3, mean)),
    aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  ggtitle("item drift")

# person parameters
ggplot(
  data.frame(true = theta,
             estimated =
               apply(samples$fit$draws("theta"), 3,
                     mean)),
  aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  labs(x = "true parameter", y = "posterior mean",
       title = "person response caution") +
  ggplot(
    data.frame(true = gamma,
               estimated =
                 apply(samples$fit$draws("gamma"), 3, mean)),
    aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  labs(x = "true parameter", y = "posterior mean",
       title = "person ability") +
  ggplot(
    data.frame(true = ndt,
               estimated =
                 apply(samples$fit$draws("tnd"), 3, mean)),
    aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  labs(x = "true parameter", y = "posterior mean",
       title = "non-decision time")


# 7) fit the data set from the diffIRT package

library(tidyverse)
data(extraversion, package = "diffIRT")
Extra <- as.data.frame(extraversion)
names(Extra)[1:10]  <- paste0("Item", 1:10, "_resp")
names(Extra)[11:20] <- paste0("Item", 1:10, "_rt")
Extra$sbj <- 1:nrow(Extra)

Extra <- tidyr::pivot_longer(
  Extra,
  cols = tidyselect::matches("^Item\\d+_(resp|rt)$"),
  names_to = c("item", ".value"),
  names_pattern = "^(Item\\d+)_(resp|rt)$"
)

Extra$item <- factor(Extra$item)
Extra$item <- factor(Extra$item)

samples2 <-
  fitBayesDiffIRT(Extra,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "d",
                  nChains = 4,  nCores = 4,
                  nWarmup =  10^3,
                  nSamples = 10^3)
summary(samples2)
checkDiagnostics(samples2)

class(samples2)
plot(samples2, parameter = "gamma", type = "interval",
           index=1:10)
plot(samples2, parameter = "nu", type = "interval")
plot(samples2, parameter = "a", type = "interval")
plot(samples2, parameter = "tnd", type = "interval")

plot(samples2, parameter = "omega_gamma", type = "density")
plot(samples2, parameter = "omega_theta", type = "density")

ppCheck(samples2, type = "response", group = "item")
ppCheck(samples2, type = "response", group = "person")
ppCheck(samples2, type = "response", group = "none")

ppCheck(samples2, type = "rtQuantile",
        minN=10, ndraws=10)
ppCheck(samples2, type = "rtQuantile", minN=10, ndraws=10,
        group="item", index=1:5)

gg <- plotResponseSurface(samples2, item = 1,
                    theta.range = NULL,
                    gamma.range = NULL,
                    grid.size = 50,
                    contours = c(.1, .3, .5, .7, .9),
                    ndraws = 200)
gg

save(samples, samples2, file = "dev/Tests-D-Diffusion.RData")
load("dev/Tests-D-Diffusion.RData")


# Rainers function #

coef <-
  function(obj) {
    stopifnot(inherits(obj,"summary.BayesDiffIRTfit"))
    n = obj$dataInfo$nPerson
    k = obj$dataInfo$nItem
    splits =   cumsum(c(n,k,1,n,k,1,n,n,n,n))
    result = list()
    result$z_theta      = data.frame(obj$variables[           1 :splits[ 1],])
    result$nu           = data.frame(obj$variables[(splits[1]+1):splits[ 2],])
    result$omega_theta  = data.frame(obj$variables[(splits[2]+1):splits[ 3],])
    result$z_gamma      = data.frame(obj$variables[(splits[3]+1):splits[ 4],])
    result$a            = data.frame(obj$variables[(splits[4]+1):splits[ 5],])
    result$omega_gamma  = data.frame(obj$variables[(splits[5]+1):splits[ 6],])
    result$tnd          = data.frame(obj$variables[(splits[6]+1):splits[ 7],])
    result$theta        = data.frame(obj$variables[(splits[7]+1):splits[ 8],])
    result$log_gamma    = data.frame(obj$variables[(splits[8]+1):splits[ 9],])
    result$gamma        = data.frame(obj$variables[(splits[9]+1):splits[10],])

    return(result)
  } # end fun coef

coef(summary(samples))

coef.BayesDiffIRTfit <- function(
    object,
    estimate = "mean",
    parameters = NULL,
    ...) {

  if (!inherits(object, "BayesDiffIRTfit")) {
    stop("`object` must inherit from class \"BayesDiffIRTfit\".",
         call. = FALSE)
  }

  modelSummary <- summary(object, ...)
  variables <- modelSummary$variables

  requiredColumns <- c("variable", estimate)

  # Remove the index from names such as theta[1].
  baseName <- sub("\\[.*$", "", variables$variable)

  if (is.null(parameters)) {
    parameters <-
      c("theta","gamma","nu","a","tnd","omega_theta","omega_gamma")
  }

  keep <- baseName %in% parameters

  result <- variables[[estimate]][keep]
  names(result) <- variables$variable[keep]

  result
}

x<- coef.BayesDiffIRTfit(samples, estimate="median")
