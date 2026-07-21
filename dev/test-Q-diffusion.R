### Test the D-Diffusion model
# 0) Simulate some data
# 1) diffIRT fits for comparison
# 2) sample from the posterior using the simulated data set
# 3) Check diagnostics
# 4) Plot traces
# 5 Plot posteriors
# 6) Plot parameter recovery
# 7) Plot posterior predictive checks

rm(list=ls())

nSbj <- 200
nItems <- 10
load("dev/Tests-Q-Diffusion.RData")

ndt <- rlnorm(nSbj,-1.25, 0.3) # non-decision time
a <- rlnorm(nItems, 0, 0.25)      # item time pressure
theta <- rlnorm(nSbj, 0, 0.5) # person drift variability
gamma <- rlnorm(nSbj, 0, 0.5) # person drift variability
nu <- rlnorm(nItems, 0, 0.25) # item drifts

quantile(gamma/a)
quantile(theta / nu)


SimData <- expand.grid(item = 1:nItems, sbj = 1:nSbj,
                       resp = NA, rt = NA)

library(RWiener)
pb <- txtProgressBar(min = 0, max = nrow(SimData), style = 3)
for (i in 1:nrow(SimData)){
  setTxtProgressBar(pb, i)
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
    delta = theta[SimData$sbj[i]] / nu[SimData$item[i]]))
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
  diffIRT::diffIRT(rt[,-1], x[,-1], model="Q", se=F)
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

# 2) sample from the posterior using the simulated data set


library(BayesDiffIRT)
samples <-
  fitBayesDiffIRT(SimData,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "q",
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
  labs(x = "true parameter",
       y = "posterior mean",
       title = "non-decision time")


# 7) Plot posterior predictive checks

yrep <- posteriorPredict(samples, ndraws=20)

ppCheck(samples, type = "response", group = "item", yrep=yrep)
ppCheck(samples, type = "response", group = "person", yrep=yrep)
ppCheck(samples, type = "response", group = "none", yrep=yrep)

ppCheck(samples, type = "rtQuantile",
        minN=10, yrep=yrep)
ppCheck(samples, type = "rtQuantile", minN=10, yrep=yrep,
        group="item", index=1:5)

plotResponseSurface(samples, items = c(1:10), nrow=2)


save(samples, SimData, file = "dev/Tests-Q-Diffusion.RData")
