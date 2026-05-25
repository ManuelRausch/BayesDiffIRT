# 0) Simulate some data

rm(list=ls())
load("dev/Tests.RData")

nSbj <- 200
nItems <- 10

ndt <- runif(nSbj, .1, .5) # non-decision time
a <- runif(nItems, 0.5, 2)      # item time pressure
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
  data = SimData,rt = "rt", resp = "resp", sbj = "sbj", item = "item",
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
                  n.chains = 3,  n.cores = 3,
                  n.warmup =  10^3,
                  n.samples = 10^3)
summary(samples)
checkDiagnostics(samples)
plot(samples, parameter = "theta", type = "trace")
plot(samples, parameter = "gamma", type = "trace")
plot(samples, parameter = "nu", type = "trace")
plot(samples, parameter = "a", type = "trace")
plot(samples, parameter = "tnd", type = "trace")
plot(samples, parameter = "omega_gamma",
     type = "trace")
plot(samples, parameter = "omega_theta",
     type = "trace")


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



# Plot parameter recovery

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
    data.frame(true = gamma,
               estimated =
                 apply(samples$fit$draws("ndt"), 3, mean)),
    aes(x = true, y = estimated)) +
  geom_smooth(method="lm", se=FALSE) +
  geom_point() +
  ggpmisc::stat_correlation(
    ggpmisc::use_label("R", "R.CI"),
    small.r=TRUE) +
  theme_classic() +
  labs(x = "true parameter", y = "posterior mean",
       title = "person ability")

save.image(file = "dev/Tests.RData")

# fit the data set from the diffIRT package

data(extraversion, package = "diffIRT")
x <- extraversion[,1:10]; rt <- extraversion[,11:20]
x <- cbind(1:nrow(x), x)
colnames(x) <- c("sbj", paste0("Item", 1:10))
x <- tidyr::pivot_longer(as.data.frame(x), cols=Item1:Item10,
                         names_to="item",
                         values_to = "resp")


rt <-  cbind(1:nrow(rt), rt)
colnames(rt) <- c("sbj", paste0("Item", 1:10))
rt <- tidyr::pivot_longer(as.data.frame(rt), cols=Item1:Item10,
                          names_to="item",
                          values_to = "rt")
Extra <- merge(x, rt)
Extra$item <- factor(Extra$item)

samples2 <-
  fitBayesDiffIRT(Extra,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "d",
                  n.chains = 3,  n.cores = 3,
                  n.warmup =  10^3,
                  n.samples = 10^3)
summary(samples2)
checkDiagnostics(samples2)

class(samples2)
plot(samples2, parameter = "gamma", type = "interval",
           index=1:10)
plot(samples2, parameter = "nu", type = "interval")
plot(samples2, parameter = "a", type = "interval")
plot(samples2, parameter = "tnd", type = "interval")

plot(samples2, parameter = "omega_gamma", type = "distribution")
plot(samples2, parameter = "omega_gamma", type = "distribution")
plot(samples2, parameter = "omega_theta", type = "distribution")

