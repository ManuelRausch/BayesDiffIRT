# 0) load example data

data(extraversion, package="diffIRT")


# 1) For comparision: diffIRT package

x=extraversion[,1:10]
rt=extraversion[,11:20]

fit_diffIRT <- diffIRT::diffIRT(rt, x, model="D")




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


# 1) test thr data preprocessing functions

data <- BayesDiffIRT:::validate_data(
  data = Extra,rt = "rt", resp = "resp", sbj = "sbj", item = "item",
  na.rm = TRUE)
head(data, n = 25)

priors <- BayesDiffIRT:::complete_priors(
  priors = list(BayesDiffIRT:::prior(lognormal(0, 1), class = "tnd")),
  model= "d")

stan_data <- BayesDiffIRT:::make_stan_data(
  data, "rt", "resp", "sbj","item",priors = priors)


# 2) sample from the posterior

samples <-
  fitBayesDiffIRT(subset(Extra, sbj %in% 1:7),
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "d",
                  seed = 0, chains = 3,
                  parallel_chains = 3,
                  iter_warmup =  10^3,
                  iter_sampling = 10^3)


test <-  fitBayesDiffIRT(subset(Extra, sbj %in% 7),
                           rt = "rt", resp = "resp", sbj = "sbj",
                           item = "item", model = "d",
                           seed = 0, chains = 3,
                           parallel_chains = 3,
                           iter_warmup =  1000,
                           iter_sampling = 1000)

}

summary(samples)

save(file = "dev/Tests.RData")
