data(rotation, package = "diffIRT")
x <- rotation[,1:10]; rt <- rotation[,11:20]
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
Data <- merge(x, rt)
Data$item <- factor(Data$item)

samples <-
  fitBayesDiffIRT(Data,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "qRV",
                  prior = prior(beta(1,1), class = "s_beta"),
                  n.chains = 3,  n.cores = 3,
                  n.warmup =  10^3,
                  n.samples = 10^3)
summary(samples)
checkDiagnostics(samples)

class(samples)
plot(samples, parameter = "gamma", type = "interval",
     index=1:10)
plot(samples, parameter = "nu", type = "interval")
plot(samples, parameter = "a", type = "interval")
plot(samples, parameter = "tnd", type = "interval")

plot(samples, parameter = "omega_gamma", type = "density")
plot(samples, parameter = "omega_theta", type = "density")
plot(samples, parameter = "s_delta", type = "density")
plot(samples, parameter = "s_beta", type = "density")

save(samples, samples, file = "dev/Tests-D-Diffusion-RV.RData")
