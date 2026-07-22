rm(list=ls())
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


priors <- list(prior(beta(1, 1), class = "s_beta"))

start.time <- Sys.time()
fitQRV <-
  fitBayesDiffIRT(Data,
                  rt = "rt", resp = "resp", sbj = "sbj",
                  item = "item", model = "qRV",
                  nChains = 10,  nCores = 10,
                  nWarmup =  1000, nSamples = 1000,
                  refresh = 100)

FittingTimeQRV <- Sys.time() - start.time
cat("Runtime QRV:", round(as.numeric(FittingTimeQRV,units = "hours"), 2), "hours\n")

gg1 <- ppCheck(fitQRV, type = "rtQuantile")

start.time <- Sys.time()
fitQfixed <- fitBayesDiffIRT(Data,
                rt = "rt", resp = "resp", sbj = "sbj",
                item = "item", model = "q",
                nChains = 10,  nCores = 10,
                nWarmup =  1000, nSamples = 1000,
                refresh = 100)
FittingTimeQ <- Sys.time() - start.time

gg2 <- ppCheck(fitQfixed, type = "rtQuantile")



save(fitQRV, FittingTimeQRV, gg1,
     fitQfixed , FittingTimeQ, gg2,

     file = "dev/Tests-Q-Diffusion-RV.RData")
