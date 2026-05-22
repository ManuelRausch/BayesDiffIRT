m <- 1
s <- 0.25

meanlog <- log(m)- 0.5* log(1 + s^2 / m^2)
meanlog
sdlog <- sqrt(log(1 + s^2/m^2))
sdlog

meanlog = 0
sdlog = 1

exp(meanlog + sdlog^2/2)
sqrt((exp(sdlog^2)-1)*exp(2*meanlog+sdlog^2))

gamma <- rlnorm(1000, 0, 0.5)
a <- rlnorm(1000, 0, 0.25)
summary(gamma/a)
