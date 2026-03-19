library(diffIRT)

data(extraversion)
x=extraversion[,1:10]
rt=extraversion[,11:20]

# fit an unconstrained D-diffusion model
 # item-specific a, item-specific v, item-specific Ter
 # pupulation soecific parameters: omega(gamma), omega 7heta
  # Why is there an item-specific Ter??? Should Ter not be a person-parameter?



res1=diffIRT(rt,x,model="D")
fs <- factest(res1, se=TRUE)
fs
Extra <- data.frame(resp = x, rt = rt)
