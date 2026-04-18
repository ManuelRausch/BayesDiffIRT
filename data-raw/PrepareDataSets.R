
load("data-raw/RDK2.RData")


RDK <- data.frame(sbj = RDK2$participant,
                  resp = RDK2$correct,
                  item = RDK2$diffCond,
                  rt = round(RDK2$rt, 3) )

usethis::use_data(RDK, overwrite = TRUE, compress = "xz")

