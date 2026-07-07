# This functions that help to create the package

# Use roxygen to write the Rd-Files for documentation.
devtools::document()

# Load complete package
devtools::load_all()

# check how the created documentation looks like
Rdpack::viewRd("man/fitBayesDiffIRT.Rd")

# Build and check a package
devtools::check()
devtools::check(manual = TRUE) # necessary for CRAN, tedious in the beginning,

# remove the old package
remove.packages("BayesDiffIRT")

# install the new package
devtools::install()



