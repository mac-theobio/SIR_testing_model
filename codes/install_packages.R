while (!require(remotes)) {
  install.packages("remotes")
}
## install development version of bbmle:
if (!require("bbmle") || packageVersion("bbmle") < "1.0.23.5") {
  remotes::install_github("bbolker/bbmle")
}
## install the target package and all its dependencies:
while (!require(McMasterPandemic)) {
  remotes::install_github("bbolker/McMasterPandemic",
                          dependencies = TRUE,
                          build_vignettes = TRUE
  )
}
if (!require("shellpipes") ) {
  remotes::install_github("dushoff/shellpipes")
}

library(shellpipes)
rpcall("SIRfunctions.Rout SIRfunctions.R")

