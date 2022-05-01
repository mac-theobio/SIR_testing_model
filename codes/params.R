library(shellpipes)
rpcall("params.Rout params.R")
# SIRfunctions.rda
loadEnvironments()

params <- c(
N0=1000000, #total population size, 
beta= 0.5, #transmission rate 0.339
gamma=1/6, #recovery rate, was 1/3
omega=0.25, #test returning rate 
rho=0.01, #testing rate percapita (also used 0.8, 1/3)
W_S=1, W_I=1, W_R=1, #Testing relative weights
p_S=0, p_I=1, p_R=0.5, #test specificity, i.e., prob of waiting for being tested positive
eta_w=0.02, eta_c=0.01 #isolation parameter (eta=0 is the perfect isolation)
## s=2 # isolation ratio parameter, i.e., s=eta_w/eta_c
            )
            
class(params) <- "params_pansim"

saveVars(params)

