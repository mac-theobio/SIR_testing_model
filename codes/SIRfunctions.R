# Functions for the SIR model
# By: A Gharouni

# install remotes package if necessary:
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

unpack <- McMasterPandemic::unpack

# Testing functions
sigma <- function(state,params){
  unpack(as.list(c(state,params)))
  W <- W_S*S_u+W_I*I_u+W_R*R_u
  testing_rate <- rho*N0/W
  return(testing_rate)
}

# SIR model
sir.model <- function(time,state,params){
  unpack(as.list(c(state,params)))
  ## Force of Infection
  Lambda <- beta * (I_u + eta_w*(I_n+I_p) + eta_c*I_t)/N0
  ## scaling the weights
  # W <- W_S*S_u+W_I*I_u+W_R*R_u
  # sigma <- rho*N0/W
  sigma <- sigma(state,params)

  # testing intensity
  F_S <- sigma*W_S
  F_I <- sigma*W_I
  F_R <- sigma*W_R
  # Equations
  dS_u.dt <- -Lambda*S_u - (1-p_S) * F_S * S_u + omega * S_n
  dS_n.dt <- -Lambda*S_n + (1-p_S) * F_S * S_u - omega * S_n
  dI_u.dt <-  Lambda*S_u + omega * I_n - F_I * I_u - gamma * I_u
  dI_n.dt <-  Lambda*S_n + (1-p_I) * F_I * I_u  - omega * I_n - gamma * I_n
  dI_p.dt <-  p_I * F_I * I_u - omega * I_p - gamma * I_p
  dI_t.dt <-  omega * I_p - gamma * I_t
  dR_u.dt <- gamma * I_u + omega * R_n - F_R * R_u
  dR_n.dt <- gamma * I_n + (1-p_R) * F_R * R_u - omega * R_n
  dR_p.dt <- gamma * I_p + p_R * F_R * R_u - omega * R_p
  dR_t.dt <- gamma * I_t + omega * R_p
  dN.dt <- omega * (S_n + I_n + R_n)
  dP.dt <- omega *(I_p + R_p)

  # return the rate of change
  dxdt <- c(dS_u.dt,dS_n.dt,dI_u.dt,dI_n.dt,dI_p.dt,dI_t.dt,dR_u.dt,dR_n.dt,dR_p.dt,dR_t.dt,dN.dt,dP.dt)
  ## }
  return(list(dxdt))
}

# TODO: fix the sigma the issue is when DFE(state,params), the multiroot() gives error.
DFE <- function(S,params){
  unpack(as.list(params))
  # S: Susceptibles
  S_u <- S[1]
  S_n <- S[2]
  I_u <-0
  R_u <-0
    #Weighted untested people
  W <- W_S*S_u+W_I*I_u+W_R*R_u
  sigma <- rho*N0/W
  F_S <- sigma*W_S

  return(c(eq1=S_u+S_n-N0,
           eq2=-F_S*S_u+omega*S_n))
}

# awaiting for negative tested Susceptible at DFE
Sn_dfe <-function(params){
  unpack(as.list(params))
  return((rho*(1-p_S)*N0)/omega)
}

# untested Susceptible at DFE
Su_dfe <-function(params){
  unpack(as.list(params))
  return(N0-Sn_dfe(params))
}

# Model simulation function
run.sir <- function(model, params,state,sim_time){
  # use update(params,beta=2)
  library(deSolve)
  unpack(as.list(c(state,params)))
  out <- as.data.frame(
    ode(
    func=model,
    y=state,
    times= sim_time,
    parms=params
  ))
  return(out)
}


F_I <- function(state,params){
  # Returns F_I at the specified state
  unpack(as.list(c(state,params)))
  sig <- sigma(state,params)
  return(sig*W_I)
}
## Fi at DFE
Fi_hat<-function(params){
  unpack(as.list(params))
  return(omega*rho/(omega-rho)*W_I/W_S)
}

# Basic Reproduction Number
R0<-function(params){
  unpack(as.list(params))
  Sn <- Sn_dfe(params)
  Su <- Su_dfe(params)
  Fi <- Fi_hat(params) #at DFE

  A<-gamma*(omega+gamma)^2+eta_w*gamma*(omega+gamma)*Fi+eta_c*omega*(omega+gamma)*p_I*Fi
  B<-gamma*omega*(omega+gamma)+eta_w*(gamma*(omega+gamma)*(Fi+gamma)+gamma*omega*p_I*Fi)+eta_c*omega^2*p_I*Fi
  C<-(gamma*(omega+gamma)+Fi*(gamma+omega*p_I))*(omega+gamma)
  return(beta/(N0*gamma*C)*(A*Su+B*Sn*eta_w))
}

## calculate the principle eigenvalue and the corresponding eigenvector.
eigvec_max<-function(params){
  unpack(as.list(params))
  Sn<-Sn_dfe(params)
  Su<-Su_dfe(params)
  Fi<-Fi_hat(params)
  colvec<-matrix(c(Su,Sn,0,0),4,1)
  rowvec<-matrix(c(1,eta_w,eta_w,eta_c),1,4)
  Fmat<-beta/N0*colvec%*%rowvec
  vmat<-matrix(c(Fi+gamma,-omega,0,0,
                 (p_I-1)*Fi,omega+gamma,0,0,
                 -p_I*Fi,0,omega+gamma,0,
                 0,0,-omega,gamma
  ),4,4,byrow = TRUE)
  G<- Fmat %*% solve(vmat)
  ev<-eigen(G)
  ev_max<- ev$vectors[,which(ev$values==max(ev$values))]
  if(all(ev_max<=0)){ev_max<--ev_max}
  return(ev_max)
}

# make grid and calc R0 for plotting
make_params_dat<-function(params,
                          rho_s,rho_e, ## range of parameters start to end
                          omega_s,omega_e,
                          eta_ws,eta_we,
                          eta_cs,eta_ce,
                          tol=1e-8,## tol: resolve the issue of very small numbers in plotting and when rho is close to omega
                          n_out=5,## n_out: facets (rows/cols)
                          n_in=41 ## n_in: grid N within facets
){
  unpack(as.list(params))
  eta_w <- if(eta_ws!=eta_we)seq(eta_ws,eta_we,length.out=n_out) else eta_ws
  eta_c <- if(eta_cs!=eta_ce)seq(eta_cs,eta_ce, length.out=n_out) else eta_cs
  rho <- if(rho_s!=rho_e)seq(rho_s,rho_e, length.out=n_in) else rho_s## note omega must be > rho
  omega <- if(omega_s!=omega_e)seq(omega_s,omega_e,length.out=n_in) else omega_s
  ## form the initial data frame with key parameters (model parameters)
  df1 <- expand.grid(N0=N0,beta=beta,gamma=gamma,
                     omega=omega,rho=rho,
                     W_S=W_S,W_I=W_I,W_R=W_R,
                     p_S=p_S,p_I=p_I,p_R=p_R,
                     eta_w=eta_w,eta_c=eta_c)
  ## principle eigenvector
  # eigvec <- t(apply(df1,1,function(params_in)eigvec_max(params=params_in)))
  dfout <-(df1 %>%
             dplyr::mutate(
               R0=apply(df1,1,function(params_in)R0(params=params_in)),
               R0_sub=ifelse(eta_w<eta_c, NA, R0),
               theta_w=1-eta_w,
               theta_c=1-eta_c,
               Delta=ifelse(eta_w<eta_c, NA, 1-(R0*gamma/beta) ),
               Delta=ifelse(abs(Delta)<tol,0,Delta)
               # I_u=eigvec[,1],
               # I_n=eigvec[,2],
               # I_p=eigvec[,3],
               # I_c=eigvec[,4]
             ))
  return(dfout)
}

catt <- function(...,file="modeldefs.tex") {
    cat(...,file=file,append=TRUE)
}

##' @param x numeric value
##' @param nm name of macro (will be \nm in LaTeX)
##' @param fmt number format (e.g. number of decimal places)
latexout <- function(x,nm,fmt="%1.1f", ...) {
    catt(sprintf("\\newcommand{\\%s}{%s}\n",nm,sprintf(fmt,x)), ...)
}

saveEnvironment()

