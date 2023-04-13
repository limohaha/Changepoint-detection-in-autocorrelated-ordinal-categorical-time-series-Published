library(mvtnorm)
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)
options("scipen"=10) # show all digits

# Self-defined function
sourceCpp("Tools/myfuncClipPaper.cpp")

Varin2006Sim <- function(par=NULL, Ts=NULL, K=NULL, myseed=NULL){
  Ts <- Ts + 100 # first 100 observations will be dropped
  len_par <- length(par)
  
  ci   <- par[1:(K-1)]
  beta <- par[K:(len_par-1)]  
  phi  <- par[len_par]
  
  if(!is.null(myseed)){set.seed(myseed)}
  
  B0 <- rep(1, Ts)
  B1 <- rnorm(Ts, mean=-1,    sd=1)
  B2 <- rnorm(Ts, mean=-0.25, sd=sqrt(0.0324))
  DesignXT <- cbind(B0, B1, B2)
  mu <- DesignXT%*%beta
  
  yt <- rep(0, Ts)
  yt[1] <- 0
  eps <- rnorm(Ts, mean=0, sd=1)
  for(t in 2:(length(yt)-1)){
    yt[t] <- mu[t] + phi*yt[t-1] + eps[t]
  }
  
  clip <- function(X=NULL){return(length(which(X>ci))+1)}
  
  # vector form
  X_hour <- sapply(X=as.vector(yt), FUN = clip)
  # matrix form
  X_hour_wide <- matrix(0, nrow=Ts, ncol=K)
  for (h in 1:Ts) {X_hour_wide[h,X_hour[h]] <- 1}
  
  # delete frist 100 observations and list outputs
  res_sim <- list(DesignXT=DesignXT[101:Ts,], yt=yt[101:Ts], 
                  X_hour=X_hour[101:Ts], X_hour_wide=X_hour_wide[101:Ts,])
  return(res_sim)
}


nSim <- 500
myseed <- 12345678

Ts <- 2000
K  <- 7
ciT    <- c(0, 1.2, 2.2, 3.1, 4.1, 5.3)
beta0T <-  2.9
beta1T <- -0.6
beta2T <-  9
phiT   <- -0.6
parT   <- c(ciT, beta0T, beta1T, beta2T, phiT)


MyRes <- matrix(0, nrow=nSim, ncol=length(parT)-1)
Comptime <- rep(0, nSim)

# par_est <- c(0.9208970,  1.6861138,  2.3216814,  3.0813160,  
#              4.0812630,  2.2443289, -0.4690252,  7.0049321,  0.4945051)

for(iii in 1:nSim){
  
  tim1 <- Sys.time()
  res_Varin2006Sim <- Varin2006Sim(par=parT, Ts=Ts, K=K, myseed=myseed+iii)
  
  res_sim <- list(Z           = res_Varin2006Sim$yt, 
                  X_hour      = res_Varin2006Sim$X_hour,
                  X_hour_wide = res_Varin2006Sim$X_hour_wide)
  DesignXEst <- res_Varin2006Sim$DesignXT
  
  X_hour <- res_sim$X_hour
  X_hour_wide <- res_sim$X_hour_wide
  Ts <- dim(X_hour_wide)[1]
  K <- dim(X_hour_wide)[2]
  
  cout <- summary(factor(X_hour))
  ci_initial <- as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
  ci_initial <- ci_initial - ci_initial[1]
  ci_initial <- ci_initial[-1]
  par_initial <- c(ci_initial,
                   rep(0, dim(DesignXEst)[2]),
                   acf(X_hour, plot=F)$acf[2,1,1])
  
  constrLSE <- matrix(0, nrow=K-2, ncol=length(par_initial))
  constrLSE[1,1] <- 1 
  for(ii in 2:(K-2)){constrLSE[ii,c(ii-1, ii)] <- c(-1,1)} # ci
  constrLSE_ci <- rep(0, dim(constrLSE)[1])
  
  OptimFirst <- constrOptim(par_initial, f=LSE_cpp1,
                            method = "Nelder-Mead", ui=constrLSE,
                            ci=constrLSE_ci, hessian=F,
                            control= list(reltol=1e-07),
                            X=X_hour_wide, DesignX=DesignXEst)
  par_est <- OptimFirst$par
  MyRes[iii,] <- par_est
  
  tim2 <- Sys.time()
  Comptime[iii] <- round(difftime(tim2, tim1, units = "secs"))
  cat("\n No.", iii, "simulation completed within", round(difftime(tim2, tim1, units = "secs")), "seconds!")
}


# MyRes <- MyRes[1:100,]

EstTable <- matrix(0, nrow=4, ncol=length(parT)-1)
colnames(EstTable) <- c(paste("C", 2:(K-1), sep=""), "beta0", "beta1", "beta2", "phi")
row.names(EstTable) <- c("True", "Mean", "Rel.bias", "sd")

EstTable[1,] <- parT[2:length(parT)]
EstTable[2,] <- colMeans(MyRes)
EstTable[3,] <- (EstTable[1,] - EstTable[2,])/EstTable[1,]
EstTable[4,] <- apply(MyRes, 2, sd)
EstTable

#------------------------------------- phi = 0.5  -------------------------------------#

#                 C2           C3           C4           C5           C6        beta0        beta1        beta2         phi
# True      1.20000000  2.200000000  3.100000000  4.100000000  5.300000000  2.900000000 -0.600000000  9.000000000  0.50000000
# Mean      1.21412752  2.207036917  3.109790280  4.116039188  5.331469948  2.915144751 -0.602267224  9.046650771  0.50057981
# Rel.bias -0.01177293 -0.003198598 -0.003158155 -0.003911997 -0.005937726 -0.005222328 -0.003778707 -0.005183419 -0.00115961
# sd        0.08154747  0.102406468  0.123189216  0.143270393  0.182315668  0.116993949  0.040701021  0.332587729  0.01461114

#                    C2           C3           C4           C5           C6        beta0        beta1        beta2          phi
# True      1.200000000  2.200000000  3.100000000  4.100000000  5.300000000  2.900000000 -0.600000000  9.000000000 0.5000000000
# Mean      1.205433670  2.207609252  3.113534177  4.118919794  5.324876285  2.918925328 -0.601951470  9.054381089 0.4998463376
# Rel.bias -0.004528059 -0.003458751 -0.004365864 -0.004614584 -0.004693639 -0.006525975 -0.003252449 -0.006042343 0.0003073248
# sd        0.074961302  0.096359703  0.114553469  0.136552358  0.172418543  0.110989487  0.038559761  0.336120806 0.0146419638


#------------------------------------- phi = -0.6  -------------------------------------#
#                    C2           C3           C4          C5          C6        beta0         beta1        beta2           phi
# True      1.200000000  2.200000000  3.100000000  4.10000000  5.30000000  2.900000000 -0.6000000000  9.000000000 -0.6000000000
# Mean      1.204534415  2.207157653  3.116591899  4.12321633  5.33761177  2.924491554 -0.6003886096  9.072304526 -0.5997616936
# Rel.bias -0.003778679 -0.003253479 -0.005352225 -0.00566252 -0.00709656 -0.008445364 -0.0006476827 -0.008033836  0.0003971773
# sd        0.062986730  0.094615435  0.122122188  0.15955544  0.21516254  0.136006488  0.0427878876  0.364974678  0.0143986391


#############################################################
# Second check #
#############################################################

#------------------------------------- phi = -0.6  -------------------------------------#
# C2           C3          C4           C5           C6        beta0        beta1        beta2          phi
# True      1.200000000  2.200000000  3.10000000  4.100000000  5.300000000  2.900000000 -0.600000000  9.000000000 -0.600000000
# Mean      1.205880170  2.209080583  3.11381813  4.115005779  5.330502923  2.918573579 -0.598733186  9.039775100 -0.598938447
# Rel.bias -0.004900141 -0.004127538 -0.00445746 -0.003659946 -0.005755269 -0.006404682  0.002111357 -0.004419456  0.001769255
# sd        0.067890204  0.096789596  0.12577157  0.159554587  0.212105168  0.140719276  0.045315477  0.380297940  0.015188333
