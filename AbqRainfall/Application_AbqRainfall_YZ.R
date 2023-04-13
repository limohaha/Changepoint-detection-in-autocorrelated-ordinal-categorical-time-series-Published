library(lubridate)
library(TSA)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(readr)
library(ggplot2)
library(xtable)
library(car)
library(ggpubr)
library(qqplotr)
library(latex2exp)
library(astsa)

options("scipen"=10) # show all digits

# Self-defined function
sourceCpp("Tools/myfuncClipPaper.cpp")


UniCpDetect <- function(res_sim = NULL, 
                        DesignXEst = NULL,
                        stepsize=NULL, 
                        BoundAdjust=NULL){
  
  #--------------- Input objects ---------------#
  X_hour <- res_sim$X_hour
  X_hour_wide <- res_sim$X_hour_wide
  Ts <- dim(X_hour_wide)[1]
  K <- dim(X_hour_wide)[2]
  
  #----------------- Estimation ----------------#
  cout <- summary(factor(X_hour))
  ci_initial <- as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
  p <- which(!is.finite(ci_initial))
  if(length(p)!=0){
    if(p==1){
      ci_initial[p] <- -6
    }else{
      ci_initial[p] <-  6
    }
  }
  if(sum(DesignXEst)==0){
    # stationary
    par_initial <- ci_initial
  }else{
    # non-stationary
    fri_initial <- rep(0, dim(DesignXEst)[2])
    par_initial <- c(ci_initial, fri_initial)
  }
  
  #------- Marginal Parameter Estimation
  res_NW <- NW_cpp(par=par_initial, X=X_hour_wide, DesignX=DesignXEst,
                   stepsize=stepsize, conv=1e-05)
  Marg_est <- res_NW$par
  MyHS <- res_NW$hs
  
  #------- Dependence Paramater Estimation
  phi_initial <- acf(X_hour, plot=F)$acf[2,1,1]
  MySampleCov <- sum(diag(t(X_hour_wide[1:(Ts-1),])%*%X_hour_wide[2:Ts,]))/Ts
  # Bound Adjustment
  myBd <- c(-0.99, 0.99)
  if(phi_initial < BoundAdjust){
    myBd <- c(-0.99, phi_initial)
  }
  OptimSecond2 <- optim(par=acf(X_hour, plot=F)$acf[2,1,1], fn=diffMOM,
                        method = "Brent", lower=-0.99, upper=0.99,
                        control=list(reltol=1e-05), Marg_est=Marg_est,
                        K=K, DesignX=DesignXEst, MySampleCov=MySampleCov)
  phi_est <- OptimSecond2$par
  par_est <- c(Marg_est, phi_est)
  
  #---------------- CUSUM on Y ----------------#
  seg_est <- c(-Inf, par_est[1:(K-1)], Inf)
  # Univariate Expectation
  if(length(par_est)==K){
    mst <- rep(0, Ts)
  }else{
    mst <- DesignXEst%*%par_est[K:(length(par_est)-1)]
  }
  EXL <- rep(0, Ts)
  for(t in 1:Ts){EXL[t] <- K - sum(pnorm(seg_est-mst[t])[2:K])}
  InnovRes <- UniInnovRcpp(EXL=EXL, mst=mst, Marg_est = Marg_est,
                           phi_est = phi_est, K=K, numCoef=1500)
  HQ <- InnovRes$HQ
  V <- InnovRes$V
  mylag <- InnovRes$lag
  rm(InnovRes)
  # One Step Ahead Prediction
  MyPred <- ClipPred(EXL=EXL, X_hour=X_hour, mylag=mylag, HQ=HQ)
  ItY    <- MyPred$Innovation/sqrt(V)
  # CUSUM Test Statistics
  MyCp   <- MyCusum(ErrErr=MyPred$Innovation, V=V)
  CpCsmY <- MyCp$Csm
  CpLocY <- MyCp$Location
  CpValY <- MyCp$Maximum
  
  #---------------- CUSUM on Z ----------------#
  # Estimation situation with tau1hat modified
  # Latent variable reconstruction
  ResLatent <- Re_latent(par = par_est, DesignX = DesignXEst, X = X_hour_wide)
  Z_expct <- ResLatent$Z_expct
  mst <- ResLatent$mst
  epsilonHat1 <- c(0,(Z_expct[2:Ts] - mst[2:Ts]) - phi_est*(Z_expct[1:(Ts-1)] - mst[1:(Ts-1)]))
  ResCusum1 <- rep(0, Ts)
  for(tau in 1:Ts){
    ResCusum1[tau] <- (sum(epsilonHat1[1:tau])-tau/Ts*sum(epsilonHat1[1:Ts]))/sqrt(Ts)
  }
  # tau1 variance approximation
  qn <- floor((Ts)^(1/3))
  epsilonbar <- mean(epsilonHat1)
  tmpsum <- 0
  for(s in 1:qn){
    tmpsum1 <- 0
    for(t in 1:(Ts-s)){
      tmpsum1 <- tmpsum1 + (epsilonHat1[t]-epsilonbar)*(epsilonHat1[t+s]-epsilonbar)
    }
    # tmpsum <- tmpsum + (1-s/(qn+1))*tmpsum1/(Ts-s)
    tmpsum <- tmpsum + (1-s/(qn+1))*tmpsum1/Ts
  }
  tau1hat <- sum((epsilonHat1-epsilonbar)^2)/Ts + 2*tmpsum
  sigmahat1 <- sqrt(tau1hat)
  ItZ <- abs(ResCusum1/sigmahat1)
  CpLocZ <- which.max(ItZ)
  CpValZ <- max(ItZ)
  
  res <- list(CpLocY, CpValY, ItY, CpCsmY, CpLocZ, CpValZ, ItZ, par_est, res_NW$hs)
  names(res) <- c("CpLocY", "CpValY", "ItY", "CpCsmY", "CpLocZ", "CpValZ", "ItZ", "ParEst", "Hessian")
  return(res)
}

# Expectation of Truncated Normal
expct_TN <- function(mu=NULL, sigma=NULL, cl=NULL, cu=NULL){
  return(mu+sqrt(sigma)*(dnorm(x=cl,mean=mu,sd=sqrt(sigma))-dnorm(x=cu,mean=mu,sd=sqrt(sigma)))/
           (pnorm(q=cu,mean=mu,sd=sqrt(sigma))-pnorm(q=cl,mean=mu,sd=sqrt(sigma))))
}

# Latent Variable Reconstruction
Re_latent <- function(par=NULL, DesignX=NULL, X=NULL){
  
  Ts <- dim(X)[1]
  K <- dim(X)[2]
  
  # Parameters Extraction
  ci <-  par[1:(K-1)]
  seg <- c(-Inf, ci, Inf)
  len_par <- length(par)
  phi <- par[len_par]
  
  if(length(par)==K){
    mst <- rep(0, Ts)
  }else{
    mst <- DesignX%*%par[K:(len_par-1)]
  }
  
  sig2 <- 1 # latent variable variance fix to be 1
  # allocate
  Z_expct <- rep(0, Ts)
  for(h in 1:Ts){
    pp <- which(X[h,]==1)
    Z_expct[h] <- expct_TN(mu=mst[h], sigma=sig2, cl=seg[pp], cu=seg[pp+1])
  }
  
  res <- list(mst, Z_expct)
  names(res) <- c("mst", "Z_expct")
  return(res)
}

##################################################
#------------ Real Data Application -------------# Discretized Rainfall
##################################################

# 0. Read and subset data
Abq <- read_csv("Abq19310301_20201211.csv", col_types = cols(DATE = col_date(format = "%Y-%m-%d")))
Abq <- as.data.frame(Abq[,c("DATE", "PRCP")])
# select data range by changing starting and ending values
pstart <- which(as.character(Abq$DATE)== "1990-01-01")
pend   <- which(as.character(Abq$DATE)== "2000-12-31")
mydat <- Abq[pstart:pend,]

# 1. Missing values Detection and Imputation

# 1.1 Date missing
FullSeq <- data.frame(DateFull=seq.Date(from = min(mydat$DATE), to = max(mydat$DATE), by = 1))
Missing <- FullSeq$DateFull[!FullSeq$DateFull %in% mydat$DATE]
Missing

# Imputation (& delete leap days)
if(length(Missing)>0){
  # Leave the missing values slots
  mydat <- merge(mydat, FullSeq, by.x= "DATE", by.y="DateFull", all.y = TRUE)
  # To obtain an exact period of 365 days, 
  # delete leap days observations (February 29).
  p2 <- which(format(mydat$DATE, format = "%m%d") == "0229")
  mydat <- mydat[-p2,]
  # calculation 365 day's average among 11 years
  DaysMean <- data.frame(date=format(seq.Date(from = as.Date("2019-01-01", "%Y-%m-%d"), 
                                              to   = as.Date("2019-12-31", "%Y-%m-%d"), 
                                              by = 1),
                                     format="%m%d"), 
                         mean=colMeans(t(matrix(mydat$PRCP, nrow=365)), na.rm = T))
  # Impute day's average for missing date
  p3 <- which(is.na(mydat$PRCP))
  # by days mean
  for(i in 1:length(p2)){
    pp <- which(DaysMean$date == format(mydat$DATE[p2[i]], format="%m%d"))
    mydat$PRCP[p2[i]] <- DaysMean$mean[pp]
  }
  # # by Overall mean
  # p4 <- which(is.na(mydat$PRCP))
  # mydat$PRCP[p3] <- mean(mydat$PRCP, na.rm = T)
}else{
  # To obtain an exact period of 365 days, 
  # delete leap days observations (February 29).
  p2 <- which(format(mydat$DATE, format = "%m%d") == "0229")
  mydat <- mydat[-p2,]
}

# 1.2 Value Missing
p4 <- which(is.na(mydat$PRCP))
p4

# 1.3 Discretized to 3 categories
mydat$Catg <- 0
mydat$Catg[which(mydat$PRCP > 0 & mydat$PRCP < 0.2)] <- 1
mydat$Catg[which(mydat$PRCP >=0.2)]                  <- 2

# # 1.3 Discretized to 5 categories
# mydat$Catg <- 0
# mydat$Catg[which(mydat$PRCP >  0   & mydat$PRCP < 0.03)] <- 1
# mydat$Catg[which(mydat$PRCP >= 0.03 & mydat$PRCP < 0.1)] <- 2
# mydat$Catg[which(mydat$PRCP >= 0.1 & mydat$PRCP < 0.2)] <- 3
# mydat$Catg[which(mydat$PRCP >= 0.2)]                    <- 4

# count of each category
summary(factor(mydat$Catg))
#    0    1    2 
# 3312  520  183 


X_hour <- mydat$Catg + 1
Ts <- length(X_hour)
K <- max(X_hour)
X_hour_wide <- matrix(0, nrow=Ts, ncol=K)
for (t in 1:Ts) {X_hour_wide[t,X_hour[t]] <- 1}
res_sim <- list(X_hour=X_hour, X_hour_wide=X_hour_wide)


##################################################
#------------- Frequency Analysis ---------------#
##################################################

# 1. Discretized Rainfall by periodogram
haha <- periodogram(y=mydat$Catg, plot = F)
p <- order(haha$spec, decreasing = T)
# related frequency
1/haha$freq[p[1:10]]
# [1] 368.181818 184.090909 119.117647  75.000000 115.714286 810.000000  21.891892
# [8]  17.608696   5.752841   3.986220

# 2. Continuous Rainfall by periodogram
haha <- periodogram(y=mydat$PRCP, plot = F)
p <- order(haha$spec, decreasing = T)
# related frequency
1/haha$freq[p[1:10]]
# [1] 368.181818   4.258675  24.695122  75.000000  94.186047  14.516129   3.986220
# [8]  21.891892  17.608696 184.090909


##################################################
#---------------- clip model fit ----------------#
##################################################
ss1 <- 365
ss2 <- 183
TrendValue <- 1:Ts/Ts
BValue <- cos(2*pi*(1:Ts)/ss1)
DValue <- sin(2*pi*(1:Ts)/ss1)
EValue <- cos(2*pi*(1:Ts)/ss2)
FValue <- sin(2*pi*(1:Ts)/ss2)

#------------------- First Fit
DesignXEst11 <- cbind(TrendValue, BValue, DValue, EValue, FValue)
FitRes11 <- UniCpDetect(res_sim = res_sim, DesignXEst = DesignXEst11, 
                        stepsize = 1, BoundAdjust = -0.5)
par_est1 <- FitRes11$ParEst

######### Parameter SE
mySE1 <- sqrt(diag(solve(FitRes11$Hessian)))
parMatrix1 <- cbind(par_est1[1:(length(par_est1)-1)],
                    par_est1[1:(length(par_est1)-1)] - 2*mySE1,
                    par_est1[1:(length(par_est1)-1)] + 2*mySE1)
parMatrix1 <- rbind(parMatrix1, c(par_est1[length(par_est1)], NA, NA))
colnames(parMatrix1) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix1) <- c(paste("C", 1:(K-1), sep=""), 
                          "Beta", "B", "D", "E", "F", "Phi")
round(parMatrix1, digits=4)
#        Param  LowerB  UpperB
# C1    0.8355  0.7311  0.9398
# C2    1.6298  1.4607  1.7989
# Beta -0.2559 -0.4684 -0.0434
# B    -0.1776 -0.2804 -0.0749
# D    -0.1729 -0.2821 -0.0637
# E     0.0835 -0.0204  0.1875
# F     0.1752  0.0720  0.2783
# Phi   0.4145      NA      NA

# xtable(t(parMatrix1), digits = 4)
# xtable(matrix(mySE1, nrow=1), digits = 4)


######### Detection Result on CUSUM_Y
round(FitRes11$CpValY, digits = 3)
# [1] 1.045
CpLocY <- FitRes11$CpLocY
# [1] 2366
CUSUM_Y <- FitRes11$CpValY
# [1] 1.044532
mydat[FitRes11$CpLocY, ]
# 1996-06-25

TestStatistic <- 1.045


######### Detection Result on CUSUM_Z
round(FitRes11$CpValZ, digits = 3)
# [1] 0.91
CpLocZ <- FitRes11$CpLocZ
# [1] 2421
CUSUM_Z <- FitRes11$CpValZ
# [1] 0.9102201
mydat[FitRes11$CpLocZ, ]
# 1996-08-19





#------------------- Model Re-fit based on CUSUM_Y
tauestY <- CpLocY
CpValueY <- c(rep(0, tauestY-1), rep(1, Ts-tauestY+1))
DesignXEst12Y <- cbind(DesignXEst11, CpValueY)
FitRes12Y <- UniCpDetect(res_sim = res_sim, DesignXEst = DesignXEst12Y, 
                         stepsize = 1, BoundAdjust = -0.5)
par_est2Y <- FitRes12Y$ParEst
mySE2Y <- sqrt(diag(solve(FitRes12Y$Hessian)))
round(mySE2Y, digits=4)

parMatrix2Y <- cbind(par_est2Y[1:(length(par_est2Y)-1)],
                     par_est2Y[1:(length(par_est2Y)-1)] - 2*mySE2Y,
                     par_est2Y[1:(length(par_est2Y)-1)] + 2*mySE2Y)
parMatrix2Y <- rbind(parMatrix2Y, c(par_est2Y[length(par_est2Y)], NA, NA))
colnames(parMatrix2Y) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix2Y) <- c(paste("C", 1:(K-1), sep=""), 
                           "Beta", "B", "D", "E", "F", "Delta", 
                           "Phi")
round(parMatrix2Y, digits=4)
#         Param  LowerB  UpperB
# C1     0.7448  0.6345  0.8551
# C2     1.5361  1.3689  1.7033
# Beta  -0.6928 -1.0858 -0.2998
# B     -0.1799 -0.2827 -0.0770
# D     -0.1660 -0.2748 -0.0573
# E      0.0795 -0.0241  0.1831
# F      0.1669  0.0637  0.2701
# Delta  0.3097  0.0463  0.5730
# Phi    0.4159      NA      NA

seg_estY <- c(-Inf, par_est2Y[1:(K-1)], Inf)
# Univariate Expectation
if(length(par_est2Y)==K){
  mstY <- rep(0, Ts)
}else{
  mstY <- DesignXEst12Y%*%par_est2Y[K:(length(par_est2Y)-1)]
}
EXLHaY <- rep(0, Ts)
for(t in 1:Ts){EXLHaY[t] <- K - sum(pnorm(seg_estY-mstY[t])[2:K])}


#------------------- Model Re-fit based on CUSUM_Z
tauestZ <- CpLocZ
CpValueZ <- c(rep(0, tauestZ-1), rep(1, Ts-tauestZ+1))
DesignXEst12Z <- cbind(DesignXEst11, CpValueZ)
FitRes12Z <- UniCpDetect(res_sim = res_sim, DesignXEst = DesignXEst12Z, 
                          stepsize = 1, BoundAdjust = -0.5)
par_est2Z <- FitRes12Z$ParEst
mySE2Z <- sqrt(diag(solve(FitRes12Z$Hessian)))
round(mySE2Z, digits=4)

parMatrix2Z <- cbind(par_est2Z[1:(length(par_est2Z)-1)],
                     par_est2Z[1:(length(par_est2Z)-1)] - 2*mySE2Z,
                     par_est2Z[1:(length(par_est2Z)-1)] + 2*mySE2Z)
parMatrix2Z <- rbind(parMatrix2Z, c(par_est2Z[length(par_est2Z)], NA, NA))
colnames(parMatrix2Z) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix2Z) <- c(paste("C", 1:(K-1), sep=""), 
                           "Beta", "B", "D", "E", "F", "Delta", 
                           "Phi")
round(parMatrix2Z, digits=4)
#         Param  LowerB  UpperB
# C1     0.7356  0.6248  0.8463
# C2     1.5311  1.3632  1.6991
# Beta  -0.7387 -1.1335 -0.3438
# B     -0.1888 -0.2924 -0.0853
# D     -0.1694 -0.2787 -0.0601
# E      0.0876 -0.0165  0.1917
# F      0.1727  0.0695  0.2760
# Delta  0.3487  0.0822  0.6152
# Phi    0.4114      NA      NA


seg_estZ <- c(-Inf, par_est2Z[1:(K-1)], Inf)
# Univariate Expectation
if(length(par_est2Z)==K){
  mstZ <- rep(0, Ts)
}else{
  mstZ <- DesignXEst12Z%*%par_est2Z[K:(length(par_est2Z)-1)]
}
EXLHaZ <- rep(0, Ts)
for(t in 1:Ts){EXLHaZ[t] <- K - sum(pnorm(seg_estZ-mstZ[t])[2:K])}






################################################
## Count summary before and after changepoint ##
################################################

#--------- tau_Y count ---------#
summary(factor(mydat$Catg[1:(tauestY-1)]))
#    0    1    2 
# 1948  323   94 
summary(factor(mydat$Catg[tauestY:Ts]))
#    0    1    2 
# 1364  197   89 
#--------- tau_Z count ---------#
summary(factor(mydat$Catg[1:(tauestZ-1)]))
#    0    1    2
# 1990  330  100
summary(factor(mydat$Catg[tauestZ:Ts]))
#    0    1    2 
# 1322  190   83 

#--------- tau_Y proportion ---------#
round(summary(factor(mydat$Catg[1:(tauestY-1)])) / sum(summary(factor(mydat$Catg[1:(tauestY-1)]))), digits = 3)
#     0     1     2 
# 0.824 0.137 0.040 
round(summary(factor(mydat$Catg[tauestY:Ts])) / sum(summary(factor(mydat$Catg[tauestY:Ts]))), digits = 3)
#     0     1     2 
# 0.827 0.119 0.054 
#--------- tau_Z proportion ---------#
round(summary(factor(mydat$Catg[1:(tauestZ-1)])) / sum(summary(factor(mydat$Catg[1:(tauestZ-1)]))), digits = 3)
#     0     1     2 
# 0.822 0.136 0.041 
round(summary(factor(mydat$Catg[tauestZ:Ts])) / sum(summary(factor(mydat$Catg[tauestZ:Ts]))), digits = 3)
#     0     1     2 
# 0.829 0.119 0.052 











##################################################
#------------- Introduction Figure 1  -----------# Continuous and Discretized Time Series
##################################################

# setEPS()
# postscript("Figure1.eps", width = 14, height = 10)
# setEPS()
# postscript("PPTFigure1.eps", width = 14, height = 10)
par(cex.lab=1.2, cex.axis=1.2, mfrow=c(2,1))
tmp <- c(1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000)
plot(x=1:Ts, y=mydat$PRCP, type="l", col="darkgrey",
     xaxt="n", xlim=c(-100, Ts+100), xaxs="i",
     xlab="Year", ylab="Rainfall in inches")
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(as.numeric(format(mydat$DATE, "%Y")) == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)

plot(x=1:Ts, y=mydat$Catg+1, type="l", cex.lab=1.2,
     xlab="Year", ylab="Rainfall categories",
     xaxt="n", yaxt="n",
     xlim=c(-100, Ts+100), ylim=c(1,3), xaxs="i")
axis(1, at = px, labels = tmp, las = 1, cex.axis=1)
axis(2, at = c(1,2,3), cex.axis=1)
# dev.off()



##################################################
#-------------- Application Figure 1  -----------# Time Series with Expected value as mean structure
##################################################
# setEPS()
# postscript("ApplyFigure1.eps", width = 14, height = 7)
par(cex.lab=1.2, cex.axis=1.2, mfrow=c(1,1))
plot(x=1:Ts, y=mydat$Catg+1, type="l", col="darkgrey",
     xlab="Year", ylab="Categorized Rainfall",
     xaxt="n", yaxt="n", 
     xlim=c(-100, Ts+100), ylim=c(1,3.2), xaxs="i")
tmp <- c(1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000)
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(as.numeric(format(mydat$DATE, "%Y")) == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)
axis(2, at = c(1,2,3))
lines(1:(tauestY-1), 
      EXLHaY[1:(tauestY-1)], 
      type="l", lty="solid", col="red", lwd=3)
lines(tauestY:Ts,
      EXLHaY[tauestY:Ts], 
      type="l", lty="dashed", col="blue", lwd=3)
legend("top", horiz = T, legend=c("Before Changepoint", "After Changepoint"),
       col=c("red", "blue"), lty=c("solid", "dashed"), lwd =c(3,3), cex=1, bty="n")
# dev.off()




####################################
#----- Application Figure 2  ------# ACF of Prediction Residuals for First fit
####################################
# setEPS()
# postscript("ApplyFigure2.eps", width = 14, height = 7)
par(cex.lab=1.2, cex.axis=1.2, mfrow=c(1,1))
plot(acf(FitRes11$ItY, plot=F, lag.max = 60), ylim=c(-0.04, 0.04), main="")
# dev.off()

####################################
#----- Application Figure 3  ------# CUSUM Test Statistics Plot on Y
####################################
# setEPS()
# postscript("ApplyFigure3.eps", width = 14, height = 7)
par(cex.lab=1, cex.axis=1.2, mfrow=c(1,1))
plot(x=1:Ts, y=FitRes11$CpCsm, type="l", col="black",
     xlab="Year", ylab=expression('CUSUM'['I'](tau)),
     xaxt="n",
     xlim=c(-100, Ts+100), xaxs="i")
abline(h=0.925, lty="dashed", col="black")
text(x=350, y=0.98, labels="95th Percentile")
tmp <- c(1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000)
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(as.numeric(format(mydat$DATE, "%Y")) == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)
# dev.off()


####################################
#----- Application Figure 4  ------# CUSUM Test Statistics Plot on Z
####################################
# setEPS()
# postscript("ApplyFigure4.eps", width = 14, height = 7)
par(cex.lab=1, cex.axis=1.2, mfrow=c(1,1))
plot(x=1:Ts, y=FitRes11$ItZ, type="l", col="black",
     xlab="Day", 
     # ylab=expression('CUSUM'['e'](tau)),
     ylab="CUSUM Statistics",
     xaxt="n",
     xlim=c(-30, Ts+30), ylim=c(0, 1), xaxs="i")
abline(h=0.910, lty="dashed", col="black")
text(x=350, y=0.96, labels="95th Percentile")
tmp <- c(1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000)
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(as.numeric(format(mydat$DATE, "%Y")) == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)
# dev.off()


