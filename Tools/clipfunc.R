# Data Sequences Simulation Function
ClipSimulation <- function(par         = NULL,
                           K           = NULL,
                           FullDesignX = NULL,
                           Ts          = NULL,
                           seed        = NULL)
{
  ci <- par[1:(K-1)]
  len_par <- length(par)
  phi <- par[len_par]
  if(!is.null(seed)){set.seed(seed)}
  
  Z <- rep(0, Ts)
  et <- arima.sim(n=Ts, list(ar=phi), sd=sqrt(1-phi^2)) # sd argument is for WN
  if(length(par)==K){
    mst <- rep(0, Ts)
  }else{
    mst <- FullDesignX%*%par[K:(len_par-1)]
  }
  Z <- et + mst
  
  clip <- function(X=NULL){return(length(which(X>ci))+1)}
  X_hour <- sapply(X=as.vector(Z), FUN = clip)
  
  X_hour_wide <- matrix(0, nrow=Ts, ncol=K)
  for (h in 1:Ts) {X_hour_wide[h,X_hour[h]] <- 1}
  
  res <- list(Z, X_hour, X_hour_wide)
  names(res) <- c("Z", "X_hour", "X_hour_wide")
  return(res)
}

UniCpDetRcpp <- function(res_sim = NULL, 
                         DesignXEst = NULL, 
                         stepsize=NULL, 
                         BoundAdjust=NULL)
{
  
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
  
  myBd <- c(-0.99, 0.99)
  if(phi_initial < BoundAdjust){
    myBd <- c(-0.99, phi_initial)
  }
  # if(phi_initial < BoundAdjust){
  #   if(phi_initial < 0){
  #     myBd <- c(-0.99, phi_initial)
  #   }else{
  #     myBd <- c(phi_initial, 0.99)
  #   }
  # }else{
  #   myBd <- c(-0.99, 0.99)
  # }
  OptimSecond2 <- optim(par=phi_initial, fn=diffMOM, method = "Brent",
                        lower=myBd[1], upper=myBd[2], control=list(reltol=1e-05),
                        Marg_est=Marg_est, K=K, DesignX=DesignXEst,
                        MySampleCov=MySampleCov)
  phi_est <- OptimSecond2$par
  par_est <- c(Marg_est, phi_est)
  
  #---------------- Innovation ----------------#
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
  
  #--------------- One Step Ahead Prediction -----------#
  MyPred <- ClipPred(EXL=EXL, X_hour=X_hour, mylag=mylag, HQ=HQ)
  It    <- MyPred$Innovation/sqrt(V)
  
  #--------------- CUSUM Test Statistics ---------------#
  MyCp  <- MyCusum(ErrErr=MyPred$Innovation, V=V)
  CpCsm <- MyCp$Csm
  CpLoc <- MyCp$Location
  CpVal <- MyCp$Maximum
  
  #--------------- Result Return ---------------#
  res <- list(MyPred$Innovation, It, MyPred$Prediction, 
              CpCsm, CpLoc, CpVal, par_est, res_NW$hs)
  names(res) <- c("Innovation", "It", "Predictions", 
                  "CpCsm", "CpLoc", "CpVal", "ParEst", "Hessian")
  return(res)
}

#----------------------------- Server Use -----------------------------#
# One job for each core
one.job <- function(job=NULL, nTasks=NULL, nCore=NULL, iii=NULL,
                    Ts=NULL,
                    K=NULL,
                    parT=NULL,
                    DesignXT=NULL,
                    DesignXEst=NULL, 
                    stepsize=NULL, 
                    BoundAdjust=NULL)
{
  
  nSubtasks <- round(nTasks/nCore)
  RES <- data.frame(job=job, task=1:nSubtasks, tim=rep(NA, nSubtasks),
                    seed=rep(NA, nSubtasks),
                    loc=rep(NA, nSubtasks), val=rep(NA, nSubtasks))
  ParamMatrix <- matrix(0, nrow=nSubtasks, ncol=K+dim(DesignXEst)[2])
  RES <- cbind(RES, ParamMatrix)
  
  for(subid in 1:nSubtasks){
    
    tim.start <- Sys.time()
    #################################################################
    #----------------------- Calculation ---------------------------#
    #################################################################
    myseed <- 20191111 + job*nSubtasks + subid + iii*100000
    res_sim <- ClipSimulation(par=parT, K=K, FullDesignX=DesignXT, Ts=Ts,
                              seed = myseed)
    UniCpRes <- UniCpDetRcpp(res_sim=res_sim, DesignXEst=DesignXEst,
                             stepsize=stepsize, BoundAdjust=BoundAdjust)
    tim.end <- Sys.time()
    timeused <- difftime(tim.end, tim.start, units="secs")
    cat("\n No.job", job, "task", subid,
        "completed within", timeused, "at seed", myseed, "!")
    #################################################################
    #---------------------- Result Return --------------------------#
    #################################################################
    RES$seed[subid]  <- myseed
    RES$tim[subid]   <- timeused
    RES$loc[subid]   <- UniCpRes$CpLoc
    RES$val[subid]   <- UniCpRes$CpVal
    RES[subid,7:(6+length(UniCpRes$ParEst))] <- UniCpRes$ParEst
  }
  return(RES)
}


# Paraell computing implementation
useMultiCore <- function(nTasks=NULL, nCore=NULL, iii=NULL,
                         Ts=NULL,
                         K=NULL,
                         parT=NULL,
                         DesignXT=NULL,
                         DesignXEst=NULL,
                         stepsize=NULL, 
                         BoundAdjust=NULL)
{
  
  cat("Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  tim.start = Sys.time()
  FinalRes <- foreach(i=1:nCore, .combine = "rbind") %dopar%
    one.job(job=i, nTasks=nTasks, nCore=nCore, iii=iii,
            Ts=Ts, K=K, parT=parT, DesignXT=DesignXT, DesignXEst=DesignXEst,
            stepsize=stepsize, BoundAdjust=BoundAdjust)
  tim.end = Sys.time()
  cat("Done.\n")
  cat("\n\n nTasks =", nTasks, "\t nCore =", nCore,
      "\t Aveg. time =", mean(FinalRes$tim),
      "\t Total Time =", difftime(tim.end, tim.start, units="hours"),
      "hours \n\n")
  return(FinalRes)
}


# Frequency analysis
myfreq <- function(x)
{
  
  require(astsa)
  
  K <- dim(x)[2]
  Var = var(x) # var-cov matrix
  xspec = mvspec(x, spans=c(7,7), plot=FALSE) 
  fxxr = Re(xspec$fxx) # fxxr is real(fxx)
  # compute Q = Var^-1/2
  ev = eigen(Var)
  Q = ev$vectors%*%diag(1/sqrt(ev$values))%*%t(ev$vectors)
  # compute spec envelope and scale vectors
  num = xspec$n.used # sample size used for FFT
  nfreq = length(xspec$freq) # number of freqs used
  specenv = matrix(0,nfreq,1) # initialize the spec envelope 
  beta = matrix(0,nfreq,K) # initialize the scale vectors 
  for (k in 1:nfreq){
    ev = eigen(2*Q%*%fxxr[,,k]%*%Q/num, symmetric=TRUE)
    specenv[k] = ev$values[1] # spec env at freq k/n is max evalue 
    b = Q%*%ev$vectors[,1] # beta at freq k/n
    beta[k,] = b/sqrt(sum(b^2)) 
  } # helps to normalize beta
  
  # output and graphics
  frequency = xspec$freq
  plot(frequency, 100*specenv, ylim=c(0,1), 
       type="l", ylab="Spectral Envelope (%)") # add significance threshold to plot
  m = xspec$kernel$m
  etainv = sqrt(sum(xspec$kernel[-m:m]^2)) 
  thresh=100*(2/num)*exp(qnorm(.9999)*etainv)
  abline(h=thresh, lty=6, col=4)
}