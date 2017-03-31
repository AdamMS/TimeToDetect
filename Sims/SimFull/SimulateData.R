library(plyr)

load("../../OVEN/parmests.Rdata") # These are parameter estimates from most models (some still running)
load("../../OVEN/Data9.Rdata")    # To get covariate values from actual values and proper data structure
load("../SimZero/Simpars.Rdata")  # I'm just taking the intercpts (plural) and shape parms for peaked distributions
  # from the SimZero values.

set.seed(NULL)

######### Keeping for historical purposes.  Once run, it needn't be run again.
##### Simulation parameters for ALL models
##### We are using the 3-interval estimates (for time-saving and detection rate reasons)
ParValTabs <- N <- vector('list',14)
ModOrder <- c("f.exp","f.expmix",rep(c("f.gamma","f.gammamix","f.logN","f.logNMix","f.weibull","f.weibullmix"),each=2))
for(j in 1:length(ParValTabs)) {
  ParValTabs[[j]] <- parmests[parmests$modeloutput==ModOrder[j],
                             match(c("X50.","parameter","modeloutput"), names(parmests))]
}

##### Simulate Abundance
Xa <- dat$Xa           # Abundance covariates
Xd <- dat$Xd           # Detection covariates
ia <- dat$ia           # Index of abund RndEff levels
id <- dat$id           # Index of detection RndEff levels... ObsID
partable <- data.frame(ttd=c(rep("Expo",2),rep("Gamma",4),rep("LogN",4),rep("Weibull",4)),
                       mix=c(F,T,F,F,T,T,F,F,T,T,F,F,T,T),
                       peak=c(F,F,rep(c(F,T),6)))
dat9 <- dat3 <- uncounted <- vector('list',nrow(partable))

# Assign 'Peak' parameters:
# These are plug-in values for Intcpt_a, Intcpt_d, and Shape parameters
# I am just borrowing the appropriate values from the SimZero simulations, so...
# ...goals: p_global = 0.8, mode = 5, n_obs = 967, gamma=0.65
sub_list <- (1:14)[partable$peak]   # index of partable rows that are peak models
sub_params <- c("intcpt_a", "intcpt_d", "alpha", "sigma_det", "shape_k")
# This code inserts the above intercept and shape values into the parameter tables
for (j in sub_list) {
  subparmsmod <- sub_params[!is.na(ParValTabs[[j]]$X50.[match(sub_params, ParValTabs[[j]]$parameter)])]
  ParValTabs[[j]]$X50.[match(subparmsmod, ParValTabs[[j]]$parameter)] <- 
    as.numeric(Simpars[j, match(subparmsmod, names(Simpars))])
  ParValTabs[[j]]$X50.[ParValTabs[[j]]$parameter == "gamma"] <- 0.65
}

# Adjustments to shrink random effects and slightly boost some intcpt_d
for(j in c(4,6,8,10,12,14)){
  ParValTabs[[j]]$X50.[ParValTabs[[j]]$parameter=="sigma_d"] <- 
    ParValTabs[[j]]$X50.[ParValTabs[[j]]$parameter=="sigma_d"]/4
}
for(j in c(6,8,10,12,14)){
  ParValTabs[[j]]$X50.[ParValTabs[[j]]$parameter=="intcpt_d"] <- 
    ParValTabs[[j]]$X50.[ParValTabs[[j]]$parameter=="intcpt_d"]*1.1
}

##### HERE'S WHERE THE LOOP STARTS
# for generating multiple (or one) replicate(s) of simulations

for(Rep in 1:1){
  
  ##### Simulate random effect values
  for(j in 1:length(ParValTabs)){
    temp <- ParValTabs[[j]] # 'Cuz 'temp' is shorter to type
    ra <- c(rnorm(dat$n_ras[1], 0, temp$X50.[temp$parameter=="sigma_a[1]"]), 
            rnorm(dat$n_ras[2], 0, temp$X50.[temp$parameter=="sigma_a[2]"]))
    rd <- rnorm(dat$n_rds, 0, temp$X50.[temp$parameter=="sigma_d"])
    
    ##### Start Generating Data
    # Abundance regression
    l.lambda <- temp$X50.[temp$parameter=="intcpt_a"] + as.matrix(Xa) %*% temp$X50.[grepl("^ba",temp$parameter)]  # log(lambda)
    for(i in 1:2) l.lambda <- l.lambda + ra[ia[,i]]            # Incorporate random effects
    
    
    # Used to be a big mistake here... I SIMULATED DATA WITH Xa and ^ba instead of Xd and ^bd
    # Detection regression
    xbeta.plus <- temp$X50.[grepl("intcpt_d",temp$parameter)] + as.matrix(Xd) %*% temp$X50.[grepl("^bd",temp$parameter)] + rd[id] 
    # Simulated site abundances
    N[[j]] <- rpois(dat$n_sites, exp(l.lambda))
    
    
    
    
    tau <- c(0, dat$tau)
    # Options
    TTD  <- partable$ttd[j]
    Mix  <- partable$mix[j]
    
    # Temporary storage object for site(row) x interval(col) detection probability
    p.vector <- data.frame(matrix(NA, nrow=dat$n_sites, ncol=9))
    
    # Calcualte interval-specific detection probabilities
    if(TTD=="Expo"){
      # intcpt_d = log(beta) from desired exponential
      for(j2 in 1:9) p.vector[,j2] <- pexp(tau[j2+1],rate=exp(xbeta.plus)) - 
        pexp(tau[j2],rate=exp(xbeta.plus))
    } 
    
    if(TTD=="Gamma"){
      alpha <- temp$X50.[temp$parameter=="alpha"]
      # intcpt_d = log(beta) from desired gamma
      # Therefore, log(beta) = Xb + RndEff
      for(j2 in 1:9) p.vector[,j2] <- pgamma(tau[j2+1],alpha,rate=exp(xbeta.plus)) - 
        pgamma(tau[j2],alpha,rate=exp(xbeta.plus))
    } 
    
    if(TTD=="LogN"){
      sigmadet <- temp$X50.[temp$parameter=="sigma_det"]
      # intcpt_d = -mu from desired logNormal
      for(j2 in 1:9) p.vector[,j2] <- plnorm(tau[j2+1],-xbeta.plus,sigmadet) - 
        plnorm(tau[j2],-xbeta.plus,sigmadet)
    } 
    
    if(TTD=="Weibull"){
      shape_k <- temp$X50.[temp$parameter=="shape_k"]
      # intcpt_d = log(beta) from desired Weibull
      for(j2 in 1:9) p.vector[,j2] <- pweibull(tau[j2+1],shape_k,scale=exp(-xbeta.plus)) - 
        pweibull(tau[j2],shape_k,scale=exp(-xbeta.plus))
    } 
    
    if(Mix) for(i in 1:(dat$n_sites)) p.vector[i,] <- 
      temp$X50.[grepl("gamma",temp$parameter)]*p.vector[i,] + 
      c(1-temp$X50.[grepl("gamma",temp$parameter)],rep(0,8))
    
    Ival9 <- data.frame(matrix(NA,nrow=length(N[[j]]),ncol=10))  # Includes uncounted in Column 10
    for(ab in 1:length(N[[j]])) Ival9[ab,] <- t(rmultinom(1, N[[j]][ab], c(as.numeric(p.vector[ab,]),1-sum(p.vector[ab,]))))
    unc <- Ival9[,10]
    Ival9 <- Ival9[,-10]
    # Ival3 <- data.frame(X1=rowSums(Ival9[,1:2]), X2=rowSums(Ival9[,3:4]), X3=rowSums(Ival9[,5:9]))
    
    datpkg <- function(n_ints, dat=dat){
      newdat <- dat
      if(n_ints==3) {
        newdat$y <- Ival3
        newdat$tau <- c(3,5,10)
        newdat$n_ints <- 3
      }
      if(n_ints==9) newdat$y <- Ival9
      return(newdat)
    }
    # dat3[[j]] <- datpkg(3, dat)
    dat9[[j]] <- datpkg(9, dat)
    uncounted[[j]] <- unc
  }

  par(mfrow=c(3,5))
  for(j in 1:14) 
    plot(c(0:9+0.5),c(rep(sum(dat9[[j]]$y[,1])/2,2),colSums(dat9[[j]]$y[,2:9])), pch=16, main=j)
  
  # dir.create(paste0("Sim_",Rep))
  # save(dat3, dat9, uncounted, file=paste0("Sim_",Rep,"/SimData.Rdata"))
  save(dat9, uncounted, file=paste0("Sim_",Rep,"/SimData.Rdata"))
  
  # Turns ParValTabs items into initial value lists for chains, as well as True Parameters for posterior coverage
  initfxn <- function(PV){
    initlist <- list()
    params <- c("intcpt_a", "ba", "intcpt_d", "bd", "sigma_a", "sigma_d", "sigma_det", "alpha", "shape_k", "gamma")
    for (Pm in 1:length(params)){
      initlist[[Pm]] <- PV$X50.[grepl(paste("^",params[Pm],sep=""),PV$parameter)]
    }
    names(initlist) <- params
    initlist$sigma_d <- initlist$sigma_d[1]          # Removes unnecessary sigma_det
    initlist <- initlist[lapply(initlist,length)>0]  # Removes empty list items
    initlist$ra <- ra
    initlist$rd <- rd
    return(initlist)
  }
  
  ParValues <- vector('list', length(ParValTabs))
  for(j in 1:length(ParValTabs)) ParValues[[j]] <- initfxn(ParValTabs[[j]])
  
  save(ParValues, file=paste0("Sim_",Rep,"/TrueParams.Rdata"))
  
}

q(ifelse(interactive(), "ask", "no"))