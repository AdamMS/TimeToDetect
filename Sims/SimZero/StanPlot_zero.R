########## This misnamed code actually does all the immediate post-model processing of data and
# returns a list object called 'outtable' for each model fit.  Outtable is saved in the '<modelname>_sum.Rdata' file.


########## Load Stan output
library(mcmcse)
library(coda)
library(rstan)
library(hexbin)

simtable <- read.csv("S0_simtable.csv")   # Simulation information is stored here

# Simplify calls to simtable
fname       <- simtable$fname
modeloutput <- as.character(simtable$modeloutput)
parmcode    <- simtable$parmcode
n_ints      <- simtable$n_ints
datacode    <- simtable$datacode
modelmix    <- simtable$modelmix
datamodelmatch <- simtable$datamodelmatch
params <- list(c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","lp__"),                     # Exponential
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","gamma","lp__"),             # ExpoMix
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","alpha","lp__"),             # Gamma
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","gamma","alpha","lp__"),     # GammaMix
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","sigma_det","lp__"),         # LogN
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","gamma","sigma_det","lp__"), # LogNMix
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","shape_k","lp__"),           # Weibull
               c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","gamma","shape_k","lp__"))   # WeibullMix
# Because it's the intercept-only model
for(l in 1:length(params)) params[[l]] <- params[[l]][-which(params[[l]] %in% c("ba", "bd", "sigma_a", "sigma_d", "lp__"))]

Reps <- 100

for(Rep in 1:Reps){
  ESS <- NULL; Geweke <- vector('list', nrow(simtable))  # Storage objects
  setwd(paste0("Sim_",Rep))
  load("SimData.Rdata")                            # Load data (dat9) for all simulations in Rep.
                                                   # Also includes the total number of birds uncounted (uncounted)
  for(i in 1:nrow(simtable)){ 
    load(paste0(fname[i],".Rdata"))                # Load Model Output
    stanobject <- eval(parse(text=modeloutput[i])) # Call the model output 'stanobject'
    
    outtable <- vector('list',5)
    
    # outtable[[1]]: rstan summary table
    outtable[[1]] <- summary(stanobject, pars=params[[parmcode[i]]])$summary
    
    # outtable[[2]]: uncounted birds and overall detection probability
    outtable[[2]] <- summary(stanobject, pars=c("uncounted", "p_global"))$summary
    
    # Read in correct simulated data
    dat <- eval(parse(text=paste("dat",n_ints[i],sep="")))[[datacode[i]]] 
    
    ### There used to be code here to run posterior predictive checks, but we have abandoned that in the large-scale simulation study
    # Posterior predictive p-value table for counts by interval
    outtable[[3]] <- NULL
    
    ##### Code to calculate posterior p-values for model parameters and quantitities
    # Note: 'p80' is the posterior p-value relative to the true parameter values.  We want this one.
    # 'p_global' is the posterior p-value relative to the realized detection probability from the dataset.
    # ALL P-VALUES ARE Pr(ACTUAL < SIMULATED)
    load("../Simpars.Rdata")   # True parameter values used in simulations
    Pval.param <- rep(NA,(nrow(outtable[[1]])+3))
    # Posterior p-values for total abundance and global detection probability
    unc.draws <- rstan::extract(stanobject, "uncounted")$uncounted  # Posterior samples of total uncounted
    counted   <- sum(dat$y)                                         # Total count from simulated dataset
    tot.unc   <- sum(uncounted[[datacode[i]]])                      # Total uncounted from simulated dataset
    p_global  <- rstan::extract(stanobject, "p_global")$p_global    # Posterior samples of detection probability
    NumDraws <- length(unc.draws)                                   # Number of posterior draws
    # Posterior Pr( E[Abundance] < posterior abundance ), where E[...] calculated from simulation parameters
    Pval.param[1] <- sum( nrow(dat$y)*exp(Simpars$intcpt_a[datacode[i]]) < counted + unc.draws ) / NumDraws
    # Posterior Pr( realized pdet from simulation < posterior p_global )
    Pval.param[2] <- sum( counted/(counted+tot.unc) < p_global ) / NumDraws
    # Posterior Pr( E[pdet] < posterior p_global ), where E[...] calculated from simulation parameters
    Pval.param[3] <- sum( Simpars$pdet[datacode[i]] < p_global ) / NumDraws
    names(Pval.param) <- c("Uncounted", "p_global", "p80")
    
    if(datamodelmatch[i]){
      ### Posterior p-values for all model parameters noted in 'params' up top
      load("TrueParams.Rdata")
      param.extract <- rstan::extract(stanobject, pars=params[[parmcode[i]]]) # Posterior draws
      counter <- 4
      for (j in 1:length(param.extract)) {
        # The if() statement is TRUE for one-dimensional parameters (e.g. 'intcpt_a') but
        # is FALSE for multi-dimensional parameters (e.g. 'ba' when there are multiple predictors)
        if(is.na(ncol(param.extract[[j]]))) {ncolsj <- 1} else {ncolsj <- ncol(param.extract[[j]])}
        for (k in 1:ncolsj){
          # Look up the true value of the parameter in question
          TrueValue <- ParValues[datacode[i]][[1]][[match(
            params[[parmcode[i]]][j], names(ParValues[datacode[i]][[1]]))]][k]
          if(ncolsj > 1) {
            Pval.param[counter] <- sum(TrueValue < param.extract[[j]][,k])/NumDraws
          } else {
            Pval.param[counter] <- sum(TrueValue < param.extract[[j]])/NumDraws}
          counter <- counter + 1
        }
      }
      names(Pval.param) <- c("Uncounted", "p_global", "p80", rownames(outtable[[1]]))
    }
    outtable[[4]] <- Pval.param
    
    ##### Calculate Geweke diagnostics
    # Extract non-permuted posterior draws
    gew.ext2 <- rstan::extract(stanobject, 
                               pars=params[[parmcode[i]]], permuted=F, inc_warmup=F)[,1,]
    # Calculate Geweke diagnostic for each parameter
    Geweke[[i]] <- geweke.diag(gew.ext2)[[1]]
    
    ##### Effective Sample Sizes
    ESS <- rbind(ESS, 
                 data.frame(parameter = rownames(outtable[[1]]),
                      n_eff = ess(gew.ext2),  # Using the 'mcmcse' package
                      fname = fname[i],
                      Rep   = paste0("Rep_", Rep))
                 )
    
    ########## DIC -- there used to be three versions of this.  But I'm only keeping the first.
    # dev1: p(n|betas, gamma, shape)
    DIC.varmeth <- function(dev) mean(dev) + 0.5*var(dev)
    DIC.minmeth <- function(dev) 2*mean(dev) - min(dev)
    DICfx <- function(dev) return(c(Var=DIC.varmeth(dev), Min=DIC.minmeth(dev)))
    
    # Deviance from each iteration of the MCMC chain 
    dev1 <- rstan::extract(stanobject, "dev1")[[1]]
    DIC <- data.frame(dev1=DICfx(dev1))
    
    meanind <- which(dimnames(outtable[[1]])[[2]]=="mean") # Which col is 'mean'
    midind <- which(dimnames(outtable[[1]])[[2]]=="50%") # Which col is '50%'
    
    ##### Calculate Dev at the median value --- this seems the most reliable approach to me, because N's are discrete
    lpn1 <- lambda <- phi <- numeric(nrow(dat$y))
    pvec <- numeric(ncol(dat$y))
    tau <- c(0,2:10)
    gamma <- shape <- NULL
    
    # Median values
    intcpt_a <- outtable[[1]][rownames(outtable[[1]])=="intcpt_a", midind]
    intcpt_d <- outtable[[1]][rownames(outtable[[1]])=="intcpt_d", midind]
    if(modelmix[i]) gamma <- outtable[[1]][rownames(outtable[[1]])=="gamma", midind]
    shape <- outtable[[1]][rownames(outtable[[1]]) %in% c("alpha", "sigma_det", "shape_k"), midind]
    unobs.med <- floor(median(unc.draws)) # The 'floor' is for the rare case when the median is non-integer
    totN.med  <- unobs.med + counted
    
    ### Calculate D(\theta_{median})
    # Functions for each TTDD -- redundant variables (gamma, shape) make it easier to call later
    # Note: I wrote these strangely... pfxn = p.det(hard-to-detect) - (1-gamma)I(t==0)
    #  but it's okay, since the function is only used for pdet differences for calculating interval-specific probs
    if (modeloutput[i]=="f.exp9") pfxn <- function(time, intcpt_d, shape, gamma) pexp(time, exp(intcpt_d))
    if (modeloutput[i]=="f.expmix9") pfxn <- function(time, intcpt_d, shape, gamma)
      gamma*pexp(time, exp(intcpt_d)) - (1-gamma)*(time==0)
    if (modeloutput[i]=="f.gamma9") pfxn <- function(time, intcpt_d, shape, gamma) pgamma(time, shape, exp(intcpt_d))
    if (modeloutput[i]=="f.gammamix9") pfxn <- function(time, intcpt_d, shape, gamma)
      gamma*pgamma(time, shape, exp(intcpt_d)) - (1-gamma)*(time==0)
    if (modeloutput[i]=="f.logN9") pfxn <- function(time, intcpt_d, shape, gamma) plnorm(time, -intcpt_d, shape)
    if (modeloutput[i]=="f.logNMix9") pfxn <- function(time, intcpt_d, shape, gamma)
      gamma*plnorm(time, -intcpt_d, shape) - (1-gamma)*(time==0)
    if (modeloutput[i]=="f.weibull9") pfxn <- function(time, intcpt_d, shape, gamma) pweibull(time, shape, exp(-intcpt_d))
    if (modeloutput[i]=="f.weibullmix9") pfxn <- function(time, intcpt_d, shape, gamma)
      gamma*pweibull(time, shape, exp(-intcpt_d)) - (1-gamma)*(time==0)
    
    # Calculate log(p(n|...))
    lambda <- exp(intcpt_a)   # + X_i %*% \beta + Z_i %*% \ksi      # i.e., expected abundance at survey j
    # Note: the following assumes an intercept-only model for detection, so all sites are identical
    for(intvl in 1:ncol(dat$y)) pvec[intvl] <- pfxn(tau[intvl+1],intcpt_d,shape,gamma) - pfxn(tau[intvl],intcpt_d,shape,gamma)
    lpn1 <- sum(dpois(colSums(dat$y), nrow(dat$y)*lambda*pvec, log=T))                     # log(p(n|int_a,int_d))

    # DEVIANCE CALCULATED DIRECTLY FROM \HAT{THETA}:
    # pD := effective number of parameters
    D.Theta1 <- -2 * sum(lpn1); pD1 <- mean(dev1)-D.Theta1 
    # DIC AS CALCULATED ABOVE:
    DIC.med1 <- 2 * mean(dev1) - D.Theta1
    
    DIC <- rbind(DIC, ThetaMed=c(DIC.med1))
    pD <- t(t(DIC) - c(mean(dev1)))
    outtable[[5]] <- list(DIC, pD)
    
    ### Save output!
    save(outtable, file=paste(fname[i],"_sum.Rdata",sep=""))
    
    
    rm(list=ls()[-match(c("i","Rep","dat9","fname","params","modeloutput","simtable","n_ints","parmcode","Geweke","ESS","datacode","modelmix","datamodelmatch","uncounted"),ls())])
    gc(verbose=T)
  }
  
  # save(quanttable, file="quanttable_oven.Rdata")
  save(Geweke, file="GewekeDiags.Rdata")
  save(ESS, file="ESS.Rdata")
  setwd("../")
}