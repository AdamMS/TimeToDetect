########## Code to gather results from across all simulations
library(reshape2)
library(ggplot2)
library(plyr)

nReps <- 100

##### Load pre-saved summary tables
simtable    <- read.csv("S0_simtable.csv")
fname       <- simtable$fname
modeloutput <- simtable$modeloutput
parmcode    <- simtable$parmcode
n_ints      <- simtable$n_ints
datacode    <- simtable$datacode
modelmix    <- simtable$modelmix
datamodelmatch <- simtable$datamodelmatch
simdat <- as.character(simtable$simdat)
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

parmests <- count.det <- posteriordist <- dic <- true.unc <- true.pglobal <- post.unc.p <- NULL

for(Rep in 1:nReps){
  setwd(paste0("Sim_",Rep))
  
  # For each of the 14 simulated datasets, collect the true number of uncounted individuals and 
  # the true global_p (proportion of simulated birds that are detected)
  load("SimData.Rdata")
  true.unc <- rbind(true.unc, data.frame(Rep=Rep, datacode=1:14, value=laply(uncounted, sum)))
  globhelp <- function(L) sum(L$y)
  true.pglobal <- rbind(
    true.pglobal, data.frame(
      Rep=Rep, datacode=1:14, value=laply(dat9, globhelp)/(laply(uncounted, sum)+laply(dat9, globhelp))
    )
  )
  
  ##### Compile output from the various models into data frames
  # A helper function to process a dataframe based on an index 'i' and the simtable
  df.process <- function(df, i, simtable){
    df$Rep <- Rep; df$indx <- i; df$parameter <- rownames(df)
    df <- cbind(df, simtable[i,])
    return(df)
  }
  
  for(i in 1:length(fname)){
    load(paste(fname[i],"_sum.Rdata",sep=""))
    # parmests is all of the summary tables appended together plus their accompanying simtable rows
    parmests <- rbind(parmests, df.process(data.frame(outtable[[1]]), i, simtable))
    # count.det is the 'uncounted' and 'p_global' summaries appended together
    count.det <- rbind(count.det, df.process(data.frame(outtable[[2]]), i, simtable))
    # In 2015, there was code here to process posterior predictive p-values (not to be confused with posterior p-values).
    # I am no longer generating that data from the no-intercept sim study.

    # Posterior p-value for parameter estimates and uncounted/p_global (NOT a posterior predictive p-value)
    # Note: using long format rather than wide, since parameters change among models
    # For parameters, we only calculate posterior p-values when the inference model matches data model
    if(datamodelmatch[i]) posteriordist <- rbind(posteriordist, cbind(pval=outtable[[4]], 
                                          Rep=Rep, indx=i, simtable[i,]))
    post.unc.p <- rbind(post.unc.p, df.process(data.frame(pval=outtable[[4]][1:3]), i, simtable))
    
    # DIC Tables
    dic <- rbind(dic, df.process(outtable[[5]][[1]], i, simtable))
  }
    setwd("../")
}

# Note: in the past, I have saved parmests to file for use in simulating datasets.
outputs <- list(parmests=parmests, count.det=count.det, posteriordist=posteriordist, 
                post.unc.p=post.unc.p, dic=dic, true.pglobal=true.pglobal)
save(outputs, file="Outputs.Rdata")

