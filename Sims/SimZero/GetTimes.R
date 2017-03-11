##### Gets Run Times
# Requires that StanPlot.R has already been run (for Geweke statistics)
library(plyr)
library(rstan)
library(mcmcse)
simtable <- read.csv("S0_simtable.csv")
nReps    <- 100
ess      <- runtimes <- data.frame(matrix(NA, nrow=nrow(simtable), ncol=nReps))
names(ess)      <- paste0("ESS_", 1:nReps)
names(runtimes) <- paste0("RunTime_", 1:nReps)
Gew.compile <- NULL
Gfxn <- function(listitem, Rep) data.frame(parameter=names(listitem), zScore=listitem, Rep=Rep)
for(Rep in 1:nReps){
  for(m in 1:nrow(simtable)){
    loadfile <- paste0("Sim_", Rep, "/", simtable[m,"fname"], ".Rdata")
    if(file.exists(loadfile)){
      load(loadfile)
      runtimes[m, Rep] <- sum(get_elapsed_time(eval(parse(text=as.character(simtable[m,"modeloutput"])))))  # Total warmup plus sampling time
      postdraws        <- as.data.frame(eval(parse(text=as.character(simtable[m,"modeloutput"]))))          # Posterior draws
      ess[m, Rep]      <- min(ess(postdraws))    # Effective Sample Size
    } else {
      runtimes[m,Rep]  <- NA
      ess[m, Rep]      <- NA
    }
  }
  if(file.exists(paste0("Sim_",Rep,"/GewekeDiags.Rdata"))){
    load(paste0("Sim_",Rep,"/GewekeDiags.Rdata"))  # Loads list 'Geweke' containing Geweke statistics for each parameter from all models in Rep
    names(Geweke) <- simtable$fname                # Label with model name (e.g. 'ExpoMix_LTF')
    Gew.compile   <- rbind(Gew.compile, ldply(Geweke, Gfxn, Rep, .id="model")) # Compile Geweke diagnostics in a single table for all Reps.
      # Long-format data frame.
  }
}
saveRDS(runtimes, file="RunTimes.rds")
saveRDS(ess, file="ESS.rds")
if(!is.null(Gew.compile)) saveRDS(Gew.compile, file="Geweke.rds")
