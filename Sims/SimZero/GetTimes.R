##### Gets Run Times
library(rstan)
library(mcmcse)
simtable <- read.csv("S0_simtable.csv")
nReps    <- 100
ess <- runtimes <- data.frame(matrix(NA, nrow=nrow(simtable), ncol=nReps))
names(ess)      <- paste0("ESS_", 1:nReps)
names(runtimes) <- paste0("RunTime_", 1:nReps)
for(Rep in 1:nReps){
  for(m in 1:nrow(simtable)){
    loadfile <- paste0("Sim_", Rep, "/", simtable[m,"fname"], ".Rdata")
    if(file.exists(loadfile)){
      load(loadfile)
      runtimes[m, Rep] <- sum(get_elapsed_time(eval(parse(text=as.character(simtable[m,"modeloutput"])))))
      postdraws        <- as.data.frame(eval(parse(text=as.character(simtable[m,"modeloutput"]))))
      ess[m, Rep]      <- min(ess(postdraws))
    } else {
      runtimes[m,Rep]  <- NA
    }
  }
}
saveRDS(runtimes, file="RunTimes.rds")

# ##### Extracts Run Times
# library(plyr)
# library(reshape2)
# simtable <- read.csv("S0_simtable.csv")
# runtimes <- readRDS("RunTimes.rds")
# runtimes <- cbind(simtable[,1:2], runtimes)
# r.melt <- melt(runtimes, id.vars=1:2, value.name="runtime")
# smurf <- ddply(r.melt, .(modelname), summarize, mu.t=mean(runtime))
# smurf$minutes <- smurf$mu.t/60
# smurf
