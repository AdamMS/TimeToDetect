########## Code to fit SimZero models on slurm
set.seed(NULL)
bundle   <- 1    # How many models will be fit by a single call to this script

# Code to replicate S0_sumtable according to the number of replicates
NumReps  <- 100
simtable <- read.csv("S0_simtable.csv")       # Load table of datasets
simtable <- do.call('rbind', rep(list(simtable), times=NumReps))
simtable$Rep <- rep(1:NumReps, each=nrow(simtable)/NumReps)

# Rows of simtable that I want to (re-)run
# indices  <- 1:nrow(simtable)   # Run the full simulation, all replicates
# Repeat models that have effective sample size less than 1000 or NA (didn't fit before)
ess      <- readRDS("ESS.rds")
indices  <- which(is.na(ess) | ess<1000)

# Collect array number from slurm and add one.  
# 'j' identifies which models will be fit based on array number and bundle number
ss = Sys.getenv("SLURM_ARRAY_TASK_ID") 
j  = as.numeric(ss) * bundle + 1

for(indx in indices[j:(j+bundle-1)]){
  Rep <- simtable$Rep[indx]
  load(paste0("Sim_", Rep, "/SimData.Rdata"))  # Load all datasets in Rep.  Object name is 'dat9'
  dataset <- dat9[[simtable$datacode[indx]]]
  load(paste0("Sim_",Rep,"/TrueParams.Rdata")) # Load True Values associated with SimData.  Object name is 'ParValues'
  
  # Provide initial values.  Since data and models are often not congruent, but the 
  # simulations use the same target values, the idea here is to use simulation parameters
  # from the most equivalent model-related family.
  # For example, if the data are LFT, and the model is a GammaMix, then I want the 
  # parameters from the GTT model
  
  # Construct a simdat-like code for the model being fit so that I can use values from the equivalent dataset
  modelcode <- paste0(substr(simtable$modelname,1,1), substr(simtable$modelmix,1,1), substr(simtable$simdat,3,3))
  modelcode[modelcode=="EFT"] <- "EFF"
  modelcode[modelcode=="ETT"] <- "ETF"
  initlist <- ParValues[[which(modelcode[indx]==levels(simtable$simdat))]] # Initial values
  
  # MCMC variables to store
  record.list <- c("intcpt_a", "intcpt_d", "uncounted", "p_global", "dev1")
  if(simtable$modelmix[indx]) record.list <- c(record.list, "gamma")
  if(simtable$parmcode[indx] %in% 3:4) record.list <- c(record.list, "alpha")
  if(simtable$parmcode[indx] %in% 5:6) record.list <- c(record.list, "sigma_det")
  if(simtable$parmcode[indx] %in% 7:8) record.list <- c(record.list, "shape_k")
  
  library(rstan)
  library(mcmcse)
  # Compile/Load Stan model
  m = stan_model(paste0(simtable$modelname[indx],".stan"), auto_write=rstan_options("auto_write" = TRUE))
  
  # For most runs, it was 60K iter.  For a few troublesome models, I boosted it.
  insuff <- TRUE
  iter   <- 60000
  fail   <- 0
  while(insuff){
    assign(as.character(simtable$modeloutput[indx]), 
           sampling(m, data=dataset, chains=1, iter=iter*(2^fail), thin=10, pars=record.list, 
                    init=list(initlist)))
    postdraws <- as.data.frame(eval(parse(text=as.character(simtable$modeloutput[indx]))))
    if(min(ess(postdraws)) > 1000){
      insuff=FALSE
    } else {
      fail <- fail + 1
    }
  }
  
  save(list=as.character(simtable$modeloutput[indx]), file=paste0("Sim_",Rep,"/",simtable$fname[indx],".Rdata"))
  
  rm(list = ls()[-which(ls() %in% c("indx", "j", "bundle", "simtable"))])
}

q(ifelse(interactive(), "ask", "no"))
