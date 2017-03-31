########## Code to fit SimZero models on slurm
set.seed(NULL)
bundle   <- 1    # How many models will be fit by a single call to this script

# Code to replicate S0_sumtable according to the number of replicates
NumReps  <- 1
simtable <- read.csv("SF_simtable.csv")       # Load table of datasets
simtable <- do.call('rbind', rep(list(simtable), times=NumReps))
simtable$Rep <- rep(1:NumReps, each=nrow(simtable)/NumReps)

# Rows of simtable that I want to (re-)run
indices  <- 1:nrow(simtable)   # Run the full simulation, all replicates
# Repeat models that have effective sample size less than 1000 or NA (didn't fit before)
# ess      <- readRDS("ESS.rds")
# indices  <- which(is.na(ess) | ess<1000)

# Collect array number from slurm and add one.  
# 'j' identifies which models will be fit based on array number and bundle number
ss = Sys.getenv("SLURM_ARRAY_TASK_ID") 
j  = as.numeric(ss) * bundle + 1

for(indx in indices[j:(j+bundle-1)]){
  Rep <- simtable$Rep[indx]
  load(paste0("Sim_", Rep, "/SimData.Rdata"))  # Load all datasets in Rep.  Object name is 'dat9'
  dataset <- dat9[[simtable$datacode[indx]]]
  load(paste0("Sim_",Rep,"/TrueParams.Rdata")) # Load True Values associated with SimData.  Object name is 'ParValues'
  
  # Provide initial values.  Since data and models are often not congruent, the idea here
  # is to take simulation parameters from the most equivalent model-related family.
  # For example, if the data are LFT, and the model is a GammaMix, then I want the 
  # parameters from the GTT model
  modelcode <- paste0(substr(simtable$modelname,1,1), substr(simtable$modelmix,1,1), substr(simtable$simdat,3,3))
  modelcode[modelcode=="EFT"] <- "EFF"
  modelcode[modelcode=="ETT"] <- "ETF"
  initlist <- ParValues[[which(modelcode[indx]==levels(simtable$simdat))]]
  
  # For model re-runs only... set the initlist based on the last not-converged / inadequate ESS fit:
  if(file.exists(paste0("Sim_1","/",simtable$fname[indx],"_sum.Rdata"))){
    load(paste0("Sim_1","/",simtable$fname[indx],"_sum.Rdata"))    # Load Previous Model Fit Summary
    for (j in 1:(length(initlist)-2)){   # The minus-2 is for the random effects
      initlist[[j]] <- outtable[[1]][,"50%"][grep(names(initlist)[j], names(outtable[[1]][,"50%"]))]
      if(names(initlist)[j]=="sigma_d") initlist[[j]] <- initlist[[j]][1]  # Because of sigma_d / sigma_det naming similarity
    }
  }
  initlist$ra <- c(rnorm(dataset$n_ras[1], 0, initlist$sigma_a[1]), rnorm(dataset$n_ras[2], 0, initlist$sigma_a[2]))
  initlist$rd <- rnorm(dataset$n_rds, 0, initlist$sigma_a[1])
  
  ##### MCMC variables to store
  # record.list <- c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","unobserved","yrep","ra","rd","lambda","totN","dev1")
  record.list <- c("intcpt_a","ba","intcpt_d","bd","sigma_a","sigma_d","uncounted","abundance","p_global","dev1")
  if(simtable$modelmix[indx]) record.list <- c(record.list, "gamma")
  if(simtable$parmcode[indx] %in% 1:2) record.list <- c(record.list)              # "rho"
  if(simtable$parmcode[indx] %in% 3:4) record.list <- c(record.list, "alpha")     # "beta"
  if(simtable$parmcode[indx] %in% 5:6) record.list <- c(record.list, "sigma_det") # "mu_det"
  if(simtable$parmcode[indx] %in% 7:8) record.list <- c(record.list, "shape_k")   # "invrate"
  
  library(rstan)
  library(mcmcse)
  m = stan_model(paste0("../../OVEN/", simtable$modelname[indx],".stan"), auto_write=rstan_options("auto_write" = TRUE))
  
  # For most runs, it was 60K iter.  For a few troublesome models, I boosted it.
  insuff <- TRUE
  iter   <- 90000
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
}