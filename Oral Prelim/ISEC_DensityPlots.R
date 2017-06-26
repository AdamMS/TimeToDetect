# #### This document comes in 3 stages:
# # (i) Some code to identify the most representative Sim # across models... 
# #     - for gamma mix/no-mix, it said Sim_2 (based on 16 sims).  For 100 sims, the order was 17, 66, 96, 2, 10...
# #     - for expo & gamma data modeled by all families, it was Sim_14.  For 100 sims, the order was 14, 44, 73, 31, 35...
# # (ii) Code to extract p_det two different ways from the MCMC chains
# # (iii) Code to generate posterior plots for these representative data
# # Only the second stage should remain uncommented!
# 
#### Identify Representative Sims
library(plyr)

### First, we consider mixture/non-mixture gamma inference models
# tempS2 comes from 'Outputs_Code_zero.R'
# It consists only of same-family inference models!  So, it's uniquely identified by simdat and modelmix
smurf <- subset(tempS2, ModelFamily=="Gamma")

plist    <- c("X2.5.", "X25.", "X50.", "X75.", "X97.5.", "mean")
simdat   <- unique(smurf$simdat)
modelmix <- c(T,F)

# Calculate differences by rank
# Create a matrix with one row for each data type and one column for each Rep
# We will then sum across data types and choose the Rep with the optimal (closest to median) value
nReps <- length(unique(smurf$Rep))
Ranx  <- data.frame(param=NA, simdat=NA, modelmix=NA, matrix(NA, nrow=48, ncol=nReps))
indx  <- 0
for(p in 1:length(plist)){     # Statistic
  for(s in 1:length(simdat)){  # Dataset
    for(m in 0:1) {            # Mixture status of inference model
      indx <- indx+1
      Ranx[indx,1:3]  <- c(plist[p], simdat[s], m)
      Ranx[indx, 4:(nReps+3)] <- rank(subset(smurf, variable==plist[p] & simdat==simdat[s] & modelmix==m)$value)
    }
  }
}
fxn   <- function(x, nReps) abs(x-(nReps+1)/2)
Ranx2 <- apply(Ranx[,-(1:3)], 2, fxn, nReps=nReps)
sort(colSums(Ranx2))


### Next, we consider Family-comparison inference with expo/gamma models fit to expo/gamma data.  All mixtures.
# tempS3 also comes from 'Outputs_Code_zero.R'
# It consists only of mixture inference models.
smurf <- subset(tempS3, ModelFamily %in% c("Exponential","Gamma") & DataFamily %in% c("Exponential","Gamma"))

plist    <- c("X2.5.", "X25.", "X50.", "X75.", "X97.5.", "mean")
simdat   <- unique(smurf$simdat)
family   <- unique(smurf$ModelFamily)

nReps  <- length(unique(smurf$Rep))
combos <- length(plist) * length(simdat) * length(family) 
Ranx <- data.frame(param=NA, simdat=NA, family=NA, matrix(NA, nrow=combos, ncol=nReps))
indx <- 0
for(p in 1:length(plist)){
  for(s in 1:length(simdat)){
    for(f in 1:length(family)) {
      indx <- indx+1
      Ranx[indx,1:3] <- c(plist[p], simdat[s], family[f])
      Ranx[indx,4:(nReps+3)] <- rank(subset(smurf, variable==plist[p] & simdat==simdat[s] & ModelFamily==family[f])$value)
    }
  }
}
fxn   <- function(x, nReps) abs(x-(nReps+1)/2)
Ranx2 <- apply(Ranx[,-(1:3)], 2, fxn, nReps=nReps)
sort(colSums(Ranx2))





# #### Code to extract MCMC draws relevant to p_det
# # setwd("C:/Users/Adam/Documents/Adam/IAState/Stat 599JN/TimeToDetect/Sims/SimZero")
# library(rstan)
# 
# # Functions for each TTDD -- redundant variables (gamma, shape) make it easier to call later
# pfxn <- list(
#   function(time, intcpt_d, shape, gamma) pexp(time, exp(intcpt_d)),
#   function(time, intcpt_d, shape, gamma) gamma*pexp(time, exp(intcpt_d)) + (1-gamma)*(time>0),
#   function(time, intcpt_d, shape, gamma) pgamma(time, shape, exp(intcpt_d)),
#   function(time, intcpt_d, shape, gamma) gamma*pgamma(time, shape, exp(intcpt_d)) + (1-gamma)*(time>0),
#   function(time, intcpt_d, shape, gamma) plnorm(time, -intcpt_d, shape),
#   function(time, intcpt_d, shape, gamma) gamma*plnorm(time, -intcpt_d, shape) + (1-gamma)*(time>0),
#   function(time, intcpt_d, shape, gamma) pweibull(time, shape, exp(-intcpt_d)),
#   function(time, intcpt_d, shape, gamma) gamma*pweibull(time, shape, exp(-intcpt_d)) + (1-gamma)*(time>0)
# )
# 
# simtable    <- read.csv("S0_simtable.csv")
# fname       <- simtable$fname
# modeloutput <- as.character(simtable$modeloutput)
# mo          <- as.numeric(simtable$modeloutput)
# 
# # Which models to collect posterior distributions from
# set <- list(which(with(simtable, grepl("G", simdat) & parmcode %in% 3:4)),  # Sim Study 1 gamma-related
#             which(with(simtable, (grepl("ET", simdat) | grepl("GT", simdat)) & modelmix)))  # Sim Study 2 expo/gamma mixture data
# 
# indx <- 0
# for(L in 1:2){
#   if(L==1) setwd("./Sim_13")  # A representative replicate
#   if(L==2) setwd("./Sim_6")   # A representative replicate
#   for(i in set[[L]])  {
#     indx <- indx + 1
#     print(indx)
#     load(paste0(fname[i],".Rdata"))                     # File containing stanfit and other objects
#     stanobject <- eval(parse(text=modeloutput[i]))      # stanfit object
#     
#     intcpt_d <- rstan::extract(stanobject, "intcpt_d")$intcpt_d
#     if(simtable$modelmix[i]) {gamma <- rstan::extract(stanobject, "gamma")$gamma
#     } else {gamma <- NULL}
#     
#     if (mo[i] %in% 3:4){
#       shape <- rstan::extract(stanobject, "alpha")$alpha
#     } else if (mo[i] %in% 5:6){
#       shape <- rstan::extract(stanobject, "sigma_det")$sigma_det
#     } else if (mo[i] %in% 7:8){
#       shape <- rstan::extract(stanobject, "shape_k")$shape_k
#     }
#     
#     pdet <- pfxn[[mo[i]]](10, intcpt_d, shape, gamma)   # Calculate pdet from the posteriors of intcpt_d, gamma, shape
#     
#     uncounted <- rstan::extract(stanobject, "unobserved")$unobserved
#     totN <- rstan::extract(stanobject, "totN")$totN
#     
#     pdet.emp <- 1 - rowSums(uncounted) / rowSums(totN)  # Calculate pdet empirically from posterior of n_unc and N
#     # Posterior predictive
#     
#     ps <- list(pdet, pdet.emp)
#     
# #    save(ps, file=paste0("../pdets/pdetPost_", i, ".Rdata"))
#     
#     rm(list=c("stanobject", modeloutput[i]))
#   }
#   setwd("..")
# }
# 
# 
# #### Code to plot posterior densities for the representative sims
# # Note: I calculated p_det two different ways, but the two are nearly identical!  Yay!
# # Version 1: n / (n + uncounted)
# # Version 2: CDF(10, ...)
# 
# # Run the first several lines of the above code in order to generate the objects 'simtable' and 'set'
# library(ggplot2)
# setwd("../../Oral Prelim/pdetPost")
# 
# ptab <- simtable[,c("fname", "parmcode", "simdat", "modelmix", "datamix", "datamodelmatch")]
# 
# pdet <- NULL
# for(j in set[[1]]){
#   load(paste0("pdetPost_",j,".Rdata"))
#   p.add <- cbind(ptab[j,], pdet=ps[[1]], row.names=NULL)
#   pdet <- rbind(pdet, p.add)
# }
# pdet$peak <- substr(pdet$simdat,3,3)
# 
# pdet$peak <- factor(pdet$peak, labels=c("Shape < 1", "Shape > 1"))
# pdet$datamix <- factor(pdet$datamix, labels=c("Nonmixture Data", "Mixture Data"))
# pdet$modelmix <- factor(pdet$modelmix, labels=c("Nonmixture", "Mixture"))
# pdet <- pdet[c(F,F,T),]    # Thinning
# 
# # Note: axis.line=element_line(color="black") didn't work, so I had to manually add
# # png("../Posteriors_Mixtures.png", height=5.75, width=8.5, units="in", res=96)
# print(
#   ggplot(pdet, aes(pdet, group=modelmix, color=modelmix)) + 
#   geom_hline(yintercept=0) + geom_vline(xintercept=0) +
#   stat_density(adjust=1, geom="line", position="identity", size=1.5) +
#   scale_color_manual(name="Inference Model:", values=c("palegreen4", "palegreen3")) +
#   geom_segment(aes(x=0.8, xend=0.8, y=0, yend=15), linetype=2, size=1.5, color="maroon2") +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position="bottom",
#         strip.text = element_text(size = 12),
#         legend.title = element_text(size = 12)) +
#   xlab("Posterior Detection Probability") + ylab("Density") +
#   facet_grid(datamix~peak) 
# )
# # dev.off()
# 
# 
# 
# set[[2]] <- c(5,12,13,40,46,47)   # Restrict to just expo and gamma models
# pdet.fam <- NULL
# for(j in set[[2]]){
#   load(paste0("pdetPost_",j,".Rdata"))
#   p.add <- cbind(ptab[j,], pdet=ps[[1]], row.names=NULL)
#   pdet.fam <- rbind(pdet.fam, p.add)
# }
# 
# pdet.fam$peak <- substr(pdet.fam$simdat,3,3)
# 
# pdet.fam$peak <- factor(pdet.fam$peak, labels=c("Shape < 1", "Shape > 1"))
# # pdet.fam$simdat <- factor(pdet.fam$simdat, 
# #                           levels=c("GTF", "ETF", "GTT"),
# #                           labels=c("Nonconstant \n Gamma; Shape < 1", "Constant \n Exponential", "Nonconstant \n Gamma; Shape > 1"))
# pdet.fam$simdat <- factor(pdet.fam$simdat, 
#                           levels=c("GTF", "ETF", "GTT"),
#                           labels=c("Gamma; Shape < 1", "Exponential", "Gamma; Shape > 1"))
# pdet.fam$parmcode <- factor(pdet.fam$parmcode, labels=c("Exponential", "Gamma"))
# pdet.fam <- pdet.fam[c(F,F,T),]    # Thinning
# #pdet.fam <- pdet.fam[c(F,F,F,F,F,F,F,F,F,T),]
# 
# # png("../Posteriors_Families_draft.png", height=3.5, width=7, units="in", res=96)
# print(
#   ggplot(pdet.fam, aes(pdet, group=parmcode, color=parmcode)) + 
#     geom_hline(yintercept=0) + geom_vline(xintercept=0) +
#     stat_density(adjust=1, geom="line", position="identity", size=1.5) +
#     scale_color_manual(name="Inference Model Family:", values=c("gray70", "gray40")) +
#     geom_segment(aes(x=0.8, xend=0.8, y=0, yend=12), linetype=2, size=1.5, color="black") +
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           legend.position="bottom",
#           strip.text = element_text(size = 11),
#           legend.title = element_text(size = 12),
#           legend.text=element_text(size=12),
#           axis.title=element_text(size=12),
#           axis.text.x=element_text(size=11)) +
#     xlab("Posterior Detection Probability") + ylab("Density") + ylim(0,12.2) + xlim(-0.05,1.05) +
#     facet_grid(~simdat)
# )
# # dev.off()
# 
# png("../Posteriors_Families_ASA.png", height=6.25, width=9, units="in", res=96)
# print(
#   ggplot(pdet.fam, aes(pdet, group=parmcode, color=parmcode)) + 
#     geom_hline(yintercept=0) + geom_vline(xintercept=0) +
#     stat_density(adjust=1, geom="line", position="identity", size=1.5) +
#     scale_color_manual(name="Inference Model Family:", values=c("palegreen4", "palegreen3")) +
#     geom_segment(aes(x=0.8, xend=0.8, y=0, yend=12), linetype=2, size=1.5, color="maroon2") +
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           legend.position="bottom",
#           strip.text = element_text(size = 20),
#           legend.title = element_text(size = 20),
#           legend.text=element_text(size=20),
#           axis.title=element_text(size=20),
#           plot.title=element_text(size=34),
#           axis.text.x=element_text(size=11)) +
#     xlab("Posterior Detection Probability") + ylab("Density") + ylim(0,12.2) +
#     facet_grid(~simdat) + ggtitle("Posteriors for Pr(detection)")
# )
# dev.off()

