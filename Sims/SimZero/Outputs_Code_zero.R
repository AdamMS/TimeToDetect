# This code generates all tables and plots for the manuscript
# As input, it needs: S0_simtable.csv, "Outputs_p##.Rdata", "Simpars_p##.Rdata"

# Specify which pdet we're analyzing
truepdet <- 0.5
fileID   <- "50"

library(ggplot2)
library(reshape2)
library(plyr)
last <- function(x) substring(x, nchar(as.character(x)))  # Function to extract last character from a string
# Lookup tables
familytab <- data.frame(Letter=c("E","G","L","W"), ModelFamily=c("Exponential", "Gamma", "Lognormal", "Weibull"))
familytab2 <- data.frame(Letter2=c("E","G","L","W"), DataFamily=c("Exponential", "Gamma", "Lognormal", "Weibull"))

data.handle <- function(df){
  df$Letter  <- substr(df$fname,1,1)   # Family of model
  df$Letter2 <- substr(df$simdat,1,1)  # Family of data
  df <- join(df, familytab, type="left", by="Letter")
  df <- join(df, familytab2, type="left", by="Letter2")
  df$familymatch <- with(df, Letter==Letter2) # Is model from the correct family (regardless of mixture)?
  return(df)
}

simtable <- read.csv("S0_simtable.csv")    # Load simtable
load(paste0("Outputs", fileID, ".Rdata"))  # Load model fits
for(i in 1:length(outputs)) assign(names(outputs)[i], outputs[[i]]) # Convert model fits from list into objects

##### There used to be a section here for posterior predictive p-values
# See First_simulations folder for code

###### Posterior p-values for global estimates of pdet
# There's a choice here. 
# Choice one: use the p-values associated with "p_global"... posterior Pr(pdet for simulated dataset < counted/total abundance)
# Choice two: use the p-values associated with "pTruth"...   posterior Pr(True pdet < counted/total abundance)
# The former is actually more flexible, because the latter is difficult to calculate across random effects.
# At large sample size, there is no practical difference between the two.
p_glob <- subset(post.unc.p, parameter=="p_global")

# Lots of data handling
names(p_glob)[1]     <- "pval_pglob"
p_glob               <- data.handle(p_glob)
p_glob$datapeaked    <- as.logical(last(p_glob$simdat))                       # Is data peaked
p_glob$datamixlabel  <- "Nomix"                                               # Label whether data mixture
p_glob$datamixlabel[p_glob$datamix] <- "Mix"
p_glob$datapeaklabel <- "Nonpeaked"
p_glob$datapeaklabel[p_glob$datapeaked] <- "Peaked"
p_glob$datapeaklabel[substr(p_glob$simdat,1,1)=="E"] <- "Peaked"              # Reckon Expos among the peaked distributions
p_glob$MixPeakLabel  <- with(p_glob, paste(datamixlabel, datapeaklabel))      # Label for mixture/peakedness

# Posterior estimates of detection probabilities -- By mix/non-mix and by family
CD.melt <- melt(subset(count.det, parameter=="p_global"), id.vars=11:ncol(count.det)) # Exclude 'uncounted' estimates
CD.melt <- CD.melt[,-which(names(CD.melt) %in% c("indx"))]
pglobs  <- join(true.pglobal, simtable)  # Construct a table of "true pdet values" with same variables as posterior estimates
                                         # These pdet values are realizations from each simulated dataset
pglobs$parameter <- "p_global"
pglobs$variable  <- "True"
CD.melt <- rbind(CD.melt, pglobs)
CD.melt <- data.handle(CD.melt)
# Label the family/peakedness of the data
CD.melt$familypeak <- "Nonpeaked"
CD.melt$familypeak[last(CD.melt$simdat)=="T"] <- "Peaked"
CD.melt$familypeak <- with(CD.melt, paste0(substr(DataFamily,1,4), " ", familypeak))
# Label whether data and model are mixtures
CD.melt$dmix.mmix <- paste0("Dat: ",substr(CD.melt$datamix,1,1), 
                            " Mod: ",substr(CD.melt$modelmix,1,1))

##### Testing whether mixtures are needed
### Quick overview histograms of Posterior p-values for each model fit.  All models are from the correct family.
# Correct model is 'blue'.  Incorrect model is 'red'.
# Facet rows    -- data model
# Facet columns -- mixture/peakedness of data
tempS2 <- subset(CD.melt, familymatch)
# pdf("pdet_post_correct.pdf", width=7.75, height=4.5)
ggplot(p_glob[p_glob$familymatch,]) + 
  geom_histogram(aes(pval_pglob, fill=datamodelmatch), binwidth=0.05) +
  facet_grid(ModelFamily+modelmix ~ MixPeakLabel, scale="free_y") + theme_bw()
# dev.off()

### Caterpillar plots for EVERY rep (quite a lot).  All models are from the correct family.
# Facet rows    -- mixture status of data and model
# Facet columns -- Family and peakedness of data 
# Background shading trick
tp <- unique(tempS2[,c("dmix.mmix","familypeak","datamodelmatch")])
# pdf(paste0("Misc/Cater",fileID,"mixtures.pdf"), width=9.25, height=6.1)
print(ggplot(data=subset(tempS2, variable %in% c("X2.5.","X97.5."))) + 
        geom_line(aes(x=value, y=as.factor(Rep), group=Rep)) + 
        geom_rect(aes(fill = datamodelmatch), data=tp,
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_fill_manual(values = c("red", "gray90")) + 
        facet_grid(dmix.mmix~familypeak) + ylab("Simulation Replicate") + xlab("Probability") + 
        geom_line(aes(x=value, y=as.factor(Rep), group=Rep), color="orange", size=I(0.6), 
                  data=subset(tempS2,variable %in% c("X25.","X75."))) + 
        geom_point(aes(x=value, y=as.factor(Rep)), color="blue", size=I(0.85), data=subset(tempS2,variable %in% "X50."))  +
        geom_vline(xintercept=truepdet, linetype=1, color="white", size=I(1.0)) +
        scale_shape_manual(values=c(2)) + theme_bw() +
        theme(axis.text.x=element_text(angle=45), axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none")
      )
# dev.off()

# There used to be a caterpillar plot of just lognormal data/models to provide an example for the Oral Prelim



##### Testing family misspecification
### Quick overview histograms of Posterior p-values for each model fit.  All models are from the correct family.
# Correct model is 'blue'.  Incorrect model is 'red'.
# Facet rows    -- data model
# Facet columns -- family and peakedness of data ('Mix' label is superfluous but correct)
# pdf("pdet_post_family.pdf", width=7.75, height=4.5)
ggplot(data=p_glob[p_glob$modelmix & p_glob$datamix,]) + 
  geom_histogram(aes(pval_pglob, fill=datamodelmatch), binwidth=0.05) +
  facet_grid(ModelFamily ~ MixPeakLabel + DataFamily, scale="free_y") + theme_bw()
# dev.off()

### Caterpillar plots of pdet for EVERY rep (quite a lot).  All models are from the correct family.
# Facet rows    -- mixture status of data and model
# Facet columns -- Family and peakedness of data 
# Background shading trick
tempS3 <- subset(CD.melt, modelmix & datamix)
tp2 <- unique(tempS3[,c("ModelFamily","familypeak","datamodelmatch")])
# pdf(paste0("Misc/Cater",fileID,"families.pdf"), width=9.25, height=6.1)
print(ggplot(data=subset(tempS3, variable %in% c("X2.5.","X97.5."))) +
        geom_line(aes(x=value, y=as.factor(Rep), group=Rep)) + 
        geom_rect(aes(fill = datamodelmatch), data=tp2,
                  xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.15) +
        scale_fill_manual(values = c("red", "gray90")) + 
        facet_grid(ModelFamily~familypeak) + ylab("Simulation Replicate") + xlab("Probability") +
        geom_line(aes(x=value, y=as.factor(Rep)), size=I(0.6),color="orange",
                  data=subset(tempS3,variable %in% c("X25.","X75."))) + 
        geom_point(aes(x=value, y=as.factor(Rep)), color="blue", data=subset(tempS3,variable %in% "X50.")) +
        geom_vline(xintercept=truepdet, linetype=1, color="white", size=I(1.0)) +
        scale_shape_manual(values=c(2)) + 
        theme(axis.text.x=element_text(angle=45), axis.ticks.y=element_blank(), axis.text.y=element_blank(), legend.position="none")
)
# dev.off()
# There used to be a caterpillar plot of just peaked mixture gamma data to provide an example for the Oral Prelim



##### Plots of posterior parameter distributions
### Caterpillar plots of parameter estimates for all Reps of all correctly specified models... takes awhile to plot
# Load true parameter values
load(paste0("Simpars_p", fileID, ".Rdata")) # Note: this used to call ParValues from the Sim_1 folder
Simpars <- Simpars[,-which(names(Simpars) %in% c("meanTTD", "pdet"))]

# Posterior estimates
PE.melt <- melt(parmests, id.vars=11:ncol(parmests))
tempS4  <- subset(PE.melt, parameter !="lp__" & datamodelmatch)   # Only looking at correctly specified models (family and mixture)
tempS4  <- data.handle(tempS4)
tempS4$parameter[tempS4$parameter %in% c("alpha", "shape_k", "sigma_det")] <- "Shape"
# Re-order simdat for plotting purposes
tempS4$simdat <- factor(tempS4$simdat, 
                        levels = c("GFF", "LFF", "WFF", "EFF", "GFT", "LFT", "WFT",
                                   "GTF", "LTF", "WTF", "ETF", "GTT", "LTT", "WTT"))

# True values
vline.data <- na.omit(melt(Simpars, id.vars=c("model", "n_ints", "modcode"), value.name="z", variable.name="parameter"))
names(vline.data)[1] <- "simdat"
# Clumsy way to rename shape parameters 'Shape'
vline.data$parameter[vline.data$parameter %in% c("alpha", "shape_k", "sigma_det")] <- "alpha"
levels(vline.data$parameter)[levels(vline.data$parameter)=="alpha"] <- "Shape"
vline.data$parameter <- factor(vline.data$parameter)

### Plots posteriors for all datasets from true models
# # pdf(paste("Zero_Posteriors.pdf",sep=""), width=6.45, height=8.5)
# print(qplot(value, as.factor(Rep), geom="line", group=Rep, data=subset(tempS4, variable %in% c("X2.5.","X97.5.")), ylab="Simulation") +
#         facet_grid(simdat~parameter, scales="free_x") +
#         geom_line(aes(x=value, y=as.factor(Rep)), size=I(0.85),color="orange",
#                   data=subset(tempS4,variable %in% c("X25.","X75.")))  +
#         geom_point(aes(x=value, y=as.factor(Rep)), data=subset(tempS4,variable %in% "X50.")) +
#         # ggtitle("Caterpillar Plots of Parameter Posteriors - Models Correctly Specified") +
#         geom_vline(aes(xintercept = z), size=1.2, color="skyblue", vline.data))
# # dev.off()

# There used to be code here to plot the posterior median TTDDs and the true TTDDs
# It borrowed TTDD density function from the compile_addendum code

##### Plot Deviances deviously
names(dic)[names(dic)=="parameter"] <- "dicmethod"
dic2 <- melt(dic, id.vars <- 2:ncol(dic), value.name="DIC")
dic2 <- dic2[dic2$dicmethod != "Var" | dic2$variable=="dev1",]
dic2 <- data.handle(dic2)

# There used to be code to generate side-by-side plots of DIC calculated each of the three ways

# Compare DIC to the true model
true.dic <- subset(dic, subset=datamodelmatch, select=c("dev1", "Rep", "simdat", "dicmethod"))
dicdiff  <- join(subset(dic2, variable=="dev1"), true.dic)  # Joining by dicmethod, simdat, and Rep
dicdiff$dicdiff <- with(dicdiff, DIC - dev1)
dicdiff$simdat  <- factor(dicdiff$simdat, levels = c("GFF", "LFF", "WFF", "EFF", "GFT", "LFT", "WFT",
                                                     "GTF", "LTF", "WTF", "ETF", "GTT", "LTT", "WTT"))
# Plot Delta_DIC for mixture vs. non-mixture models:
# Only reliably identifies cases when data are mixture peaked
# pdf("DIC_Diff_mix.pdf", height=4)
print(
  ggplot(subset(dicdiff, familymatch & dicmethod=="ThetaMed" & !datamodelmatch)) +
    geom_jitter(aes(x=simdat, y=dicdiff), width=0.1, height=0) +
    ylab("DIC Difference") + xlab("Dataset") + 
    theme(axis.text.x=element_text(angle=45)) +
    geom_hline(yintercept=c(-10, -5, 5, 10), linetype="dashed") +
    geom_hline(yintercept=0)
  )
# dev.off()

# Plot Delta_DIC across families:
dicdiff$familypeak <- "No Peak"
dicdiff$familypeak[as.logical(last(dicdiff$simdat))] <- "Peak"
# pdf(paste("DIC_Diff_fam.pdf",sep=""))
print(
  ggplot(subset(dicdiff, datamix & dicmethod=="ThetaMed" & modelmix)) +
    geom_jitter(aes(x=ModelFamily, y=dicdiff), width=0.1, height=0) +
    ylab("DIC Difference") + xlab("Inference Family") + 
    theme(axis.text.x=element_text(angle=45)) +
    geom_hline(yintercept=c(-10, -5, 5, 10), linetype="dashed") +
    geom_hline(yintercept=0) +
    facet_grid(DataFamily~familypeak, scale="free")
)
# dev.off()



# ########## Generate Summary Tables -- Mix versus Non-Mix
# # Table of coverages for each dataset and inference model combo (averaged over reps) plus average posterior p-value
# cover.p <- ddply(p_glob, .(simdat, ModelFamily, DataFamily, modelmix), summarize, 
#                  C50=sum(abs(pval_pglob-0.5)<0.25)/length(unique(p_glob$Rep)), 
#                  C90=sum(abs(pval_pglob-0.5)<0.45)/length(unique(p_glob$Rep)), avg.postp=mean(pval_pglob))
# 
# Tbl.mix.m  <- subset(cover.p, DataFamily==ModelFamily & modelmix)[,c("simdat", "C50", "C90", "avg.postp")]  # Table for mixture inference models
# Tbl.mix.nm <- subset(cover.p, DataFamily==ModelFamily & !modelmix)[,c("simdat", "C50", "C90", "avg.postp")] # Table for non-mixture inference models
# # Mix.mean and NM.mean are average p.values, not average detection probabilities
# # Mix50, Mix90, etc. are coverage proportions
# names(Tbl.mix.m)[2:4] <- c("Mix50", "Mix90", "Mix.avg.postp")
# names(Tbl.mix.nm)[2:4] <- c("NM50", "NM90", "NM.avg.postp")
# Tbl.mix <- join(Tbl.mix.nm, Tbl.mix.m, by=c("simdat"))
# 
# # Average median(gamma) across Reps.
# # This table consists only of mixture inference models (but all data models). 
# # Note: the true value of gamma is 5/6 for the p50 simulations
# parmests  <- data.handle(parmests)
# avggamma  <- ddply(subset(parmests, parameter=="gamma" & DataFamily==ModelFamily),
#                   .(simdat), summarize, avg.med.gamma = mean(X50.))
# # Tbl.mix <- join(Tbl.mix, avggamma) # No longer reporting this optional column
# 
# ### Average posterior median estimate of pdet ###
# count.det <- data.handle(count.det)
# avgpdet <- ddply(subset(count.det, parameter=="p_global" & DataFamily==ModelFamily),
#                   .(simdat, modelmix), summarize, avg.med.pdet = mean(X50.))
# avgpdet <- join(subset(avgpdet, !modelmix), subset(avgpdet, modelmix), by="simdat")[,-c(2,4)] # Discard 'modelmix' columns
# names(avgpdet) <- c("simdat", "NM_pdet", "Mix_pdet")
# # Note that 'NM_pdet' and 'Mix_pdet' are equal to: mean(median(pdet))
# Tbl.mix <- join(Tbl.mix, avgpdet)[,c("simdat", "NM50", "NM90", "NM.avg.postp", "NM_pdet", 
#                                      "Mix50", "Mix90", "Mix.avg.postp", "Mix_pdet")] # Sort columns appropriately
# 
# # Instead of standard error (code deleted), we can alternatively use ratios of credible interval widths for p_det.  
# # The results are about the same.
# # Only for models where family matches data (mixture may or may not)
# Tbl.ciw <- subset(count.det, parameter=="p_global" & DataFamily==ModelFamily)[, c("X2.5.", "X97.5.", "Rep", "simdat", "modelmix", "datamodelmatch")]
# Tbl.ciw$width <- Tbl.ciw$X97.5. - Tbl.ciw$X2.5.
# # Join mixture and non-mixture fits onto same rows
# Tbl.ciw.wide <- join(Tbl.ciw[!Tbl.ciw$modelmix,-c(1,2,5)], Tbl.ciw[Tbl.ciw$modelmix,-c(1,2,5)], by=c("Rep", "simdat"))[,-3]
# names(Tbl.ciw.wide)[3:5] <- c("ciw.nomix", "datamix", "ciw.mix")
# Tbl.ciw <- ddply(Tbl.ciw.wide, .(simdat, datamix), summarize, m.div.nm = mean(ciw.mix/ciw.nomix), 
#                  nm.div.m = mean(ciw.nomix/ciw.mix))
# # Arrange so that ciw.ratio is the ratio of incorrect:correct
# Tbl.ciw$ciw.ratio <- Tbl.ciw$m.div.nm
# Tbl.ciw$ciw.ratio[Tbl.ciw$datamix] <- Tbl.ciw$nm.div.m[Tbl.ciw$datamix]
# # Joining the incorrect:correct ratio.  Can also join the mix:non-mix ratio or vice versa.
# Tbl.mix <- join(Tbl.mix, Tbl.ciw[,c("simdat","ciw.ratio")])
# 
# Tbl.dic <- subset(dicdiff, familymatch & dicmethod=="ThetaMed" & !datamodelmatch)
# Tbl.dic <- ddply(Tbl.dic, .(simdat), summarize, Dic_diff=mean(dicdiff))
# Tbl.mix <- join(Tbl.mix, Tbl.dic)
# 
# # Reorder...
# names(Tbl.mix) <- c("Dataset", "50%.nm", "90%.nm", "p_det quantile.nm", "Mean pdet.nm", 
#                     "50%.m", "90%.m", "p_det quantile.m", "Mean pdet.m", "CI Ratio", "DIC Diff")
# Tbl.mix <- Tbl.mix[match(c("GFF", "LFF", "WFF", "EFF", "GFT", "LFT", "WFT", "GTF", "LTF", "WTF", "ETF", "GTT", "LTT", "WTT"), Tbl.mix$Dataset), 
#                    c(1,5,4,2,3,9,8,6,7,10,11)]
# Tbl.mix$Dataset <- rep(c("Gamma", "Lognormal", "Weibull", "Exponential", "Gamma", "Lognormal", "Weibull"),2)
# 
# library(xtable)
# # The blanks pad the table for easier cut-and-paste into LaTex
# print(xtable(cbind(a="", b="", c="", Tbl.mix[,1:(ncol(Tbl.mix)-2)])), include.rownames=FALSE)
# 
# 
# 
# ########## Generate Summary Tables -- Compare Families
# cover.p2 <- ddply(subset(p_glob, datamix & modelmix), .(simdat, ModelFamily), summarize, 
#                   C50=sum(abs(pval_pglob-0.5)<0.25)/length(unique(p_glob$Rep)), #C75=sum(abs(pval_pglob-0.5)<0.375),
#                   C90=sum(abs(pval_pglob-0.5)<0.45)/length(unique(p_glob$Rep)), average=mean(pval_pglob))
# cover.p2 <- melt(cover.p2, id.vars=.(simdat, ModelFamily))
# 
# # Average median(gamma) across Reps.  True value is 5/6.
# # No surprise, this is most volatile for exponential mixture models.
# tempparm <- subset(parmests, parameter=="gamma" & modelmix & datamix)
# tempparm <- data.handle(tempparm)
# avggamma2 <- ddply(tempparm, .(simdat, ModelFamily, parameter), summarize, value = mean(X50.))
# # cover.p2 <- rbind(cover.p2, avggamma2)
# 
# # Average estimate of pdet for mixture data/models
# avgpdet.mix <- ddply(subset(count.det, parameter=="p_global" & modelmix & datamix), .(simdat, ModelFamily), 
#                      summarize, variable="median.pdet", value = mean(X50.))
# cover.p2 <- rbind(cover.p2, avgpdet.mix)
# 
# # Create table of credible interval widths
# Tbl.ciw2 <- subset(count.det, parameter=="p_global" & datamix & modelmix)[,
#                               c("fname","X2.5.", "X97.5.", "Rep", "simdat", "modelmix", "datamodelmatch")]
# Tbl.ciw2$width <- Tbl.ciw2$X97.5. - Tbl.ciw2$X2.5.
# Tbl.ciw2 <- data.handle(Tbl.ciw2)
# Tbl.ciw.true <- Tbl.ciw2[Tbl.ciw2$datamodelmatch,]
# names(Tbl.ciw.true)[match("width", names(Tbl.ciw.true))] <- "truewidth"
# Tbl.ciw2 <- join(Tbl.ciw2, Tbl.ciw.true[,c("simdat","Rep","truewidth")], by=c("simdat", "Rep"))
# Tbl.ciw2 <- ddply(Tbl.ciw2, .(fname, simdat, datamodelmatch, ModelFamily, DataFamily), summarize, ciw.ratio = mean(width/truewidth))
# 
# # Calculate Delta DIC
# Tbl.dic2 <- subset(dicdiff, modelmix & datamix & dicmethod=="ThetaMed")
# Tbl.dic2 <- ddply(Tbl.dic2, .(simdat, ModelFamily), summarize, variable="dicdiff", value=mean(dicdiff))
# cover.p2 <- rbind(cover.p2, Tbl.dic2)
# 
# # Combine the above into a single table
# Tbl.fam <- dcast(cover.p2, simdat ~ ModelFamily + variable, fun.aggregate=mean, value.var="value")
# 
# # Reorder rows and columns...
# Tbl.fam <- Tbl.fam[match(c("GTF", "LTF", "WTF", "ETF", "GTT", "LTT", "WTT"), Tbl.fam$simdat),]
# # Doing this now avoids some re-naming that R does with duplicate col. names and slicing...
# # Tbl.1a: Exponential and Gamma inference models
# # Tbl.1b: Lognormal and Weibull inference models
# Tbl.1a  <- Tbl.fam[, which(substr(names(Tbl.fam),1,1) %in% c("s", "E", "G"))]   
# Tbl.1b  <- Tbl.fam[, which(substr(names(Tbl.fam),1,1) %in% c("s", "L", "W"))]   
# names(Tbl.1a) <- names(Tbl.1b) <- c("Dataset", rep(c("50%", "90%", "Q(p_det)", "Median p_det", "DIC Diff"),2))
# # Reorder columns
# Tbl.1a <- Tbl.1a[,c(1,5,4,2,3,6,10,9,7,8,11)]; Tbl.1b <- Tbl.1b[,c(1,5,4,2,3,6,10,9,7,8,11)]
# Tbl.1a$Dataset <- Tbl.1b$Dataset <- c("Gamma","Lognormal","Weibull","Exponential","Gamma","Lognormal","Weibull")
# 
# library(xtable)
# print(xtable(cbind(a="", b="", Tbl.1a[,!grepl("DIC", names(Tbl.1a))])), include.rownames=FALSE) #, floating.environment="sidewaystable")
# print(xtable(cbind(a="", b="", Tbl.1b[,!grepl("DIC", names(Tbl.1a))])), include.rownames=FALSE)





########## Code to combine all simulation studies
# Make tables summarizing output across all 4 levels of pdet
library(ggplot2)
library(reshape2)
library(plyr)

# Combine tables of estimates of pdet & posterior p-values from across pdets
oest <- op <- NULL
fileID   <- c("50", "65", "80", "95")
for (j in 1:4){
  load(paste0("Outputs", fileID[j], ".Rdata"))
  outputs$count.det$pdet  <- fileID[j]
  outputs$post.unc.p$pdet <- fileID[j]
  oest <- rbind(oest, subset(outputs$count.det, parameter=="p_global")) # Overall (across true pdet values) estimates of pdet
  op   <- rbind(op, outputs$post.unc.p)                                 # Overall summary of posterior p-values
}

# Compile values across reps and combine tables
# The tables are 'rbind'-ed instead of 'join'-ed, because that facilitates the 'dcast' step
oci  <- ddply(oest, .(fname, simdat, modelmix, datamix, datamodelmatch, pdet), summarize,  # Used later to explore CIs
              q2.5 = round(mean(X2.5.),2), q25 = round(mean(X25.),2), q50 = round(mean(X50.),2), q75 = round(mean(X75.),2), q97.5 = round(mean(X97.5.),2)) 
oest <- ddply(oest, .(fname, simdat, modelmix, datamix, datamodelmatch, pdet), summarize, value = round(mean(X50.),2))  # Average of medians
op   <- ddply(subset(op, parameter=="p_global"), .(fname, simdat, modelmix, datamix, datamodelmatch, pdet), summarize, 
              value = sum(abs(pval-0.5) < 0.45)/length(pval))  # Coverage of 90% credible intervals
oest$statistic <- "Med.p"
op$statistic   <- "C90"
ov   <- rbind(oest, op)  
ov.c <- dcast(ov, fname+simdat+modelmix+datamix+datamodelmatch ~ statistic + pdet, value.var="value")

# ##### Make Output Tables
# # Data handling
# ov.c <- data.handle(ov.c)
# ov.c$datapeaked <- as.logical(last(ov.c$simdat))
# ov.c$datapeaked[ov.c$DataFamily=="Exponential"] <- TRUE
# ov.c <- ov.c[with(ov.c, order(datamix, datapeaked, DataFamily)),]  # Sort data in correct order for output tables
# 
# # Subset the ov.c to obtain correctly arranged chunks of the final publishable table
# ov.subs <- list(
#   mix.left.cov   = subset(ov.c, subset = familymatch & !modelmix, select = c(DataFamily, grep("^C90", names(ov.c)))),
#   mix.left.medp  = subset(ov.c, subset = familymatch & !modelmix, select = c(DataFamily, grep("^Med", names(ov.c)))),
#   mix.right.cov  = subset(ov.c, subset = familymatch & modelmix,  select = c(DataFamily, grep("^C90", names(ov.c)))),
#   mix.right.medp = subset(ov.c, subset = familymatch & modelmix,  select = c(DataFamily, grep("^Med", names(ov.c)))),
#   fam.ul.cov     = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Exponential", select = c(DataFamily, grep("^C90", names(ov.c)))),
#   fam.ul.medp    = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Exponential", select = c(DataFamily, grep("^Med", names(ov.c)))),
#   fam.ll.cov     = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Gamma", select = c(DataFamily, grep("^C90", names(ov.c)))),
#   fam.ll.medp    = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Gamma", select = c(DataFamily, grep("^Med", names(ov.c)))),
#   fam.ur.cov     = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Lognormal", select = c(DataFamily, grep("^C90", names(ov.c)))),
#   fam.ur.medp    = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Lognormal", select = c(DataFamily, grep("^Med", names(ov.c)))),
#   fam.lr.cov     = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Weibull", select = c(DataFamily, grep("^C90", names(ov.c)))),
#   fam.lr.medp    = subset(ov.c, subset = datamix & modelmix & ModelFamily=="Weibull", select = c(DataFamily, grep("^Med", names(ov.c))))
# )
# 
# # Bind the chunks together to create the final tables
# tbl1  <- with(ov.subs, cbind(mix.left.medp, mix.left.cov[,-1], mix.right.medp[,-1], mix.right.cov[,-1]))
# tbl2a <- with(ov.subs, cbind(fam.ul.medp, fam.ul.cov[,-1], fam.ur.medp[,-1], fam.ur.cov[,-1]))
# tbl2b <- with(ov.subs, cbind(fam.ll.medp, fam.ll.cov[,-1], fam.lr.medp[,-1], fam.lr.cov[,-1]))
# 
# library(xtable)
# print(xtable(cbind(a="", b="", c="", tbl1)), include.rownames=FALSE)
# 
# print(xtable(cbind(a="", b="", tbl2a)), include.rownames=FALSE)
# print(xtable(cbind(a="", b="", tbl2b)), include.rownames=FALSE)

##### Make Output Plots
ov          <- data.handle(ov)
ov$peaked   <- substr(ov$simdat,3,3)
ov$peaked2  <- factor(ov$peaked, labels=c("Nonpeaked", "Peaked"))
ov$peaked[ov$DataFamily=="Exponential"] <- "E"
ov$peaked2[ov$DataFamily=="Exponential"] <- "Nonpeaked"
ov$datamix  <- factor(ov$datamix, labels=c("Non-Mixed", "Mixed"))
ov$peaked   <- factor(ov$peaked, labels=c("Exponential", "Nonpeaked", "Peaked"))
ov$peaked   <- factor(ov$peaked, levels=levels(ov$peaked)[c(2,1,3)])
ov$mt       <- with(ov, paste0(modelmix, ModelFamily))
ov$mt       <- factor(ov$mt, labels=c("Non-mix Exponential", "Non-mix Gamma", "Non-mix Lognormal", "Non-mix Weibull",
                                      "Mix Exponential", "Mix Gamma", "Mix Lognormal", "Mix Weibull"))

# Mix - NoMix plots
library(gridExtra)
library(grid)
# Following code from https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height) # Legend height
  lwidth <- sum(legend$width)   # Legend width
  gl <- lapply(plots, function(x) x + theme(legend.position="none")) # List of 'plots' with legends removed
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- switch(position, # It's like an IF statment to run arrangeGrob 
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)
                                           ),
                     "right"  = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 2,
                                            widths = unit.c(unit(1, "npc") - lwidth, lwidth)
                                           )
  )
  grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}

lp <- ggplot(subset(ov, statistic=="Med.p" & familymatch)) + 
  geom_point(aes(x=as.numeric(pdet), y=value, shape=mt), size=2.25) +
  geom_segment(aes(x=50 , y=0.5 , xend=95, yend=0.95), color="gray60", linetype=2, size=0.5) +
  facet_grid(datamix~peaked) +
  scale_shape_manual(values=c(15,16,17,18,0,1,2,5), name="Inference Model") +
  scale_x_continuous(breaks = seq(50, 95, 15)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw() + 
  ylab("Estimated Detection Probability") + xlab("Actual Detection Probability")

rp <- ggplot(subset(ov, statistic=="C90" & familymatch)) + 
  geom_hline(aes(yintercept=0.835), linetype=2, size=1, color="gray60") +
  geom_hline(aes(yintercept=0.955), linetype=2, size=1, color="gray60") +
  geom_point(aes(x=as.numeric(pdet), y=value, shape=mt), size=2.25) +
  facet_grid(datamix~peaked) +
  scale_shape_manual(values=c(15,16,17,18,0,1,2,5), name="Inference Model") +
  scale_x_continuous(breaks = seq(50, 95, 15)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits=c(-0.05, 1.05)) +
  theme_bw() + ylab("Coverage") + xlab("Actual Detection Probability")

pdf("mixture_fig.pdf", width=10, height=4.7)
grid_arrange_shared_legend(lp, rp, position="bottom")
dev.off()

# Family comparison plots
lf <- ggplot(subset(ov, statistic=="Med.p" & datamix=="Mixed" & modelmix)) + 
  geom_segment(aes(x=50 , y=0.5 , xend=95, yend=0.95), color="gray60", linetype=2, size=0.5) +
  geom_point(aes(x=as.numeric(pdet), y=value, shape=ModelFamily), size=2.25) +
  facet_grid(DataFamily~peaked2) +
  scale_shape_manual(values=c(0,1,2,5), name="Inference Model") +
  scale_x_continuous(breaks = seq(50, 95, 15)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw() + ylab("Estimated Detection Probability") + xlab("Actual Detection Probability")

rf <- ggplot(subset(ov, statistic=="C90" & datamix=="Mixed" & modelmix)) + 
  geom_hline(aes(yintercept=0.835), linetype=2, size=1, color="gray60") +
  geom_hline(aes(yintercept=0.955), linetype=2, size=1, color="gray60") +
  geom_point(aes(x=as.numeric(pdet), y=value, shape=ModelFamily), size=2.25) +
  facet_grid(DataFamily~peaked2) +
  scale_shape_manual(values=c(0,1,2,5), name="Inference Model") +
  scale_x_continuous(breaks = seq(50, 95, 15)) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits=c(-0.05, 1.05)) +
  theme_bw() + ylab("Coverage") + xlab("Actual Detection Probability")

pdf("family_fig.pdf", width=10, height=4.7)
grid_arrange_shared_legend(lf, rf, position="bottom")
dev.off()
 
 
########## Looking at CI Widths

# RESULTS:
# (1) CI widths don't much change b/t pdet50 and pdet 65.  NP/NM, NP/Mix, and Pk/Mix
#     In these scenarios, Mix CI's at low pdets are wider by about 5, 5, and 10 percent, respectively
# (2) For expo and Peaked distributions, the widths do decrease.
#     There's little difference between Mix and NM widths, except in the Expo/Mix data scenario

oci   <- data.handle(oci)
oci$peaked   <- substr(oci$simdat,3,3)
oci$peaked2  <- factor(oci$peaked, labels=c("Nonpeaked", "Peaked"))
oci$peaked[oci$DataFamily=="Exponential"] <- "E"
oci$peaked2[oci$DataFamily=="Exponential"] <- "Peaked"
oci$datamix  <- factor(oci$datamix, labels=c("Non-Mixed", "Mixed"))
oci$peaked   <- factor(oci$peaked, labels=c("Exponential", "Nonpeaked", "Peaked"))
oci$peaked   <- factor(oci$peaked, levels=levels(oci$peaked)[c(2,1,3)])
oci$mt       <- with(oci, paste0(modelmix, ModelFamily))
oci$mt       <- factor(oci$mt, labels=c("NM Expo", "NM Gamma", "NM Lognormal", "NM Weibull",
                                        "Mix Expo", "Mix Gamma", "Mix Lognormal", "Mix Weibull"))
oci$width95  <- oci$q97.5 - oci$q2.5
oci$width50  <- oci$q75 - oci$q25

ggplot(subset(oci, familymatch)) + 
  # geom_point(aes(x=as.numeric(pdet), y=width50, shape=mt), size=2) +
  geom_point(aes(x=as.numeric(pdet), y=q2.5, shape=mt), size=2) +
  geom_point(aes(x=as.numeric(pdet), y=q50, shape=mt), size=2, color="orange") +
  geom_point(aes(x=as.numeric(pdet), y=q97.5, shape=mt), size=2) +
  geom_segment(aes(x=50 , y=0.5 , xend=95, yend=0.95), linetype=2) +
  facet_grid(datamix~peaked) +
  scale_shape_manual(values=c(15,16,17,18,0,1,2,5), name="Inference Model") +
  theme_bw() + ylab("Estimated Detection Probability") + xlab("Actual Detection Probability")




########## Looking at posterior p-values associated with p= seq(0.4, 0.8, 0.1)

oqs <- ops <- oest <- NULL
fileID     <- c("50", "65", "80", "95")
for (j in 1:4){
  load(paste0("Outputs", fileID[j], ".Rdata"))
  outputs$post.unc.p$pdet <- fileID[j]
  outputs$count.det$pdet  <- fileID[j]
  oqs      <- rbind(oqs, subset(outputs$post.unc.p, parameter %in% c("Pval40", "Pval50", "Pval60", "Pval70", "Pval80")))
  ops      <- rbind(ops, subset(outputs$post.unc.p, parameter %in% c("p_global")))
  oest     <- rbind(oest, subset(outputs$count.det, parameter=="p_global"))
}
rownames(oqs) <- rownames(ops) <- rownames(oest) <- NULL

library(plyr)
library(ggplot2)
# Here's a plot of posterior p-values for p=(40, 50, 60, 70, 80), but I don't know what to make of it
# or even how to define what is a good pattern
# Problem is that the y-axis is a posterior p-value that should correlate with the x-axis
oqs.sum <- ddply(oqs, .(fname, pdet, simdat, modelmix, datamix, datamodelmatch, parameter), summarize, M = mean(pval))
qplot(pdet, M, geom="line", data=oqs.sum, group=parameter, color=parameter, facets=~fname)

ops$TTDD <- substr(ops$fname,1,1)
ops$TTDD[!ops$modelmix] <- "Non"
ops      <- join(ops, subset(oest, select=c("fname", "Rep", "pdet", "X50.")), by=c("Rep", "fname", "pdet"))
ops2     <- ddply(ops, .(simdat, TTDD, pdet), summarize, med=mean(X50.), pv=mean(pval))
# Boxplots of posterior pvalues (across Reps) faceted by dataset and inference model... quite an eyeful
qplot(pdet, pval, geom="boxplot", data=ops, facets=simdat~TTDD)
# Point plot of pval vs. estimate of pdet
# As usual, it shows that exponential models are not robust,
# Peaked data are more accurately estimated,
# and Nonpeaked data are a mess.
# It's pretty, but kind of difficult to interpret.
# Can change the 'data' option to look at peaked, nonpeaked, or both
xpts <- c(0.5, 0.65, 0.8, 0.95)
ypts <- rep(0.5, 4)
qplot(X50., pval, geom="point", facets=simdat~TTDD, color=pdet, size=I(0.5), data=subset(ops, substr(simdat,3,3)=="F")) +
  geom_point(aes(x=med, y=pv), data=subset(ops2, substr(simdat,3,3)=="F"), size=I(3)) +
  geom_point(aes(x=xpts[1], y=ypts[1]), size=I(3), shape=I(4), color=1) + geom_point(aes(x=xpts[2], y=ypts[2]), size=I(3), shape=I(4), color=2) +
  geom_point(aes(x=xpts[3], y=ypts[3]), size=I(3), shape=I(4), color=3) + geom_point(aes(x=xpts[4], y=ypts[4]), size=I(3), shape=I(4), color=4) +
  geom_hline(aes(yintercept=0.5)) +
  geom_hline(aes(yintercept=0.025), linetype=2) +
  geom_hline(aes(yintercept=0.975), linetype=2)
