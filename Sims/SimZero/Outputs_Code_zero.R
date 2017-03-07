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

truepdet <- 0.5

simtable <- read.csv("S0_simtable.csv")  # Load simtable
load("Outputs_p50.Rdata")  # Load model fits
for(i in 1:length(outputs)) assign(names(outputs)[i], outputs[[i]]) # Convert model fits from list into objects

##### There used to be a section here for posterior predictive p-values
# See First_simulations folder for code

###### Posterior p-values for global estimates of pdet
# There's a choice here. 
# Choice one: use the p-values associated with "p_global" and make comparisons to the data realized in simulation.
# Choice two: use the p-values associated with (misnamed) "p80".  This seems more correct.
p_glob <- subset(post.unc.p, parameter=="p80")

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
pglobs  <- join(true.pglobal, simtable)  # Construct a table of true pdet values with same variables as posterior estimates
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
# pdf("pdet_cater_correct.pdf", width=7.75, height=5.1)
print(ggplot(data=subset(tempS2, variable %in% c("X2.5.","X97.5."))) + 
        geom_line(aes(x=value, y=as.factor(Rep), group=Rep)) + 
        geom_rect(aes(fill = datamodelmatch), data=tp,
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.15) +
        scale_fill_manual(values = c("red", "gray90")) + 
        facet_grid(dmix.mmix~familypeak) + ylab("Simulation") + xlab("Probability") +
        geom_line(aes(x=value, y=as.factor(Rep), group=Rep), size=I(0.85), color="orange",
                  data=subset(tempS2,variable %in% c("X25.","X75."))) + 
        geom_point(aes(x=value, y=as.factor(Rep)), data=subset(tempS2,variable %in% "X50.")) +
        theme(axis.text.x=element_text(angle=45)) + geom_vline(xintercept=truepdet, linetype=1, color="blue") +
        scale_shape_manual(values=c(2)) + theme(legend.position="none")
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
# pdf("pdet_cater_family.pdf", width=7.75, height=5.1)
print(ggplot(data=subset(tempS3, variable %in% c("X2.5.","X97.5."))) +
        geom_line(aes(x=value, y=as.factor(Rep), group=Rep)) + 
        geom_rect(aes(fill = datamodelmatch), data=tp2,
                  xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf,alpha = 0.15) +
        scale_fill_manual(values = c("red", "gray90")) + 
        facet_grid(ModelFamily~familypeak) + ylab("Simulation") + xlab("Probability") +
        geom_line(aes(x=value, y=as.factor(Rep)), size=I(0.85),color="orange",
                  data=subset(tempS3,variable %in% c("X25.","X75."))) + 
        geom_point(aes(x=value, y=as.factor(Rep)), data=subset(tempS3,variable %in% "X50.")) +
        # geom_point(aes(x=value, y=as.factor(Rep), shape=variable), color="blue", data=subset(tempS3,variable %in% "True")) +
        # ggtitle("Posterior Estimates of Detection Probability") +
        theme(axis.text.x=element_text(angle=45)) + geom_vline(xintercept=truepdet, linetype=1, color="blue") +
        scale_shape_manual(values=c(2)) + theme(legend.position="none"))
# dev.off()

# There used to be a caterpillar plot of just peaked mixture gamma data to provide an example for the Oral Prelim



##### Plots of posterior parameter distributions
### Caterpillar plots of parameter estimates for all Reps of all correctly specified models... takes awhile to plot
# Load true parameter values
load("Simpars_p50.Rdata") # Note: this used to call ParValues from the Sim_1 folder
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



########## Generate Summary Tables -- Mix versus Non-Mix
# Table of coverages for each dataset and inference model combo (averaged over reps) plus average posterior p-value
cover.p <- ddply(p_glob, .(simdat, ModelFamily, DataFamily, modelmix), summarize, 
                 C50=sum(abs(pval_pglob-0.5)<0.25)/length(unique(p_glob$Rep)), 
                 C90=sum(abs(pval_pglob-0.5)<0.45)/length(unique(p_glob$Rep)), average=mean(pval_pglob))

Tbl.mix.m  <- subset(cover.p, DataFamily==ModelFamily & modelmix)[,c("simdat", "C50", "C90", "average")]  # Table for mixture inference models
Tbl.mix.nm <- subset(cover.p, DataFamily==ModelFamily & !modelmix)[,c("simdat", "C50", "C90", "average")] # Table for non-mixture inference models
# Mix.mean and NM.mean are average p.values, not average detection probabilities
# Mix50, Mix90, etc. are coverage proportions
names(Tbl.mix.m)[2:4] <- c("Mix50", "Mix90", "Mix.mean")
names(Tbl.mix.nm)[2:4] <- c("NM50", "NM90", "NM.mean")
Tbl.mix <- join(Tbl.mix.nm, Tbl.mix.m, by=c("simdat"))

# Average median(gamma) across Reps.
# This table consists only of mixture inference models (but all data models). 
# Note: the true value of gamma is 5/6 for the p50 simulations
parmests  <- data.handle(parmests)
avggamma  <- ddply(subset(parmests, parameter=="gamma" & DataFamily==ModelFamily),
                  .(simdat), summarize, avg.med.gamma = mean(X50.))
# Tbl.mix <- join(Tbl.mix, avggamma) # No longer reporting this optional column

### Average posterior median estimate of pdet ###
count.det <- data.handle(count.det)
avgpdet <- ddply(subset(count.det, parameter=="p_global" & DataFamily==ModelFamily),
                  .(simdat, modelmix), summarize, avg.med.pdet = mean(X50.))
avgpdet <- join(subset(avgpdet, !modelmix), subset(avgpdet, modelmix), by="simdat")[,-c(2,4)] # Discard 'modelmix' columns
names(avgpdet) <- c("simdat", "NM_pdet", "Mix_pdet")
# Note that 'NM_pdet' and 'Mix_pdet' are equal to: mean(median(pdet))
Tbl.mix <- join(Tbl.mix, avgpdet)[,c(1:4,8,5:7,9)] # Sort columns appropriately

# Instead of standard error (code deleted), we can alternatively use ratios of credible interval widths for p_det.  
# The results are about the same.
# Only for models where family matches data (mixture may or may not)
Tbl.ciw <- subset(count.det, parameter=="p_global" & DataFamily==ModelFamily)[, c("X2.5.", "X97.5.", "Rep", "simdat", "modelmix", "datamodelmatch")]
Tbl.ciw$width <- Tbl.ciw$X97.5. - Tbl.ciw$X2.5.
# Join mixture and non-mixture fits onto same rows
Tbl.ciw.wide <- join(Tbl.ciw[!Tbl.ciw$modelmix,-c(1,2,5)], Tbl.ciw[Tbl.ciw$modelmix,-c(1,2,5)], by=c("Rep", "simdat"))[,-3]
names(Tbl.ciw.wide)[3:5] <- c("ciw.nomix", "datamix", "ciw.mix")
Tbl.ciw <- ddply(Tbl.ciw.wide, .(simdat, datamix), summarize, m.div.nm = mean(ciw.mix/ciw.nomix), 
                 nm.div.m = mean(ciw.nomix/ciw.mix))
# Arrange so that ciw.ratio is the ratio of incorrect:correct
Tbl.ciw$ciw.ratio <- Tbl.ciw$m.div.nm
Tbl.ciw$ciw.ratio[Tbl.ciw$datamix] <- Tbl.ciw$nm.div.m[Tbl.ciw$datamix]
# Joining the incorrect:correct ratio.  Can also join the mix:non-mix ratio or vice versa.
Tbl.mix <- join(Tbl.mix, Tbl.ciw[,c("simdat","ciw.ratio")])

Tbl.dic <- subset(dicdiff, familymatch & dicmethod=="ThetaMed" & !datamodelmatch)
Tbl.dic <- ddply(Tbl.dic, .(simdat), summarize, Dic_diff=mean(dicdiff))
Tbl.mix <- join(Tbl.mix, Tbl.dic)

# Reorder...
names(Tbl.mix) <- c("Dataset", "50%", "90%", "p_det quantile", "Mean pdet", "50%.m", "90%.m", "p_det quantile.m", "Mean pdet.m", "CI Ratio", "DIC Diff")
Tbl.mix <- Tbl.mix[match(c("GFF", "LFF", "WFF", "EFF", "GFT", "LFT", "WFT", "GTF", "LTF", "WTF", "ETF", "GTT", "LTT", "WTT"), 
                         Tbl.mix$Dataset), c(1,5,4,2,3,9,8,6,7,10,11)]

library(xtable)
print(xtable(Tbl.mix), include.rownames=FALSE)













##### Generate Summary Tables -- Compare Families
cover.p2 <- ddply(p_glob[p_glob$datamix & p_glob$modelmix,], .(simdat, ModelFamily), summarize, 
                  C50=sum(abs(pval_pglob-0.5)<0.25)/length(unique(p_glob$Rep)), #C75=sum(abs(pval_pglob-0.5)<0.375),
                  C90=sum(abs(pval_pglob-0.5)<0.45)/length(unique(p_glob$Rep)), average=mean(pval_pglob))
cover.p2 <- melt(cover.p2, id.vars=.(simdat, ModelFamily))

# Mean(median(gamma)) across Reps
# tempparm <- parmests[parmests$parameter=="gamma" & parmests$modelmix & parmests$datamix,]
# tempparm$Letter <- substr(tempparm$fname,1,1); tempparm <- join(tempparm, familytab); tempparm$variable <- "gamma"
# avggamma2 <- ddply(tempparm, .(simdat, ModelFamily, variable), summarize, value = mean(X50.))
# cover.p2 <- rbind(cover.p2, avggamma2)
count.det$Letter <- substr(count.det$fname,1,1); count.det <- join(count.det, familytab)
avgpdet.mix <- ddply(count.det[count.det$parameter=="p_global" & count.det$modelmix & count.det$datamix,],
                 .(simdat, ModelFamily), summarize, variable="median.pdet", value = mean(X50.))
cover.p2 <- rbind(cover.p2, avgpdet.mix)

# Tbl.se2 <- CD.melt[CD.melt$datamix & CD.melt$modelmix & CD.melt$variable=="se",][,c("Rep", "simdat", "value", "ModelFamily", "DataFamily")]
# Tbl.se2.wide <- join(Tbl.se2[,-5], Tbl.se2[Tbl.se2$ModelFamily==Tbl.se2$DataFamily,-(4:5)], by=c("Rep", "simdat"))
# names(Tbl.se2.wide)[3:5] <- c("se.model", "ModelFamily", "se.correct")
# Tbl.se2 <- ddply(Tbl.se2.wide, .(simdat, ModelFamily), summarize, variable="se.ratio", value = mean(se.model/se.correct))
# cover.p2 <- rbind(cover.p2, Tbl.se2)

# For the moment, these are not in the published table.  I can do that, but it requires
# subsetting columns and rbind-ing with cover.p2
Tbl.ciw2 <- count.det[count.det$parameter=="p_global" & count.det$datamix & count.det$modelmix,
                     c("fname","X2.5.", "X97.5.", "Rep", "simdat", "modelmix", "datamodelmatch")]
Tbl.ciw2$width <- Tbl.ciw2$X97.5. - Tbl.ciw2$X2.5.
Tbl.ciw2$Letter <- substr(Tbl.ciw2$simdat,1,1); Tbl.ciw2$Letter2 <- substr(Tbl.ciw2$fname,1,1)
Tbl.ciw2 <- join(join(Tbl.ciw2, familytab, by="Letter"), familytab2, by="Letter2")
Tbl.ciw.true <- Tbl.ciw2[Tbl.ciw2$datamodelmatch,]; names(Tbl.ciw.true)[match("width", names(Tbl.ciw.true))] <- "truewidth"
Tbl.ciw2 <- join(Tbl.ciw2, Tbl.ciw.true[,c("simdat","Rep","truewidth")], by=c("simdat", "Rep"))
Tbl.ciw2 <- ddply(Tbl.ciw2, .(fname, simdat, datamodelmatch, ModelFamily, DataFamily), summarize, ciw.ratio = mean(width/truewidth))

Tbl.dic2 <- dicdiff[dicdiff$modelmix & dicdiff$datamix & dicdiff$dicmethod=="Median",]
Tbl.dic2 <- ddply(Tbl.dic2, .(simdat, ModelFamily), summarize, variable="dicdiff", value=mean(dicdiff))
cover.p2 <- rbind(cover.p2, Tbl.dic2)

Tbl.fam <- dcast(cover.p2, simdat ~ ModelFamily + variable, fun.aggregate=mean, value.var="value")

# Reorder...
Tbl.fam <- Tbl.fam[match(c("GTF", "LTF", "WTF", "ETF", "GTT", "LTT", "WTT"), Tbl.fam$simdat),]
Tbl.1a <- Tbl.fam[,c(1:11)]   # Doing this now avoids some re-naming that R does with duplicate col. names and slicing...
Tbl.1b <- Tbl.fam[,c(1,12:21)]
names(Tbl.1a) <- names(Tbl.1b) <- c("Dataset", rep(c("50%", "90%", "Q(p_det)", "Median p_det", "DIC Diff"),2))
# Reorder columns
Tbl.1a <- Tbl.1a[,c(1,5,4,2,3,6,10,9,7,8,11)]; Tbl.1b <- Tbl.1b[,c(1,5,4,2,3,6,10,9,7,8,11)]
Tbl.1a$Dataset <- Tbl.1b$Dataset <- c("Gamma","Lognormal","Weibull","Exponential","Gamma","Lognormal","Weibull")
# This used to have se-widths & DIC.  I haven't updated it yet to do CI-widths, so DIC is smushed in above.
# Tbl.2 <- Tbl.fam[,c(1,seq(6,24,6),seq(7,25,6))]
#names(Tbl.2) <- c("Dataset",rep(as.character(familytab$ModelFamily),2))

library(xtable)
print(xtable(Tbl.1a), include.rownames=FALSE) #, floating.environment="sidewaystable")
print(xtable(Tbl.1b), include.rownames=FALSE) #, floating.environment="sidewaystable")
#print(xtable(Tbl.2), include.rownames=FALSE) #, floating.environment="sidewaystable")




# Intervals for Abundance
# I'm comparing to the ACTUAL abundance, not the expected, because these are based on the 'uncounted' credible intervals
# If I wanted expected counts, I ought to compare to estimates of Intercept^A
# Frankly, the differences are negligible
Abund <- NULL
cu <- function(u,p) p*u/(1-p)
exp.N <- 381 * exp(ParValues[[1]]$intcpt_a)   # 381 sites.  Same intcpt_a used in all simulations
# Code to compile all datasets along with actual counts, actual abundance, and expected abundance
for(Rep in 1:length(unique(count.det$Rep))){
  for(dc in unique(count.det$datacode)){
    df <- count.det[count.det$Rep==Rep & count.det$datacode==dc,]
    actual.count <- with(df, cu(X50.[1], X50.[2]))
    actual.N <- actual.count / (true.pglobal$value[true.pglobal$Rep==Rep & true.pglobal$datacode==dc])
    df[,c("mean", "X2.5.", "X25.", "X50.", "X75.", "X97.5.")] <- 
      (df[,c("mean", "X2.5.", "X25.", "X50.", "X75.", "X97.5.")] + actual.count) / actual.N
    Abund <- rbind(Abund, cbind(
      df[df$parameter=="uncounted",], actual.count=actual.count, actual.N=actual.N, exp.N=exp.N/actual.N))
  }
}
Abund$mixpeak <- substr(Abund$simdat,2,3)
Abund$datfamily <- substr(Abund$simdat,1,1)
# Melt data
Ab.long <- melt(Abund, id.vars=c("Rep", "simdat", "mixpeak", "ModelFamily", "datfamily", "modelmix", "datamix", "datamodelmatch"),
                measure.vars=c("mean", "X2.5.", "X25.", "X50.", "X75.", "X97.5.", "exp.N"))
# Summarize by data/model across Reps
Ab.long.ply <- ddply(Ab.long, .(simdat, mixpeak, ModelFamily, datfamily, modelmix, datamix, datamodelmatch, variable), summarize,
                     Ratio=round(mean(value),2), logRatio=round(mean(log(value)),2))
# Add a factor level for 'NonMix' inference models
Ab.long.ply$ModelFamily <- factor(Ab.long.ply$ModelFamily, levels=c(levels(Ab.long.ply$ModelFamily), "NonMix"))
Ab.long.ply$ModelFamily[!Ab.long.ply$modelmix] <- "NonMix"

# Plots of CIs, median, and expected (abundance scale and log scale)
ggplot() + 
  geom_line(aes(ModelFamily, Ratio), data=Ab.long.ply[Ab.long.ply$variable %in% c("X2.5.", "X97.5."),], color="black", size=1.25) +
  geom_line(aes(ModelFamily, Ratio), data=Ab.long.ply[Ab.long.ply$variable %in% c("X25.", "X75."),], color="orange", size=1.25) +
  geom_point(aes(ModelFamily, Ratio), data=Ab.long.ply[Ab.long.ply$variable %in% c("X50."),], color="dodgerblue3", size=2.5) +
  geom_point(aes(ModelFamily, Ratio), data=Ab.long.ply[Ab.long.ply$variable %in% c("mean"),], color="chocolate2", size=2.5) +
  geom_point(aes(ModelFamily, Ratio), data=Ab.long.ply[Ab.long.ply$variable %in% c("exp.N"),], color="maroon2", size=2.5) +
  facet_grid(datfamily~mixpeak) + coord_flip()

ggplot() + 
  geom_line(aes(ModelFamily, logRatio), data=Ab.long.ply[Ab.long.ply$variable %in% c("X2.5.", "X97.5."),], color="black", size=1.25) +
  geom_line(aes(ModelFamily, logRatio), data=Ab.long.ply[Ab.long.ply$variable %in% c("X25.", "X75."),], color="orange", size=1.25) +
  geom_point(aes(ModelFamily, logRatio), data=Ab.long.ply[Ab.long.ply$variable %in% c("X50."),], color="dodgerblue3", size=2.5) +
  geom_point(aes(ModelFamily, logRatio), data=Ab.long.ply[Ab.long.ply$variable %in% c("mean"),], color="chocolate2", size=2.5) +
  geom_point(aes(ModelFamily, logRatio), data=Ab.long.ply[Ab.long.ply$variable %in% c("exp.N"),], color="maroon2", size=2.5) +
  facet_grid(datfamily~mixpeak)

Ab.wide.med   <- dcast(Ab.long.ply[Ab.long.ply$variable=="X50.",],   simdat ~ ModelFamily, value.var="Ratio")
Ab.wide.lower <- dcast(Ab.long.ply[Ab.long.ply$variable=="X2.5.",],  simdat ~ ModelFamily, value.var="Ratio")
Ab.wide.upper <- dcast(Ab.long.ply[Ab.long.ply$variable=="X97.5.",], simdat ~ ModelFamily, value.var="Ratio")

# A function to plot caterpillar plots for Scaled Abundance... Jarad felt it was switching topics, but keep it in case
# X-axis is scaled to the ACTUAL abundance for each Rep.  So, the actual value is at X=1 (or X=0 of LT=log-transform=TRUE)
ab.plot <- function(Sdat, LT=FALSE){
  d.Sdat <- subset(Ab.long, simdat %in% Sdat & modelmix)
  levels(d.Sdat$variable)[c(1,4)] <- c("Mean", "Median")
  names(d.Sdat)[names(d.Sdat)=="variable"] <- "Estimator"
  Lab <- "Scaled Abundance"
  if(LT){d.Sdat$value <- log(d.Sdat$value); Lab="log(Scaled Abundance)"}
  print(
    ggplot() + ylab("Inference Model") + xlab(Lab) +
      geom_line(aes(x=value, y=Rep, group=as.factor(Rep)), 
                data=subset(d.Sdat, Estimator %in% c("X2.5.", "X97.5."))) +
      geom_line(aes(x=value, y=Rep, group=as.factor(Rep)), 
                data=subset(d.Sdat, Estimator %in% c("X25.", "X75.")), color="orange", size=0.85) + 
      geom_point(aes(x=value, y=Rep, color=Estimator, shape=Estimator), 
                 data=subset(d.Sdat, Estimator %in% c("Median")), size=I(2)) + 
      facet_grid(ModelFamily~simdat) + scale_color_manual(values=c("blue", "orchid3")) + theme_bw() +
      ggtitle("Posterior Estimates of log10(Abundace)") + theme(legend.position='none')
  )
}
pdf("../../Oral Prelim/abund_cater_family_prelim.pdf", width=7.5, height=7)
ab.plot(c("ETF","GTF", "GTT"), LT=T)
dev.off()

# Plots from all Reps and Models... mixpeakmix = datamix--peak--modelmix
Abund$ModelFamily <- factor(Abund$ModelFamily, levels=c(levels(Abund$ModelFamily), "NonMix"))
Abund$ModelFamily[!Abund$modelmix] <- "NonMix"
Abund$mixpeakmix <- paste0(Abund$mixpeak, substr(Abund$modelmix,1,1))
qplot(X50., Rep, geom="point", data=Abund[Abund$datfamily==Abund$Letter,], color=modelmix) + 
  geom_point(aes(x=X2.5., y=Rep), shape=I(4)) +
  geom_point(aes(x=X97.5., y=Rep), shape=I(4)) +
  facet_grid(datfamily~mixpeak, scale="free_x")

Abund$ciw <- Abund$X97.5. - Abund$X2.5.
smurf.mix <- Abund[Abund$modelmix & Abund$Letter==Abund$datfamily,]
smurf.NM <- Abund[!Abund$modelmix,]
names(smurf.mix)[29] <- "ciw.mix"
names(smurf.NM)[29] <- "ciw.NM"
smurf <- join(smurf.mix, smurf.NM[,c("Rep", "simdat", "ciw.NM")])
smurf$Ratio <- smurf$ciw.NM/smurf$ciw.mix
qplot(Ratio, Rep, geom="point", data=smurf) + facet_grid(datfamily~mixpeak, scale="free_x")
ddply(smurf, .(datfamily, mixpeak), summarize, M=mean(Ratio))
