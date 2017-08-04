########## This script uses output from StanPlot_<name>.R.  Specifically, it loads the list object 'outtable' for each model fit.  'Outtable' has three parts (four for simulated datasets): (1) The Stan summary table, (2) summary-style tables for uncounted and overall p_det, and (3) posterior predictive p-values.  This script then:
# (a) Collects parameter estimates (object: parmests) and makes caterpillar plots
#     - also does this for uncounted and overall p_det
# (b) Collects interval-by-interval posterior predictive p-values (pval3 and pvals9)

setwd("oven_sum")

##### Load pre-saved summary tables
library(plyr)
oventable <- read.csv("../oventable.csv")
fname <- oventable$fname
n_ints <- oventable$n_ints
modeloutput <- oventable$modeloutput
load("../Data9.Rdata")

##### Compile output from the various models into three data frames
parmests <- count.det <- pvals3 <- pvals9 <- dic <- NULL
for(i in 1:length(fname)){
  load(paste(fname[i],"_sum.Rdata",sep=""))
  smurf <- data.frame(outtable[[1]])
  smurf$model <- fname[i]
  smurf$n_ints <- n_ints[i]
  smurf$modeloutput <- modeloutput[i]
  smurf$parameter <- dimnames(outtable[[1]])[[1]]
  # parmests is all of the summary tables appended together plus a couple of ID columns
  parmests <- rbind(parmests, smurf)

  # Read in the posterior estimates for uncounted birds and p_global for model fname[i]
  smurf <- data.frame(outtable[[2]])
  # Convert from uncounted to total abundance
  smurf$model <- fname[i]
  smurf$parameter <- dimnames(outtable[[2]])[[1]]
  # count.det is the 'uncounted' and 'p_global' summaries appended together
  count.det <- rbind(count.det, smurf)
  
  # pvals3 and pvals9 separately append all of the poster p-values together
  if(length(outtable[[3]])==3) {pvals3 <- rbind(pvals3, outtable[[3]])
  } else {
    pvals9 <- rbind(pvals9, outtable[[3]])
  }
  
  # DIC Tables
  smurf <- c(dic=outtable[[5]][[1]], parmcode=i)
  dic <- rbind(dic, smurf)
}
dic <- join(oventable, as.data.frame(dic))
dic$dicdiff <- dic$dic - min(dic$dic)
# save(parmests, file="../parmests.Rdata")

# ### Posterior predictive plots of counts by interval for each model
# # pvals3 <- data.frame(pvals3)
# pvals9 <- data.frame(pvals9)
# # pvals3$model <- fname[oventable$n_ints==3]
# pvals9$model <- fname[oventable$n_ints!=3]
# par(mfrow=c(1,2))
# plot(c(1,(2:9)+0.5), seq(0,1,length=9), type='n', main="Posterior predictive p-values by interval")
# abline(h=c(0,0.05,0.95,1), lty=3)
# for(i in c(1,3,5,7)) lines(c(1,(2:9)+0.5), pvals9[i,1:9], col=(i+1)/2, lty=(i+1)/2)
# legend("right",legend=pvals9$model[c(1,3,5,7)], cex=0.5, lty=1:4, col=1:4)
# plot(c(1,(2:9)+0.5), seq(0,1,length=9), type='n', main="Posterior predictive p-values by interval")
# for(i in c(2,4,6,8)) lines(c(1,(2:9)+0.5), pvals9[i,1:9], col=i/2, lty=i/2)
# abline(h=c(0,0.05,0.95,1), lty=3)
# legend("right",legend=pvals9$model[c(2,4,6,8)], cex=0.5, lty=1:4, col=1:4)

# ### Summarize posterior p-values in tabular form
# library(xtable)
# intvl.counts <- colSums(dat$y)
# ICtbl <- rbind(intvl.counts, round(pvals9[,1:9],3))
# names(ICtbl) <- paste0(c(0,dat$tau[-9]),"-",dat$tau)
# levels(pvals9$model) <- c("Expo","Exponential Mix","Gamma","Gamma Mix","LogNormal","LogNormal Mix","Weibull","Weibull Mix")
# rownames(ICtbl) <- c("Count",as.character(pvals9$model))
# print(xtable(ICtbl[c(1,2,4,8,6,3,5,9,7),], digits=3), include.rownames=T)

##### Reshape and plot results
# parmests
library(reshape2)
library(ggplot2)
library(plyr)
# Melt parameter and count/pdet estimate tables
PE.melt <- melt(parmests[,-match("modeloutput",names(parmests))], id.vars=c("model","parameter"))
count.det[grepl("uncounted", rownames(count.det)),c("mean","X2.5.","X25.","X50.","X75.","X97.5.")] <- 
  sum(dat$y) + count.det[grepl("uncounted", rownames(count.det)),c("mean","X2.5.","X25.","X50.","X75.","X97.5.")]
CD.melt <- melt(count.det, id.vars=c("model","parameter"))

# # To pare away the Gamma's or lognormals for simplified plotting
# CD.melt <- CD.melt[!(CD.melt$model %in% c("Gamma3", "Gamma9")),]
# CD.melt <- CD.melt[!(CD.melt$model %in% c("LogNormal3", "LogNormal9")),]

# Plot abundcance on the log-scale
CD.melt$value[CD.melt$parameter=="uncounted"] <- log10(CD.melt$value[CD.melt$parameter=="uncounted"])
CD.melt$parameter[CD.melt$parameter=="uncounted"] <- "Abundance_log10" 

# Combine parameters and count/pdet into a single melted table; remove "lp__"
S.melt <- rbind(PE.melt, CD.melt)
S2.melt <- S.melt[S.melt$parameter!="lp__",]

### Re-labelling
# Parameter Names
# -- beware: the labels are hard-coded... twice.  If I change names, I change the alphabetical order, and I mess up the labels...
S2.melt$label <- as.factor(S2.melt$parameter)
levels(S2.melt$label) <- c("log[10](Abundance)", "alpha", "Site~Age~(A)", "Year~(A)", "Stock~Density~(A)", "Logging~(A)",
                           "Julian~Date~(D)", "Time~of~Day~(D)", "Temperature~(D)", "First~Year~(D)", "JulDate~x~FirstYr~(D)",
                           "gamma", "Intercept~(A)", "Intercept~(D)", "p^(det)", "k", "Year~Std~Dev~(A)", "Stand~Std~Dev~(A)", 
                           "Observer~Std~Dev~(D)", "sigma_det")
# Re-order the labels.  Then merge all shape parameters
S2.melt$label <- factor(S2.melt$label, levels(S2.melt$label)[c(13,3:6,17:18,14,7:11,19,12,2,16,20,15,1)])
levels(S2.melt$label)[levels(S2.melt$label) %in% c("alpha", "k", "sigma_det")] <- "alpha"
# Model Names
levels(S2.melt$model) <- c("Exponential", "Exponential Mix", "Gamma", "Gamma Mix", "Lognormal", "Lognormal Mix", "Weibull", "Weibull Mix")
S2.melt$model <- factor(S2.melt$model, levels(S2.melt$model)[c(6,8,4,2,5,7,3,1)])

# If I want only mixture models... currently redundant
S2.melt$mix <- grepl("Mix", S2.melt$model)
S2.melt <- S2.melt[S2.melt$mix,]

#pdf(paste("OVEN_Posteriors.pdf",sep=""), width=7.5, height=6.5)
jpeg(paste("OVEN_Posteriors.jpg",sep=""), width=7.5, height=6.5, units="in", res=300)
# Dot is median; diamond is mean (commented out for now)
print(qplot(value, model, geom="line", group=model, data=subset(S2.melt, variable %in% c("X2.5.","X97.5."))) + 
          facet_wrap(~label, scales="free_x", labeller=labeller(label=label_parsed,.multi_line=T), nrow=5) + 
          xlab("Posterior Parameter Estimates") + ylab("Model") + theme_bw() +
          geom_line(aes(x=value, y=model),size=I(1.15),color="gray60", # "orange",
                    data=subset(S2.melt,variable %in% c("X25.","X75."))) + 
          geom_point(aes(x=value, y=model),data=subset(S2.melt,variable %in% "X50.")) + 
          #geom_point(aes(x=value, y=model, color="blue"), shape=I(5),data=subset(S2.melt,variable %in% "mean")) +
          theme(legend.position="none", strip.text.x=element_text(size=9.5), 
                axis.text.x = element_text(size=6.5))
      )
dev.off()

##### Divide Above Plot into Two for ISEC Talk
effects <- c("ba[1]", "ba[2]", "ba[3]", "ba[4]", "bd[1]", "bd[2]", "bd[3]", "bd[4]", "bd[5]", "sigma_a[1]", 
             "sigma_a[2]", "sigma_d")
modelz <- c("Exponential Mix", "Gamma Mix", "Lognormal Mix", "Weibull Mix")
levels(S2.melt$label)[c(4,7,13:14)] <- c("Stock \n Density (A)", "Stand \n Std Dev (A)",
                                         "JulDate x\n FirstYr (D)", "Observer\n Std Dev (D)")
png("../../Oral Prelim/PosteriorOvenEffects.png", height=5.75, width=9, units="in", res=96)
# Dot is median; diamond is mean (commented out for now)
print(
  qplot(value, model, geom="line", group=model, size=I(1.8),
        data=subset(S2.melt, variable %in% c("X2.5.","X97.5.") & parameter %in% effects & model %in% modelz)) +
    geom_vline(xintercept=0, linetype=2, color="deepskyblue2", size=1.25)+ 
    facet_wrap(~label, scales="free_x", nrow=2) + xlab("Parameter Estimate") + 
    ylab("Detection Parameters                Abundance Parameters\n Model") + theme_bw() +
    geom_line(aes(x=value, y=model),size=I(2.35),color="orange",
              data=subset(S2.melt, variable %in% c("X25.","X75.") & parameter %in% effects & model %in% modelz)) + 
    geom_point(aes(x=value, y=model), size=I(2.5), 
               data=subset(S2.melt, variable %in% "X50." & parameter %in% effects & model %in% modelz)) + 
    #geom_point(aes(x=value, y=model, color="blue"), shape=I(5),data=subset(S2.melt,variable %in% "mean")) +
    theme(legend.position="none", strip.text.x=element_text(size=11), 
          panel.grid.major = element_blank(), axis.text.y=element_text(size=12.5), 
          panel.grid.minor = element_blank(), axis.title=element_text(size=14.5))
)
dev.off()

png("../../Oral Prelim/PosteriorOvenEffects_abund.png", height=2.5, width=9.5, units="in", res=96)
# Dot is median; diamond is mean (commented out for now)
S2.A <- subset(S2.melt, parameter %in% effects & model %in% modelz & label %in% c("Site Age (A)", "Year (A)", "Logging (A)", "Stand \n Std Dev (A)"))
print(
  qplot(value, model, geom="line", group=model, size=I(1.8), data=subset(S2.A, variable %in% c("X2.5.","X97.5."))) +
    geom_vline(xintercept=0, linetype=2, color="deepskyblue2", size=1.25)+ 
    facet_wrap(~label, scales="free_x", nrow=1) + xlab("Posterior Parameter Estimate") + 
    ylab("Model") + theme_bw() +
    geom_line(aes(x=value, y=model), size=I(2.35), color="orange",
              data=subset(S2.A, variable %in% c("X25.","X75."))) + 
    geom_point(aes(x=value, y=model), size=I(2.5), 
               data=subset(S2.A, variable %in% "X50.")) + 
    theme(legend.position="none", strip.text.x=element_text(size=15), plot.title=element_text(size=24), 
          panel.grid.major = element_blank(), axis.text.y=element_text(size=15), 
          panel.grid.minor = element_blank(), axis.title=element_text(size=17.5)) 
    # + ggtitle("Posterior Abundance Effect Estimates")
)
dev.off()
png("../../Oral Prelim/PosteriorOvenEffects_det.png", height=2.5, width=9.5, units="in", res=96)
# Dot is median; diamond is mean (commented out for now)
S2.D <- subset(S2.melt, parameter %in% effects & model %in% modelz & label %in% c("Julian Date (D)", "Time of Day (D)", "Temperature (D)", "Observer\n Std Dev (D)"))
print(
  qplot(value, model, geom="line", group=model, size=I(1.8), data=subset(S2.D, variable %in% c("X2.5.","X97.5."))) +
    geom_vline(xintercept=0, linetype=2, color="deepskyblue2", size=1.25)+ 
    facet_wrap(~label, scales="free_x", nrow=1) + xlab("Posterior Parameter Estimate") + 
    ylab("Model") + theme_bw() +
    geom_line(aes(x=value, y=model), size=I(2.35), color="orange",
              data=subset(S2.D, variable %in% c("X25.","X75."))) + 
    geom_point(aes(x=value, y=model), size=I(2.5), 
               data=subset(S2.D, variable %in% "X50.")) + 
    theme(legend.position="none", strip.text.x=element_text(size=15), 
          panel.grid.major = element_blank(), axis.text.y=element_text(size=15), plot.title=element_text(size=24), 
          panel.grid.minor = element_blank(), axis.title=element_text(size=17.5))
    # + ggtitle("Posterior Detection Effect Estimates")
)
dev.off()

vline <- data.frame(label=c("gamma", "p^(det)", "gamma", "p^(det)"), Z = c(0,0,1,1))
png("../../Oral Prelim/OVEN_Posteriors_detabund.png", height=3.75, width=8.75, units="in", res=96)
print(
  qplot(value, model, geom="line", group=model, size=I(1.8),
        data=subset(S2.melt, variable %in% c("X2.5.","X97.5.") & 
                      parameter %in% c("gamma", "Abundance_log10", "p_global") & model %in% modelz)) +
    geom_vline(data=vline, aes(xintercept=Z), linetype=2, color="deepskyblue2", size=1.25)+ 
    facet_wrap(~label, scales="free_x", nrow=1, labeller=label_parsed) + xlab("Parameter Estimate") + ylab("Model") + theme_bw() +
    geom_line(aes(x=value, y=model),size=I(2.35),color="orange",
              subset(S2.melt, variable %in% c("X25.","X75.") & 
                       parameter %in% c("gamma", "Abundance_log10", "p_global") & model %in% modelz)) + 
    geom_point(aes(x=value, y=model), size=I(2.5), data=subset(S2.melt, variable == "X50." & 
                              parameter %in% c("gamma", "Abundance_log10", "p_global") & model %in% modelz)) + 
    #geom_point(aes(x=value, y=model, color="blue"), shape=I(5),data=subset(S2.melt,variable %in% "mean")) +
    theme(legend.position="none", strip.text.x=element_text(size=17), 
          panel.grid.major = element_blank(), axis.text.y=element_text(size=17), 
          panel.grid.minor = element_blank(), axis.title=element_text(size=17),
          plot.title=element_text(size=22), axis.text.x=element_text(size=13))) #+
#    ggtitle(bquote("Posteriors: Mixing, Pr(Detection), log" [10] (Abundance))))
dev.off()




##### Plot TTDDs at median posterior values

# ##### Um... this ignores the role of random effects on expected values...

# source("../Compile_addendum.R")
# parmvals <- parmests[,c("X50.","model","n_ints","parameter")]
# parmvals <- dcast(parmvals, model + n_ints ~ parameter, value.var="X50.")
# parmvals$modcode <- 1:8
# parmvals$indx <- 1     # It's a necessary input for plotmedians but meaningless for a single-rep sim study
# parmvals$intcpt_d <- parmvals$intcpt_d + parmvals[,match("ba[3]",names(parmvals))] # Set I(first-year)=1, since most are
# curve(dexp(x,0,0.5), xlim=c(0,10), ylim=c(0,0.35), type='n', ylab="Density")
# plotmedians(parmvals, adz=T)
# legend("topright", legend=parmvals$model, lty=c(1,2,1,3,1,4,1,5), col=c(2,3,2,4,2,1,2,6))
# # NEED TO CHANGE COLORS AND ADD A LEGEND!




# # Summary of R.hat values
# R.hat <- PE.melt[PE.melt$variable=="Rhat",]
# R.hat[R.hat$value > 1.1,]  # All model paramters with inadequate R.hat

# Summary of sample sizes --- outdated.  Go see Geweke.R
# Neff <- PE.melt[PE.melt$variable=="n_eff" & PE.melt$parameter!="lp__",]
# Neff[Neff$value < 1000,]  # All model paramters with inadequate R.hat
# Neff[order(Neff$value, decreasing=T),]

##### WIDTHS

widths <- parmests[parmests$parameter != "lp__",names(parmests) %in% c("mean", "sd", "X2.5.", "X25.", "X75.", "X97.5.", "model", "parameter")]
widths$CI95 <- widths$X97.5. - widths$X2.5.
widths$CI50 <- widths$X75. - widths$X25.
par(mfrow=c(4,6), mar = c(4, 4, 0.5, 0.6) + 0.1, oma = c(0,0,2,0), mgp=c(2,0.8,0))
# X-Axis order:
# widths$model[widths$parameter=="lp__"]
# [1] "Expo3"         "Expo9"         "ExpoMix3"      "ExpoMix9"      "Gamma3"        "Gamma9"    
# [7] "GammaMix3"     "GammaMix9"     "LogNormal3"    "LogNormal9"    "LogNormalMix3" "LogNormalMix9"
# [13] "Weibull3"     "Weibull9"      "WeibullMix3"   "Weibull9"

# Plot of 95% CI widths for each model - for root parameters only; no p_det, no uncounted
for (varble in unique(widths$parameter)){
  print(widths[widths$parameter==varble, names(widths) %in% c("sd", "CI95", "CI50", "model")])
  plot(widths$CI95[widths$parameter==varble], ylab=varble, ylim=c(0,max(widths$CI95[widths$parameter==varble],
                                              3.92*widths$sd[widths$parameter==varble])), pch=16)
  # points(3.92*widths$sd[widths$parameter==varble], col="orange")  # Plot scaled SD instead of CI-widths
}
par(mfrow=c(1,1))
