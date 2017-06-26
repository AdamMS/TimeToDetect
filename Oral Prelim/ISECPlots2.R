#### Code to extract MCMC draws relevant to p_det
# setwd("C:/Users/Adam/Documents/Adam/IAState/Stat 599JN/TimeToDetect/OVEN")
library(rstan)

oventable   <- read.csv("oventable.csv")
fname       <- oventable$fname
modeloutput <- as.character(oventable$modeloutput)
mo          <- oventable$parmcode

for(i in 1:nrow(oventable))  {
  load(paste0(fname[i],".Rdata"))                     # File containing stanfit and other objects
  stanobject <- eval(parse(text=modeloutput[i]))      # stanfit object
  
  intcpt_a <- rstan::extract(stanobject, "intcpt_a")$intcpt_a
  intcpt_d <- rstan::extract(stanobject, "intcpt_d")$intcpt_d
  ba <- rstan::extract(stanobject, "ba")$ba
  bd <- rstan::extract(stanobject, "bd")$bd
  sigma_a <- rstan::extract(stanobject, "sigma_a")$sigma_a
  sigma_d <- rstan::extract(stanobject, "sigma_d")$sigma_d
  if(!mo[i] %% 2) {gamma <- rstan::extract(stanobject, "gamma")$gamma
  } else {gamma <- NULL}
  
  if (mo[i] %in% 1:2) shape <- NULL
  if (mo[i] %in% 3:4){
    shape <- rstan::extract(stanobject, "alpha")$alpha
  } else if (mo[i] %in% 5:6){
    shape <- rstan::extract(stanobject, "sigma_det")$sigma_det
  } else if (mo[i] %in% 7:8){
    shape <- rstan::extract(stanobject, "shape_k")$shape_k
  }
  
  uncounted <- rstan::extract(stanobject, "unobserved")$unobserved
  totN <- rstan::extract(stanobject, "totN")$totN
  
  pdet <- 1 - rowSums(uncounted) / rowSums(totN)  # Calculate pdet empirically from posterior of n_unc and N
  # Posterior predictive
  
  out <- list(intcpt_a, intcpt_d, ba, bd, gamma, shape, pdet, sigma_a, sigma_d)
  
  save(out, file=paste0("pdets/", fname[i], "_post.Rdata"))
  
  rm(list=c("stanobject", modeloutput[i]))
}





##### Generate posteriors for ISEC presentation
library(ggplot2)
library(reshape2)
setwd("../Oral Prelim/pdet_Oven")
oventable$n_ints <- NULL
oventable$mix <- c("NonMix", "Mix")

pdet <- NULL
for(j in 1:8){
  load(paste0(oventable$fname[j],"_post.Rdata"))
  if(is.null(out$gamma)) out[[5]] <- rep(NA, length(out[[1]]))
  if(is.null(out$shape)) out[[6]] <- rep(NA, length(out[[1]]))
  outtab <- NULL
  for(k in 1:9) outtab <- cbind(outtab, out[[k]], row.names=NULL)
  p.add <- data.frame(cbind(oventable[j,], outtab, row.names=NULL))
  names(p.add)[5:21] <- c("intcpt_a", "intcpt_d", paste0("ba_",1:4), paste0("bd_",1:5), "gamma", "shape", "pdet", paste0("sigma_a",1:2), "sigma_d")
  pdet <- rbind(pdet, p.add)
  rm(outtab)
}

pdet$fname <- factor(pdet$fname, labels=c("Expo", "Expo_Mix", "Gamma", "Gamma_Mix",
                                          "LogNormal", "LogNormal_Mix", "Weibull", "Weibull_Mix"))
pdet <- pdet[c(F,F,T),]    # Thinning
#pdet <- pdet[c(F,F,F,F,F,F,T,F,F,F),]

pdet <- melt(pdet, id.vars=1:4, variable.name="Parameter", value.name="Estimate")
pdet <- na.omit(pdet)

#pdf("../Posteriors_Mixtures.pdf", height=5.75, width=8.5)
print(
  ggplot(pdet, aes(Estimate, group=fname, color=fname)) + 
    geom_hline(yintercept=0) + geom_vline(xintercept=0) +
    stat_density(adjust=1, geom="line", position="identity", size = 0.75, aes(linetype=mix)) +
    scale_color_manual(name="Inference Model:", values=c("palegreen4", "palegreen3", "palegreen2", "palegreen1", "maroon4", "maroon3", "maroon2", "maroon1")) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="bottom",
          strip.text = element_text(size = 12),
          legend.title = element_text(size = 12)) +
    xlab("Estimate") + ylab("Density") +
    facet_wrap(~Parameter, scales="free") 
)
#dev.off()


load("../../OVEN/Data9.Rdata")
y <- colSums(dat$y); y <- unname(c(y[1]/2, y[2:9]))  # per minute
y2 <- NULL
for(j in 1:length(y)) y2 <- c(y2, 0, y[j], y[j]); y2 <- c(y2, 0)
tau <- c(0,2:10)
Time <- rep(tau, each=3)[-c(1,30)]
oven <- data.frame(Time, y2)
png("../Ovenbirds.png", height=4.75, width=11.5, units="in", res=96)
print(
  ggplot(oven, aes(x=Time, y=y2)) + 
    geom_hline(yintercept=0, color="palegreen3", size=1.5) + 
    geom_line(size=1.5, color="palegreen4") +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="none",
          plot.title = element_text(size=37),
          axis.title = element_text(size=28),
          axis.text.x = element_text(size=28),
          axis.text.y = element_text(size=28))  +
    xlab("Time") + ylab("Counts per Minute") + ggtitle("Ovenbirds - All Sites")
)
dev.off()
