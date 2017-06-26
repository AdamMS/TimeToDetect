########## Generate figures of mixture expo and mixture gamma distributions to illustrate
# how these models fit the data

setwd("C:/Users/Adam/Documents/Adam/IAState/Stat 599JN/TimeToDetect/Oral Prelim")
#load("../OVEN/Data9.Rdata")
load("AMRO.Rdata") # There's a sharper contrast when we use AMRO

# Read in and restructure data
counts    <- colSums(dat$y)
cntbyint  <- counts
counts[1] <- counts[1]/2
names(counts) <- paste0(c(0,2:9)," - ", 2:10)

### MLEs for TTDDs

# Exponential MLE -- outputs estimate as (phi, gamma)
phi8.llik <- function(parameters, cnt, tau=c(0,2:10)){
  ints  <- length(cnt)
  phi   <- parameters[1]
  gamma <- parameters[2]
  llk <- cnt[1] * log(1-gamma + gamma * (pexp(tau[2],phi) - pexp(tau[1],phi)))
  for (j in 2:ints) llk <- llk + cnt[j]*log(gamma * (pexp(tau[j+1],phi) - pexp(tau[j],phi)))
  llk <- llk - sum(cnt)*log(1-gamma + gamma*pexp(max(tau),phi))
  return(unname(llk))
}
efit <- optim(c(0.13, 0.6), fn=phi8.llik, gr=NULL, cnt=cntbyint, control=list(fnscale=-1))$par
p.e <- efit[2]*pexp(10,efit[1]) + (1-efit[2])

# Gamma MLE -- outputs estimate as (alpha, phi, gamma)
# Ignore the warnings
gamma.llik <- function(parameters, cnt, tau=c(0,2:10)){
  ints  <- length(cnt)
  alpha <- parameters[1]
  phi   <- parameters[2]
  gamma <- parameters[3]
  llik <- cnt[1] * log(1-gamma + gamma * (pgamma(tau[2],alpha,phi) - pgamma(tau[1],alpha,phi)))
  for(j in 2:ints) llik <- llik + 
      cnt[j] * log(gamma * (pgamma(tau[j+1],alpha,phi) - pgamma(tau[j],alpha,phi)))
  llik <- llik - sum(cnt)*log(1-gamma + gamma*pgamma(max(tau),alpha,phi))
  return(unname(llik))
}
# Fit model to data
gfit <- optim(c(1,0.3,0.5), fn=gamma.llik, gr=NULL, cnt=cntbyint, control=list(fnscale=-1))$par
p.g <- gfit[3]*pgamma(10,gfit[1],gfit[2]) + (1-gfit[3])

ntot <- sum(cntbyint)    # Total birds observed

### Barplot of the data
pdf(file="AMROobs.pdf", width=5.25, height=3.25)
par(mfrow=c(1,1))
barplot(counts, width=c(2,rep(1,8)), space=0, xlab="Time to first detection", ylab="Counts per Minute", col="gray75", 
        #main="American Robins - 159 birds", cex.main=1.3, 
        cex.names=1.25, cex.lab=1.2, xaxt='n', yaxt='n')
axis(1, at=c(1,(2:9)+0.5), labels=FALSE)
text(x=c(1,(2:9)+0.6), y=par("usr")[3] - 4.5, labels=names(counts), srt=45, adj=1, xpd=T, cex=1.25)
axis(2, at=seq(0,35,5), cex=1.25)
dev.off()

# # Reproduce above barplot on proportion scale with: (i) TTDD curve fit and (ii) mixture
# # Easy-to-detect
# rect(0, efit[2]*pexp(2,efit[1])/2/p.e, 2, 
#      efit[2]*pexp(2,efit[1])/2/p.e + (1-efit[2])/p.e/2, col="palegreen3")
# # Hard-to-detect
# curve(efit[2]*dexp(x, efit[1])/p.e, add=T, lwd=2, col="blue")
# ### Now with the gamma distribution
# barplot(counts/ntot, width=c(2,rep(1,8)), space=0, xlab="Time to first detection", col="palegreen4", font.lab=2, cex.lab=1.6)
# # Easy-to-detect
# rect(0, gfit[3]*pgamma(2,gfit[1],gfit[2])/2/p.g, 2, 
#      gfit[3]*pgamma(2,gfit[1],gfit[2])/2/p.g + (1-gfit[3])/p.g/2, col="palegreen3")
# # Hard-to-detect
# curve(gfit[3]*dgamma(x,gfit[1],gfit[2])/p.g, add=T, lwd=2, col="blue")

### Exponential
# Hard-to-detect
pdf(file="AMROfit.pdf", width=7.5, height=3.15)
par(mfrow=c(1,2))
curve(efit[2]*dexp(x, efit[1])/p.e, xlim=c(0,10), ylim=c(0,0.32),lwd=2, col="palegreen4",
      main="Constant Detection", ylab="Density", xlab="Time")
# Easy-to-detect
top.e <- efit[2]*pexp(2,efit[1])/2/p.e + (1-efit[2])/p.e/2
segments(0, efit[2]*dexp(0, efit[1])/p.e, 0, top.e, lwd=2, col="palegreen3")
segments(0, top.e, 2, top.e, lwd=2, col="palegreen3")
segments(2, top.e, 2, efit[2]*dexp(2, efit[1])/p.e, lwd=2, col="palegreen3")
text(7,0.25,paste("Pr(det) =",round(p.e,2)), font=2)

### Gamma
# Hard-to-detect
curve(gfit[2]*dgamma(x, gfit[1], gfit[2])/p.g, xlim=c(0,10), ylim=c(0,0.32),lwd=2, 
      col="palegreen4", main="Non-constant Detection", ylab="Density", xlab="Time")
# Easy-to-detect
top.g <- gfit[3]*pgamma(2,gfit[1],gfit[2])/2/p.g + (1-gfit[3])/p.g/2
segments(0, gfit[3]*dgamma(0, gfit[1], gfit[2])/p.g, 0, top.g, lwd=2, col="palegreen3")
segments(0, top.g, 2, top.g, lwd=2, col="palegreen3")
segments(2, top.g, 2, gfit[2]*dgamma(2, gfit[1], gfit[2])/p.g, lwd=2, col="palegreen3")
text(7,0.25,paste("Pr(det) =",round(p.g,2)), font=2)
dev.off()




### Generate similar plots of continuous and interval-censored gamma distributions
# Using the above example.  Find same code using Ovenbird estimates below
png(file="../../../Dissertation/figs/GammaTTDDex_Dissertation.pdf", width=540, height=230)
# efit <- c(phi, gamma); gfit <- c(alpha, phi, gamma)
gamma <- gfit[3]; alpha <- gfit[1]; beta <- gfit[2]
# Hard-to-detect
p.g <- (1-gamma) + gamma*pgamma(10, alpha, beta)
par(mfrow=c(1,2))
plot(0,0,type='n', xlim=c(0,15), ylim=c(0,0.25), main="Mixture Gamma TTDD (Continuous)", ylab="Density", xlab="Time")
# Easy-to-detect
# I've decided to NOT scale to the conditional probabilities...
p.g <- 1  # If I want to scale, then reset this to its proper values
top.g <- gamma*pgamma(3, alpha, beta)/3/p.g + (1-gamma)/p.g/3
# segments(0, 0, 0, 0.24, lwd=1, lty=2, col="black")
segments(3, 0, 3, gamma*dgamma(3, alpha, beta)/p.g, lwd=1, lty=2, col="black")
segments(5, 0, 5, gamma*dgamma(5, alpha, beta)/p.g, lwd=1, lty=2, col="black")
segments(10, 0, 10, gamma*dgamma(10, alpha, beta)/p.g, lwd=1, lty=2, col="black")
segments(0, 0, 15, 0, lwd=2, col="palegreen4")  # not drawn, b/c dgamma(0,...) = Inf
segments(0, gamma*dgamma(0, alpha, beta)/p.g, 0, top.g, lwd=2, col="palegreen3")  # not drawn, b/c dgamma(0,...) = Inf
segments(0.0846, top.g, 3, top.g, lwd=2, col="palegreen3")   # 0.1868 calculated by optimization
segments(3, top.g, 3, gamma*dgamma(3, alpha, beta)/p.g, lwd=2, col="palegreen3")
curve(gamma*dgamma(x, alpha, beta)/p.g, n=1001, lwd=2, col="palegreen4", add=T)

# Interval-censored
# p.03 <- gamma*(pgamma(3, alpha, beta) - pgamma(0, alpha, beta)) / p.g / 3
# p.35 <- gamma*(pgamma(5, alpha, beta) - pgamma(3, alpha, beta)) / p.g / 2
# p.510 <- gamma*(pgamma(10, alpha, beta) - pgamma(5, alpha, beta)) / p.g / 5
# plot(c(0,0,3,3), c(p.03,top.g,top.g,p.03), type='l', lwd=2,
#      col="palegreen3", xlim=c(0,10), ylim=c(0,0.4), main="Mixture Gamma TTDD (Censored)", ylab="Density", xlab="Time")
# lines(c(0,0,3,3,3,5,5,5,10,10,0), c(0,p.03,p.03,0,p.35,p.35,0,p.510,p.510,0,0), lwd=2, col="palegreen4")

p.03     <- gamma*(pgamma(3, alpha, beta) - pgamma(0, alpha, beta)) / p.g
p.35     <- gamma*(pgamma(5, alpha, beta) - pgamma(3, alpha, beta)) / p.g
p.510    <- gamma*(pgamma(10, alpha, beta) - pgamma(5, alpha, beta)) / p.g
p.10plus <- gamma*(1 - pgamma(10, alpha, beta)) / p.g
ps <- c(p.03, p.35, p.510, p.10plus)
plot(1:5, 1:5, type='n', xlim=c(0.5,4.5), ylim=c(0,p.03+(1-gamma)+0.05), xaxt='n',
     main="Mixture Gamma TTDD (Censored)", ylab="Density", xlab="Observation Interval")
segments(1,p.03,1,p.03+(1-gamma),lwd=8,col="palegreen3")
for(j in 1:4) segments(j, 0, j, ps[j], lwd=8, col="palegreen4")
axis(1, at=1:4, labels=c("0-3", "3-5", "5-10", "Not \n detected"), cex.axis=0.8)
dev.off()

##### The original version!!! Uses actual Ovenbird mean posterior estimates
# pdf(file="GammaTTDDex.pdf", width=7.5, height=3.15)
# gamma <- 0.777; alpha <- 0.727; beta <- exp(-2.914)
# # Hard-to-detect
# p.g <- (1-gamma) + gamma*pgamma(10, alpha, beta)
# par(mfrow=c(1,2))
# plot(0,0,type='n', xlim=c(0,15), ylim=c(0,0.25), main="Mixture Gamma TTDD (Continuous)", ylab="Density", xlab="Time")
# # Easy-to-detect
# # I've decided to NOT scale to the conditional probabilities...
# p.g <- 1  # If I want to scale, then reset this to its proper values
# top.g <- gamma*pgamma(3, alpha, beta)/3/p.g + (1-gamma)/p.g/3
# segments(0, 0, 0, 0.24, lwd=1, lty=2, col="black")
# segments(3, 0, 3, gamma*dgamma(3, alpha, beta)/p.g, lwd=1, lty=2, col="black")
# segments(5, 0, 5, gamma*dgamma(5, alpha, beta)/p.g, lwd=1, lty=2, col="black")
# segments(10, 0, 10, gamma*dgamma(10, alpha, beta)/p.g, lwd=1, lty=2, col="black")
# segments(0, 0, 15, 0, lwd=2, col="palegreen4")  # not drawn, b/c dgamma(0,...) = Inf
# segments(0, gamma*dgamma(0, alpha, beta)/p.g, 0, top.g, lwd=2, col="palegreen3")  # not drawn, b/c dgamma(0,...) = Inf
# segments(0.0846, top.g, 3, top.g, lwd=2, col="palegreen3")   # 0.1868 calculated by optimization
# segments(3, top.g, 3, gamma*dgamma(3, alpha, beta)/p.g, lwd=2, col="palegreen3")
# curve(gamma*dgamma(x, alpha, beta)/p.g, n=1001, lwd=2, col="palegreen4", add=T)
# 
# # Interval-censored
# # p.03 <- gamma*(pgamma(3, alpha, beta) - pgamma(0, alpha, beta)) / p.g / 3
# # p.35 <- gamma*(pgamma(5, alpha, beta) - pgamma(3, alpha, beta)) / p.g / 2
# # p.510 <- gamma*(pgamma(10, alpha, beta) - pgamma(5, alpha, beta)) / p.g / 5
# # plot(c(0,0,3,3), c(p.03,top.g,top.g,p.03), type='l', lwd=2,
# #      col="palegreen3", xlim=c(0,10), ylim=c(0,0.4), main="Mixture Gamma TTDD (Censored)", ylab="Density", xlab="Time")
# # lines(c(0,0,3,3,3,5,5,5,10,10,0), c(0,p.03,p.03,0,p.35,p.35,0,p.510,p.510,0,0), lwd=2, col="palegreen4")
# 
# p.03     <- gamma*(pgamma(3, alpha, beta) - pgamma(0, alpha, beta)) / p.g
# p.35     <- gamma*(pgamma(5, alpha, beta) - pgamma(3, alpha, beta)) / p.g
# p.510    <- gamma*(pgamma(10, alpha, beta) - pgamma(5, alpha, beta)) / p.g
# p.10plus <- gamma*(1 - pgamma(10, alpha, beta)) / p.g
# ps <- c(p.03, p.35, p.510, p.10plus)
# plot(1:5, 1:5, type='n', xlim=c(0.5,4.5), ylim=c(0,p.03+(1-gamma)+0.05), xaxt='n',
#      main="Mixture Gamma TTDD (Censored)", ylab="Density", xlab="Observation Interval")
# segments(1,p.03,1,p.03+(1-gamma),lwd=8,col="palegreen3")
# for(j in 1:4) segments(j, 0, j, ps[j], lwd=8, col="palegreen4")
# axis(1, at=1:4, labels=c("0-3", "3-5", "5-10", "Not \n detected"), cex.axis=0.8)
# dev.off()



########## ASA Slides

pdf(file="AMRO_ASA1.pdf", width=11.25, height=4.25)
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.1, mgp=c(2.3,1,0))
curve(efit[2]*dexp(x, efit[1])/p.e, xlim=c(0,10), ylim=c(0,0.32), lwd=2, col="palegreen4",
      main="Constant Detection", ylab="Density", xlab="Time", xaxt='n', yaxt='n', cex.lab=2.5, cex.main=2.5)
# Easy-to-detect
top.e <- efit[2]*pexp(2,efit[1])/2/p.e + (1-efit[2])/p.e/2
segments(0, efit[2]*dexp(0, efit[1])/p.e, 0, top.e, lwd=2, col="palegreen3")
segments(0, top.e, 2, top.e, lwd=2, col="palegreen3")
segments(2, top.e, 2, efit[2]*dexp(2, efit[1])/p.e, lwd=2, col="palegreen3")
text(7,0.25,paste("Pr(det) =",round(p.e,2)), font=2, cex=2)

# Plot data
barplot(counts, width=c(2,rep(1,8)), space=0, xlab="Time to first detection", ylab="Counts per Minute", 
        col="palegreen4", main="American Robins \n 159 birds", cex.main=2.5, cex.names=2.5, cex.lab=2.5, 
        xaxt='n', yaxt='n')
# axis(1, at=c(1,(2:9)+0.5), labels=FALSE)
# text(x=c(1,(2:9)+0.5), y=par("usr")[3] - 1.65, labels=names(counts), srt=45, adj=1, xpd=T, cex=2.5)
# axis(2, at=seq(0,35,5), cex=2.5)
dev.off()




pdf(file="AMRO_ASA2.pdf", width=11.25, height=4.25)
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.1, mgp=c(2.3,1,0))
curve(efit[2]*dexp(x, efit[1])/p.e, xlim=c(0,10), ylim=c(0,0.32), lwd=2, col="palegreen4",
      main="Constant Detection", ylab="Density", xlab="Time", xaxt='n', yaxt='n', cex.lab=2.5, cex.main=2.5)
# Easy-to-detect
top.e <- efit[2]*pexp(2,efit[1])/2/p.e + (1-efit[2])/p.e/2
segments(0, efit[2]*dexp(0, efit[1])/p.e, 0, top.e, lwd=2, col="palegreen3")
segments(0, top.e, 2, top.e, lwd=2, col="palegreen3")
segments(2, top.e, 2, efit[2]*dexp(2, efit[1])/p.e, lwd=2, col="palegreen3")
text(7,0.25,paste("Pr(det) =",round(p.e,2)), font=2, cex=2)

# Plot data
barplot(counts, width=c(2,rep(1,8)), space=0, xlab="Time to first detection", ylab="Counts per Minute", 
        col="palegreen4", main="American Robins \n 159 birds", cex.main=2.5, cex.names=2.5, cex.lab=2.5, 
        xaxt='n', yaxt='n')

curve(gfit[2]*dgamma(x, gfit[1], gfit[2])/p.g, xlim=c(0,10), ylim=c(0,0.32), lwd=2, col="palegreen4",
      main="Non-constant Detection", ylab="Density", xlab="Time", xaxt='n', yaxt='n', cex.lab=2.5, cex.main=2.5)
# Easy-to-detect
top.g <- gfit[3]*pgamma(2,gfit[1],gfit[2])/2/p.g + (1-gfit[3])/p.g/2
segments(0, gfit[3]*dgamma(0, gfit[1], gfit[2])/p.g, 0, top.g, lwd=2, col="palegreen3")
segments(0, top.g, 2, top.g, lwd=2, col="palegreen3")
segments(2, top.g, 2, gfit[2]*dgamma(2, gfit[1], gfit[2])/p.g, lwd=2, col="palegreen3")
text(7,0.25,paste("Pr(det) =",round(p.g,2)), font=2, cex=2)
dev.off()



##### Modify above plots for JABES
pdf(file="AMRO_JABES.pdf", width=8.25, height=3.15)
par(mfrow=c(1,3), mar=c(5,4,4,2)+0.1, mgp=c(2.3,1,0))
tauplot <- c(0,2:10,15)
### Exponential
plot(0,0,xlim=c(0,15),ylim=c(0,0.22), type='n', ann=F, xaxt='n', yaxt='n')
mtext(side = 1, text = "Time           ", line = 2.4, cex=1.2)
mtext(side = 2, text = "Density", line = 1, cex=1.2)
axis(1, at=c(0,10), labels=c(0,"C"), cex.axis=1.2)
# Hard to detect mixture
rect(0, efit[2]*pexp(2,efit[1])/2, 2, 
     efit[2]*pexp(2,efit[1])/2 + (1-efit[2])/p.e/2, col="gray90")
# Multinomial p_i with heights adjusted for interval width
for(j in 1:9) rect(tauplot[j], 0, tauplot[j+1], 
                   efit[2]*(pexp(tauplot[j+1],efit[1])-pexp(tauplot[j],efit[1])) / (tauplot[j+1]-tauplot[j]), col="gray65")
# Easy-to-detect continuous-time TTDD
curve(efit[2]*dexp(x, efit[1]), lwd=3, col="gray10", add=T, xlim=c(0,22))
text(10,0.135,paste("p(det) =",round(p.e,2)), font=2, cex=1.6)
text(12,0.065,"1-p(det)", cex=1.6)
segments(-2,0,22,0) # x-axis
segments(11.25,0.015,12.5,0.045) # Pointer line

# Plot data
# barplot(counts, width=c(2,rep(1,8)), space=0, xlab=NULL, ylab="Counts per Minute", col="gray65", 
#         #main="American Robins - 159 birds", cex.main=1.3, 
#         cex.names=1.25, cex.lab=1.6, xaxt='n', yaxt='n', ylim=c(0,40))
# axis(1, at=c(1,(2:9)+0.5), labels=FALSE)
# text(x=c(1,(2:9)+0.6), y=par("usr")[3] - 4.5, labels=names(counts), srt=45, adj=1, xpd=T, cex=1.25)
barplot(counts, width=c(2,rep(1,8)), space=0, xlab=NULL, ylab="Counts per Minute", col="gray65", 
        #main="American Robins - 159 birds", cex.main=1.3, 
        cex.names=1.25, cex.lab=1.6, xaxt='n', yaxt='n', ylim=c(0,40))
mtext(side = 1, text = "Time", line = 2.4, cex=1.2)
axis(1, at=c(0,(2:10)), labels=FALSE)
text(x=c(0.25, (2:9)+0.25, 10.5), y=par("usr")[3] - 3.5, labels=c(0,2:10), adj=1, xpd=T, cex=1.3)
axis(2, at=seq(0,35,5), cex=1.25)

### Gamma
plot(0,0,xlim=c(0,15),ylim=c(0,0.22), type='n', ann=F, xaxt='n', yaxt='n')
mtext(side = 1, text = "Time           ", line = 2.4, cex=1.2)
mtext(side = 2, text = "Density", line = 1, cex=1.2)
axis(1, at=c(0,10), labels=c(0,"C"), cex.axis=1.2)
# Hard to detect mixture
rect(0, gfit[3]*pgamma(2,gfit[1],gfit[2])/2, 2, 
     gfit[3]*pgamma(2,gfit[1],gfit[2])/2 + (1-gfit[2])/p.g/2, col="gray90")
# Multinomial p_i with heights adjusted for interval width
for(j in 1:9) rect(tauplot[j], 0, tauplot[j+1], 
                   gfit[3]*(pgamma(tauplot[j+1],gfit[1],gfit[2])-pgamma(tauplot[j],gfit[1],gfit[2])) / (tauplot[j+1]-tauplot[j]), col="gray65")
# Easy-to-detect continuous-time TTDD
curve(gfit[3]*dgamma(x, gfit[1], gfit[2]), lwd=3, col="gray10", add=T, xlim=c(0,22))
text(10,0.135,paste("p(det) =",round(p.g,2)), font=2, cex=1.6)
text(5,0.185,expression(paste("1 - ", gamma)), cex=1.6)
segments(-2,0,22,0) # x-axis
segments(1,0.18,3,0.185)
dev.off()