setwd(".../Sims/SimZero") # Run from the TimeToDetect directory

########## Helper equations for each distribution (homogeneous surveys)
### Exponential Distribution
# First  function calculates intcpt_d to obtain a target pdet (no mixture)
# Second function calculates intcpt_d to obtain a target pdet given gamma (mixture model)
e.id.pdet    <- function(pdet) log(-log(1-pdet)/10)
emix.id.pdet <- function(pdet, gamma=0.65) e.id.pdet((pdet-(1-gamma))/gamma)

### Gamma Distribution
# First  function calculates intcpt_d & alpha to obtain a target pdet and mode for f_T(t) (peaked, no mixture)
# Second function calculates intcpt_d & alpha to obtain a target pdet and Pr(t<2)/pdet    (nonpeaked, no mixture)
# Third  function is like the first  but conditioned on gamma 
#         (technically solves the first function for intcpt_d & alpha to obtain target p(det|bird is hard-to-detect)
# Fourth function is like the second but conditioned on gamma
g.id.pdet.p <- function(pdet, mode){
  whelp <- function(alpha, pdet, mode) abs(pdet - pgamma(10, alpha, (alpha-1)/mode))
  alpha <- optimize(whelp, interval=c(0,10), pdet, mode)$minimum
  int.d <- log((alpha-1)/mode)
  return(c(int.d=int.d, alpha=alpha))
}
g.id.pdet.n <- function(pdet, p2){
  # Helper fxn: for optimizing both F_T(10) = pdet AND F_T(2) = Pr(t<2)
  whelp2 <- function(parms, pdet, p2) abs(pgamma(10,parms[1],exp(parms[2]))-pdet) +
    abs(pgamma(2,parms[1],exp(parms[2]))/pgamma(10,parms[1],exp(parms[2]))-p2)
  return(optim(c(2,-2), whelp2, gr=NULL, pdet, p2)$par)
}
gmix.id.pdet.p <- function(pdet, mode, gamma=0.65) g.id.pdet.p((pdet-(1-gamma))/gamma, mode)
gmix.id.pdet.n <- function(pdet, p2, gamma=0.65) g.id.pdet.n((pdet-(1-gamma))/gamma, 
                                                             (p2*pdet-(1-gamma))/(pdet-(1-gamma)))

### Lognormal Distribution
# First  function calculates intcpt_d & alpha to obtain a target pdet and mode for f_T(t) (peaked, no mixture)
# Second function calculates intcpt_d & alpha to obtain a target pdet and Pr(t<2)/pdet    (nonpeaked, no mixture)
# Third  function is like the first  but conditioned on gamma
# Fourth function is like the second but conditioned on gamma
l.id.pdet.p <- function(pdet, mode){
  lhelp <- function(sigma.det, pdet, mode) abs(pdet - plnorm(10,sigma.det^2+log(mode),sigma.det))
  sig.det <- optimize(lhelp, interval=c(0,10), pdet, mode)$minimum
  int.d <- -(sig.det^2 + log(mode))
  return(c(int.d=int.d, sig.det=sig.det))
}
l.id.pdet.n <- function(pdet, p2){
  lhelp2 <- function(parms, pdet, p2) abs(plnorm(2,-parms[1],parms[2])/plnorm(10,-parms[1],parms[2])-p2) + 
                                     abs(plnorm(10,-parms[1],parms[2])-pdet)
  return(optim(c(-1.5,0.75), lhelp2, gr=NULL, pdet, p2)$par)
}
lmix.id.pdet.p <- function(pdet, mode, gamma=0.65) l.id.pdet.p((pdet-(1-gamma))/gamma, mode)
lmix.id.pdet.n <- function(pdet, p2, gamma=0.65) l.id.pdet.n((pdet-(1-gamma))/gamma, 
                                                             (p2*pdet-(1-gamma))/(pdet-(1-gamma)))

### Weibull Distribution
# First  function calculates intcpt_d & alpha to obtain a target pdet and mode for f_T(t) (peaked, no mixture)
# Second function calculates intcpt_d & alpha to obtain a target pdet and Pr(t<2)/pdet    (nonpeaked, no mixture)
# Third  function is like the first  but conditioned on gamma
# Fourth function is like the second but conditioned on gamma
w.id.pdet.p <- function(pdet, mode){
  whelp <- function(shape_k, pdet, mode) abs(pdet - pweibull(10, shape_k, mode*(1-1/shape_k)^(-1/shape_k)))
  shape_k <- optimize(whelp, interval=c(0,10), pdet, mode)$minimum
  int.d <- log(1-1/shape_k)/shape_k - log(mode)
  return(c(int.d=int.d, shape_k=shape_k))
}
w.id.pdet.n <- function(pdet, p2){
  whelp2 <- function(parms, pdet, p2) abs(pweibull(10,parms[1],exp(-parms[2]))-pdet) +
    abs(pweibull(2,parms[1],exp(-parms[2]))/pweibull(10,parms[1],exp(-parms[2]))-p2)
  return(optim(c(2,-2), whelp2, gr=NULL, pdet, p2)$par)
}
wmix.id.pdet.p <- function(pdet, mode, gamma=0.65) w.id.pdet.p((pdet-(1-gamma))/gamma, mode)
wmix.id.pdet.n <- function(pdet, p2, gamma=0.65) w.id.pdet.n((pdet-(1-gamma))/gamma, 
                                                             (p2*pdet-(1-gamma))/(pdet-(1-gamma)))

########## Start constructing Simpars for the models of interest to me
### Set simulation targets
pdet.cand  <- c(0.5, 0.6, 0.8, 0.95)
gamma.cand <- c(5/6, 4/5, 0.65, 1/2)  # To get smaller pdet, gamma must be larger
cand  <- 1    # Choose which set of simulations to establish
pdet  <- pdet.cand[cand]
gamma <- gamma.cand[cand]
mode  <- 5
p2    <- 0.7  # Target Pr(t<2)/pdet for nonpeaked dist'n.  Value from Ovenbird data: p2 = 596/947    
datatype <- c("EFF","ETF","GFF","GFT","GTF","GTT","LFF","LFT","LTF","LTT","WFF","WFT","WTF","WTT")
modcode  <- c(1,2,rep(3:8,each=2))  # Numeric code associated with each datatype
Simpars  <- data.frame(matrix(NA, nrow=14, ncol=5), model=datatype, n_ints=9, modcode=modcode)
names(Simpars)[1:5] <- c("intcpt_d", "gamma", "alpha", "sigma_det", "shape_k")

# Calculate necessary parameter values to attain targets
options(warn=-1) # Suppress warnings
Simpars[1,1]         <- e.id.pdet(pdet) # EFF
Simpars[2,c(1,2)]    <- c(emix.id.pdet(pdet,gamma), gamma) # ETF
Simpars[3,c(3,1)]    <- g.id.pdet.n(pdet, p2) # GFF
Simpars[4,c(1,3)]    <- g.id.pdet.p(pdet, mode) # GFT
Simpars[5,c(3,1,2)]  <- c(gmix.id.pdet.n(pdet, p2, gamma), gamma) # GTF
Simpars[6,c(1,3,2)]  <- c(gmix.id.pdet.p(pdet, mode, gamma), gamma) # GTT
Simpars[7,c(1,4)]    <- l.id.pdet.n(pdet, p2) # LFF
Simpars[8,c(1,4)]    <- l.id.pdet.p(pdet, mode) # LFT
Simpars[9,c(1,4,2)]  <- c(lmix.id.pdet.n(pdet, p2, gamma), gamma) # LTF
Simpars[10,c(1,4,2)] <- c(lmix.id.pdet.p(pdet, mode, gamma), gamma) # LTT
Simpars[11,c(5,1)]   <- w.id.pdet.n(pdet, p2) # WFF
Simpars[12,c(1,5)]   <- w.id.pdet.p(pdet, mode) # WFT
Simpars[13,c(5,1,2)] <- c(wmix.id.pdet.n(pdet, p2, gamma), gamma) # WTF
Simpars[14,c(1,5,2)] <- c(wmix.id.pdet.p(pdet, mode, gamma), gamma) # WTT
options(warn=0) # Turn warnings back on

# par(mfrow=c(3,5)) # If we plot each distribution separately
source("../../OVEN/Compile_addendum.R")
Simpars <- plotmedians(Simpars, T)
Simpars <- Simpars[,-which(names(Simpars)=="indx")]

# Calculate the mean TTD conditional upon the bird being observed
Simpars$meanTTD <- NA
modcode <- Simpars$modcode
for(j in 1:nrow(Simpars)){
  if(modcode[j] %in% 1:2) smurf <- function(x, blank, id) x*dexp(x,exp(id))/
    pexp(10,exp(id))
  if(modcode[j] %in% 3:4) smurf <- function(x, alpha, id) x*dgamma(x,alpha,exp(id))/
    pgamma(10,alpha,exp(id))
  if(modcode[j] %in% 5:6) smurf <- function(x, sigma_det, id) x*dlnorm(x,-id,sigma_det)/
    plnorm(10,-id,sigma_det)
  if(modcode[j] %in% 7:8) smurf <- function(x, shape_k, id) x*dweibull(x,shape_k,exp(-id))/
    pweibull(10,shape_k,exp(-id))
  Simpars$meanTTD[j] <- integrate(smurf, lower=0, upper=10, Simpars[j,floor((j+9)/4)], 
                                  Simpars$intcpt_d[j])$value
}

# Add intcpt_a -- to obtain about as many counts as are in our data
Simpars$intcpt_a <- log((1/pdet)*(947/381))


save(Simpars, file="../Simpars.Rdata")
