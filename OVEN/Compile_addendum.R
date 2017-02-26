### Dual-purpose function (why did I write it this way?)
# 'simvals' is a data.frame with appropriately named columns for known parameter values
# 'adz' is a graphical parameter indicating whether to overlay the density plot... currently disabled
# Plots interval-censored (discrete) distributions of f(t|det)
# AND
# Modifies the data.frame ?!?

# Trick: if I only want mixtures, change col=j to col=j/2
plotmedians <- function(simvals, adz=T){
  taus <- c(0,2:10)   # Interval cutoffs
  # For plotting intervals... duplicates all interior cutoffs:
  ends <- c(taus[1], rep(taus[-c(1,length(taus))],each=2), taus[length(taus)])
  ys <- numeric(length(ends))
  modcode <- simvals$modcode
  simvals$pdet <- NA
  simvals$indx <- simvals$modcode%%2
  par(mfrow=c(4,4))
  # Loop row-by-row through the table
  # Calculate pdet from the appropriate distribution (identified by 'modcode' variable)
  for(j in 1:nrow(simvals)){
    plot(1, 1, type='n', xlim=c(0,10), ylim=c(0,0.5), ylab=paste0("f(t|det): ",simvals$model[j]), xlab="Time")
    if(modcode[j] == 1) {
      simvals$pdet[j] <- pexp(10,exp(simvals$intcpt_d[j]))
      ys <- rep((pexp(taus[-1],exp(simvals$intcpt_d[j])) - pexp(taus[-length(taus)],exp(simvals$intcpt_d[j])))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*ys[1:2]
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="Expo")
#       if(simvals$n_ints[j]==9) curve(dexp(x,exp(simvals$intcpt_d[j]))/simvals$pdet[j],
#                                        add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="Expo")
    }
    if(modcode[j] == 2){
      simvals$pdet[j] <- 1-simvals$gamma[j]*(1-pexp(10,exp(simvals$intcpt_d[j])))
      ys <- simvals$gamma[j]*rep((pexp(taus[-1],exp(simvals$intcpt_d[j])) - 
                                    pexp(taus[-length(taus)],exp(simvals$intcpt_d[j])))/simvals$pdet[j] , each=2)
      ys[1:2] <- 0.5*(ys[1:2] + (1 - simvals$gamma[j])/simvals$pdet[j])
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{3}, lty=2, ylab="ExpoMix")
#       if(simvals$n_ints[j]==9) curve(((x<=2)*(1-simvals$gamma[j]*(1-pexp(2,exp(simvals$intcpt_d[j]))))/2 +
#                                        (x>2)*simvals$gamma[j]*dexp(x,exp(simvals$intcpt_d[j])))/
#                                        simvals$pdet[j], 
#                                        add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{3}, lty=2, ylab="ExpoMix")
    }
    if(modcode[j] == 3){
      simvals$pdet[j] <- pgamma(10,simvals$alpha[j],exp(simvals$intcpt_d[j]))
      ys <- rep((pgamma(taus[-1],simvals$alpha[j],exp(simvals$intcpt_d[j])) - 
                   pgamma(taus[-length(taus)],simvals$alpha[j],exp(simvals$intcpt_d[j])))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*ys[1:2]
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="Gamma")
#       if(simvals$n_ints[j]==9) curve(dgamma(x,simvals$alpha[j],exp(simvals$intcpt_d[j]))/simvals$pdet[j], 
#                                      add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="Gamma")
    }
    if(modcode[j] == 4){
      simvals$pdet[j] <- 1-simvals$gamma[j]*(1-pgamma(10,simvals$alpha[j],exp(simvals$intcpt_d[j])))
      ys <- simvals$gamma[j]*rep((pgamma(taus[-1],simvals$alpha[j],exp(simvals$intcpt_d[j])) - 
                                    pgamma(taus[-length(taus)],simvals$alpha[j],exp(simvals$intcpt_d[j])))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*(ys[1:2] + (1 - simvals$gamma[j])/simvals$pdet[j])
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{4}, lty=3, ylab="GammaMix")
#       if(simvals$n_ints[j]==9) curve(((x<=2)*(1-simvals$gamma[j]*(1-dgamma(2,simvals$alpha[j],exp(simvals$intcpt_d[j]))))/2 +
#                                        (x>2)*simvals$gamma[j]*dgamma(x,simvals$alpha[j],exp(simvals$intcpt_d[j])))/
#                                        simvals$pdet[j], 
#                                        add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{4}, lty=3, ylab="GammaMix")
    }
    if(modcode[j] == 5){
      simvals$pdet[j] <- plnorm(10,-simvals$intcpt_d[j],simvals$sigma_det[j])
      ys <- rep((plnorm(taus[-1],-simvals$intcpt_d[j],simvals$sigma_det[j]) - 
                   plnorm(taus[-length(taus)],-simvals$intcpt_d[j],simvals$sigma_det[j]))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*ys[1:2]
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="LogN")
#       if(simvals$n_ints[j]==9) curve(dlnorm(x,-simvals$intcpt_d[j],simvals$sigma_det[j])/simvals$pdet[j], 
#                                      add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="LogN")
    }
    if(modcode[j] == 6){
      simvals$pdet[j] <- 1-simvals$gamma[j]*(1-plnorm(10,-simvals$intcpt_d[j],simvals$sigma_det[j]))
      ys <- simvals$gamma[j]*rep((plnorm(taus[-1],-simvals$intcpt_d[j],simvals$sigma_det[j]) - 
                                    plnorm(taus[-length(taus)],-simvals$intcpt_d[j],simvals$sigma_det[j]))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*(ys[1:2] + (1 - simvals$gamma[j])/simvals$pdet[j])
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{1}, lty=4, ylab="LogNMix")
#       if(simvals$n_ints[j]==9) curve(((x<=2)*(1-simvals$gamma[j]*(1-plnorm(2,-simvals$intcpt_d[j],simvals$sigma_det[j])))/2 +
#                                        (x>2)*simvals$gamma[j]*dlnorm(x,-simvals$intcpt_d[j],simvals$sigma_det[j]))/
#                                        simvals$pdet[j],add=adz, 
#                                        xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{1}, lty=4, ylab="LogNMix")
    }
    if(modcode[j] == 7){
      simvals$pdet[j] <- pweibull(10,simvals$shape_k[j],exp(-simvals$intcpt_d[j]))
      ys <- rep((pweibull(taus[-1],simvals$shape_k[j],exp(-simvals$intcpt_d[j])) - 
                   pweibull(taus[-length(taus)],simvals$shape_k[j],exp(-simvals$intcpt_d[j])))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*ys[1:2]
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="Weibull")
#       if(simvals$n_ints[j]==9) curve(dweibull(x,simvals$shape_k[j],exp(-simvals$intcpt_d[j]))/simvals$pdet[j], 
#                                      add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{2}, lty=1, ylab="Weibull")
    }
    if(modcode[j] == 8){
      simvals$pdet[j] <- 1-simvals$gamma[j]*(1-pweibull(10,simvals$shape_k[j],exp(-simvals$intcpt_d[j])))
      ys <- simvals$gamma[j]*rep((pweibull(taus[-1],simvals$shape_k[j],exp(-simvals$intcpt_d[j])) - 
                                    pweibull(taus[-length(taus)],simvals$shape_k[j],exp(-simvals$intcpt_d[j])))/simvals$pdet[j],each=2)
      ys[1:2] <- 0.5*(ys[1:2] + (1 - simvals$gamma[j])/simvals$pdet[j])
      lines(ends, ys, col=if(simvals$indx[j]==0){"orange"}else{6}, lty=5, ylab="WeibullMix")
#       if(simvals$n_ints[j]==9) curve(((x<=2)*(1-simvals$gamma[j]*(1-pweibull(2,simvals$shape_k[j],exp(-simvals$intcpt_d[j]))))/2 +
#                                        (x>2)*simvals$gamma[j]*dweibull(x,simvals$shape_k[j],exp(-simvals$intcpt_d[j])))/
#                                        simvals$pdet[j],
#                                        add=adz, xlim=c(0,10), col=if(simvals$indx[j]==0){"orange"}else{6}, lty=5, ylab="WeibullMix")
    }
  }
  par(mfrow=c(1,1))
  return(simvals)
}