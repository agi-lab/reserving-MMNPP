#Deterministic periodic regime changes for order 3 MMNPP
# Set Regime Order
reg.order <- c(1,2,3,1,3,2,1,2,2,3)
exposure <- c(1,2,2,1,3,8,5,4,6,2)
base.freq <- c(50,100,200)
length.periods <- 100

############ No other inputs required ##################
set.seed(1)
period.freq <- array(data = 0, dim = length(reg.order))
claimtimes <- c()
gamma <- c()
for (i in 1:length(reg.order)){
  period.freq[i] <- rpois(1,exposure[i]*base.freq[reg.order[i]])
  claimtimes<-c(claimtimes,sort(runif(n = period.freq[i],min = (i-1)*length.periods,max = i*length.periods)))
  gamma <- c(gamma,rep(x=exposure[i],times=period.freq[i]))  
}

numclaims <- sum(period.freq)
