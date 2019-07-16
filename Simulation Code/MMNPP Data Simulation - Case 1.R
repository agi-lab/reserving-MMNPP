############################################
### Simulated data using a order 3 MMNPP ###
############################################

#required inputs
set.seed(1000)
dim = 3
trueLambda = c(5,10,20)
trueQ = array(data=0,dim=c(3,3))
trueQ[1,2]=0.2
trueQ[1,3]=0.3
trueQ[2,1]=0.4
trueQ[2,3]=0.2
trueQ[3,1]=0.5
trueQ[3,2]=0.3
for (i in 1:dim){
  trueQ[i,i] = - (sum(trueQ[i,]) - trueQ[i,i])
}

modQ<-trueQ
for (i in 1:dim){
  modQ[i,i] = 0
}

sumQ1 <- sum(modQ[1,])
sumQ2 <- sum(modQ[2,])
sumQ3 <- sum(modQ[3,])

time=0

#Exposure-related variables
timelength<-50
exposure.values <-c(1,1.5,2,2.5,1.5,3,2.5,1,0.5,2)
exposure.times <-seq(1:10)*timelength


### No further inputs required ###
#simulate starting state
startstate <- ceiling(3*runif(1))
regimes <- c(startstate)
reg.time<- c(0)
reg.count=1
#simulate regimes
while(time < 10*timelength){
  jump.time <- rexp(n=1,sum(modQ[regimes[reg.count],]))
  #determine which state to jump to
  P <-runif(1)
  if(regimes[reg.count]==1){
    if(P<=trueQ[1,2]/sumQ1){
      next.regime = 2
    }
    else{
      next.regime = 3
    }
  }
  if(regimes[reg.count]==2){
    if(P<=trueQ[2,1]/sumQ2){
      next.regime = 1
    }
    else{
      next.regime = 3
    }
  }
  if(regimes[reg.count]==3){
    if(P<=trueQ[3,1]/sumQ3){
      next.regime = 1
    }
    else{
      next.regime = 2
    }
  }
  regimes<-c(regimes,next.regime)
  time = time + jump.time
  reg.time <- c(reg.time,time)
  reg.count = reg.count + 1
}


#Find indexes of elements before exposure changes
oldregindex <-vector(mode = "integer",length = 10)
for (i in 1:10){
oldregindex[i] <- max(which(reg.time<i*timelength))
reg.time<-append(reg.time,values=i*timelength,after=oldregindex[i])
regimes<-append(regimes,values=regimes[oldregindex[i]],after=oldregindex[i])
}

totallength<-length(reg.time)

reg.exposures <-vector(mode = "integer",length = length(reg.time))
#Add exposure numbers
for (counter in 1:(length(reg.time)-2)){
  reg.exposures[counter] <- exposure.values[floor(reg.time[counter]/timelength)+1]
}
reg.exposures[totallength-1]<-tail(exposure.values,1)


reg.time<-reg.time[1:(totallength-1)]
reg.exposures<-reg.exposures[1:(totallength-1)]
regimes<-regimes[1:(totallength-1)]

regime.df<-data.frame(reg.time,regimes,reg.exposures)
diff.times<-diff(reg.time)
numclaimssim<-vector(mode="integer",length=length(diff.times))
#simulate claims
for(i in 1:length(diff.times)){
  numclaimssim[i]<-rpois(n = 1,lambda = reg.exposures[i]*trueLambda[regimes[i]]*diff.times[i])
}

numclaims<-sum(numclaimssim)

#simulate claim interarrival times
arrtimes<-c()
gamma<-c()
for(i in 1:length(diff.times)){
  arrtimes<-c(arrtimes,sort(runif(n = numclaimssim[i],min=reg.time[i],max=reg.time[i+1])))
  gamma<-c(gamma,rep(x=reg.exposures[i],times=numclaimssim[i]))
}
gamma<-c(gamma[1],gamma)

