###Runs test against a Poisson process###
library(randtests)
#To check that means and medians are close
median.1order <-median(rawdata.df$FREQ/rawdata.df$EXP)
mean.1order <-mean(rawdata.df$FREQ/rawdata.df$EXP)
#Runs test
runs.test <- runs.test(x=(rawdata.df$FREQ/rawdata.df$EXP),alternative = "l")
runs.test$parameter
runs.test$p.value



###Residual Test###
library(randtests)
library(stats)

#Determine daily regimes based on mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

regimes <- max.col(L.recur,'first')
regimechanges.num <- sum(diff(regimes)!=0)
cum.claimcounts <- cumsum(rawdata.df$FREQ)
daily.regime <- c(Mode(regimes[1:cum.claimcounts[1]]))

for (i in 2:length(cum.claimcounts)){
  daily.regime<-c(daily.regime,Mode(regimes[cum.claimcounts[i-1]:cum.claimcounts[i]]))
}

#Calculate predicted frequency
freq.pred <- array(data=0,dim=c(length(daily.regime)))
for (i in 1:length(daily.regime)){
  freq.pred[i] <- Lambda[daily.regime[i],daily.regime[i],numiter+1]*rawdata.df$EXP[i]
}

#Calculate residuals
freq.res <- array(data=0,dim=c(length(daily.regime)))
for (i in 1:length(daily.regime)){
  freq.res[i] <- rawdata.df$FREQ[i] - freq.pred[i]
}

#Calculate standardised residuals
freq.stdres <- array(data=0,dim=c(length(daily.regime)))
for (i in 1:length(daily.regime)){
  freq.stdres[i] <- (rawdata.df$FREQ[i] - freq.pred[i])/sqrt(freq.pred[i])
}

#Runs Test
runstest <- runs.test(x=as.vector(freq.res),alternative = "l",threshold=0)
runstest$p.value

runstest <- runs.test(x=as.vector(freq.stdres),alternative = "l",threshold=0)
runstest$p.value

#Ljung Box test for autocorrelation
autotest <- Box.test(x=as.vector(freq.res),lag = 120, type="Ljung-Box")
pacf(as.vector(freq.res),lag=120)
autotest$p.value

autotest <- Box.test(x=as.vector(freq.stdres),lag = 120, type="Ljung-Box")
pacf(as.vector(freq.stdres),lag=120)
autotest$p.value

#Bartlett B test for white noise
library(hwwntest)
bartlettB.test(as.vector(freq.res))
bartlettB.test(as.vector(freq.stdres))
hwwn.test(as.vector(freq.res[1:2048]))
hwwn.test(as.vector(freq.stdres[1:2048]))
