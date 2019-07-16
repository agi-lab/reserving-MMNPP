###################################################
#### MMNPP calibration from simulated examples ####
###################################################

### Define parameters to be updated ###
library(expm)

# Set dimensions for the MMNPP
dim = 3

# For illustratory purposes, you can instead set the number of EM iterations. Otherwise, this controls the maximum number of iterations.
numiter = 100

# Otherwise, set tolerance for algorithm termination
# tol = 1e-3

Q <- array(data=0,dim = c(dim,dim,numiter+1))
Q[,,1] <- t(matrix(rep(0.1,times=dim^2),nrow=dim,ncol=dim))
for (i in 1:dim){
  Q[i,i,1] = -(sum(Q[i,,1]) - Q[i,i,1])
}

Lambda <- array(data=0,dim=c(dim,dim,numiter+1))
Lambda[,,1] <- diag(seq(from=5, by=5, length.out=dim),nrow=dim,ncol=dim)


### Define density functions ###
f_density<-function(y,Q,Lambda,gamma,gamma.next){
  return((expm((Q-Lambda*gamma)*y,method = c("Pade"))%*%Lambda)*gamma.next)
}

fbar_density<-function(y,Q,Lambda,gamma){
  return(expm((Q-Lambda*gamma)*y,method = c("Pade"))%*%Lambda)
}


### Define forward/backward recursion vectors ###
L.recur <-matrix(data=0,nrow=numclaims+1,ncol=dim)
R.recur <-matrix(data=0,nrow=dim,ncol=numclaims+1)
c.recur <-array(data=0,dim=c(1,numclaims))
L.recur[1,] = rep(x=1/dim,times=dim)
R.recur[,numclaims+1] = t(rep(x=1,times=dim))


### Initialise vector for loglikelihood values
loglikelihood <- c()

### Calculate interarival times ###
intarrtimes <- c(arrtimes[1],diff(arrtimes))

### Set final gamma equal to last gamma for recursion purposes ###
gamma <- c(gamma,gamma[length(gamma)])

### EM Algorithm ###
for (loopcount in 1:numiter) {

  for (i in 1:(numclaims)) {
    c.recur[i] = (L.recur[i,]%*%f_density(y= intarrtimes[i],Q = Q[,,loopcount],Lambda = Lambda[,,loopcount],gamma = gamma[i],gamma.next = gamma[i+1])%*%rep(x=1,times=dim))
    L.recur[i+1,] = (L.recur[i,]%*%f_density(y = intarrtimes[i],Q = Q[,,loopcount],Lambda = Lambda[,,loopcount],gamma = gamma[i],gamma.next = gamma[i+1])) / c.recur[i]
  }
  for (i in numclaims:1){
    R.recur[,i] = (f_density(y = intarrtimes[i],Q = Q[,,loopcount],Lambda = Lambda[,,loopcount],gamma = gamma[i],gamma.next = gamma[i+1]) %*% R.recur[,i+1]) / c.recur[i]
  }


  ### Calculate M-step parameters ###
  # m_{i,j}
  m.hat <- array(data=0,dim=c(dim,dim))

  # Calculate I_k
  C.k.mat <- array(data=0,dim=c(dim*2,dim*2,numclaims))
  I.k     <- array(data=0,dim=c(dim,dim,numclaims))
  for (k in 1:numclaims) {
    C.k.mat[,,k] = cbind(rbind(Q[,,loopcount]-Lambda[,,loopcount]*gamma[k],matrix(data=0,nrow=dim,ncol=dim)),rbind((Lambda[,,loopcount]*gamma[k+1])%*%R.recur[,k+1]%*%L.recur[k,],Q[,,loopcount]-Lambda[,,loopcount]*gamma[k]))
    I.k[,,k] <- (expm(C.k.mat[,,k]*intarrtimes[k],method=c("Pade")))[1:dim,(dim+1):(dim*2)]
  }

  I.k.sum <-array(data=0,dim=c(dim,dim))
  for (k in 1:numclaims) {
    I.k.sum <- I.k.sum + (I.k[,,k] / c.recur[k])
  }
  m.hat <- Q[,,loopcount] * t(I.k.sum)

  # n_{i}
  n.hat <- array(data=0,dim=c(dim,1))
  for (k in 1:numclaims) {
    n.hat <- n.hat + (L.recur[k+1,] * R.recur[,k+1])
  }

  # T_{i}
  t.hat <- array(data=0,dim=c(dim,1))
  t.hat <- diag(I.k.sum)

  # T*_{i}
  t.hat.star <- array(data=0,dim=c(dim,1))
  I.k.sum.star <-array(data=0,dim=c(dim,dim))
  for (k in 1:numclaims) {
    I.k.sum.star <-I.k.sum.star + (gamma[k] * t(I.k[,,k]) / c.recur[k])
  }
  t.hat.star <- diag(I.k.sum.star)

  #Sub in estimators into MLE functions
  Q[,,loopcount+1] <- m.hat / t.hat


  for (i in 1:dim){
    Q[i,i,loopcount+1]<-0
    Q[i,i,loopcount+1] <- -(sum(Q[i,,loopcount+1]))
  }

  Lambda[,,loopcount+1]    <- diag(as.vector(n.hat) / as.vector(t.hat.star))
  loglikelihood <- c(loglikelihood,sum(log(c.recur)))
  print(loopcount)
  loopcount <- loopcount+1
  
  # # Terminate if loglikelihood difference is below tolerance
  # if(abs(diff(tail(loglikelihood,n=2)))<tol){
  #   break
  # }
}

#Clear elements to reduce file size
rm("C.k.mat")
rm("I.k")