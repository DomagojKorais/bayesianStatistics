#Lab 3 
library(gclus)
data(bank)
bank <- as.matrix(bank)

y <- bank[,1]
X <- cbind(1,bank[,2:5])
q <- dim(X)[2] 
n <- dim(X)[1]



posterior = function(b,X,y){
  pi = pnorm(X %*% matrix(b,ncol=1))
  print(pi)
  return (prod((pi^y)*(1-pi)^y))
}

library(mvtnorm) 
mh = function(Nsim, tau, y, X){
  randomValues = runif(Nsim)
  beta = matrix(0, nrow=Nsim,ncol=ncol(X))
  sigma.asymp = summary(glm(y~X,family = binomial(link = "probit")))$cov.unscaled
  beta[1,] = summary(glm(y~X,family = binomial(link = "probit")))$coefficients[,1]
 
  for (i in 2:Nsim) {
    beta[i,] = rmvnorm(1,beta[i-1,], tau^2*sigma.asymp)
    ro = min(1,posterior(beta[i,],X[i,],y)/posterior(beta[i-1,],X[i-1,],y))
    if (randomValues[i]< ro) 
      {beta[i, ] = beta[i, ]}
    else
      {beta[i, ] = beta[i-1, ]}
  }
  return(beta)
}  


mhv=mh(1000,0.1,y,X)

##traceplot
par(mfrow=c(2,2))
plot(mhv, type="l", ylab="Trace", xlab="Iteration")
##histogram
hist(mhv, breaks=100, border="gray40",freq=F, main="")
lines(seq(-5,5, by=.01), dt(seq(-5,5, by=.01),4),type="l", col="gray30")
abline(v=mean(mhv), col="firebrick3", lty=2)
##cumulative mean
a <- cumsum(mhv)/1:Nsim
plot(a, type="l",ylab="Cumulative mean plot", xlab="Iteration")
abline(h=mean(mhv), col="firebrick3", lty=2)
##autocorrelation
acf(mhv, main="",ylab="Autocorrelation")


#Ghibbs
library(truncnorm)

ghibbs = function(Nsim, y, X){
  
  beta = matrix(0, nrow=Nsim,ncol=ncol(X))
  sigma.asymp = summary(glm(y~X,family = binomial(link = "probit")))$cov.unscaled
  beta[1,] = summary(glm(y~X,family = binomial(link = "probit")))$coefficients[,1]
  
 for (i in 2:Nsim) {
      z=0
    if (y[i] == 1) {
      z = rtruncnorm(1,a=0,mean=X[i,]*beta[i,])
    }
    else {
      z = rtruncnorm(1,b=0,mean=X[i,]*beta[i,])
    }
  
  beta[i,] = rmvnorm(1,solve(t(X[i,])%*%X[i,]) %*%t(X[i,])%*%z, solve(t(X[i,])%*%X[i,]) )
  
  }
  return(beta)
}  



ghibbs(1000,y,X)