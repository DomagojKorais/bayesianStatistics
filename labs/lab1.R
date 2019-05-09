require(tidyverse)

set.seed(123)


#creo hyperparameters
h1 = c(1,1,1,2,.5)  
h2 = c(0.5,2,10,2,0)  



priordist <- function(theta, k ){
  return(dgamma(theta,a[k],b[k]))
}
#punto 1: 
curve(dgamma(x,1,0.5), from=0, to=10,ylim=c(0,3) , xlab="x", ylab="y")
curve(dgamma(x,1,2), from=0, to=10, , xlab="x", ylab="y",add = TRUE,col = "red")
curve(dgamma(x,1,10), from=0, to=10, , xlab="x", ylab="y",,add = TRUE,col = "blue")
curve(dgamma(x,2,2), from=0, to=10, , xlab="x", ylab="y",,add = TRUE,col = "green")
curve(x**(-1/2), from=0, to=10, , xlab="x", ylab="y",,add = TRUE,col = "yellow")

##quindi la meno informativa è la Jeffrey (code più lunghe), gamma 1,10 più informativa 
##perchè è la più piccata.
##ggplot2 equivalent

#punto 2:
a=rpois(20,0.5)
b=rpois(100,0.5)
c=rpois(1000,0.5)

sum(a)
#posterior per la prior Poisson
posterior <- function(x,arrayPoisson,alpha,beta){
  alpha=alpha + sum(arrayPoisson)
  beta=beta+length(arrayPoisson)
  return (dgamma(x,alpha,beta))
}

#posterior per la Jeffrey
posteriorJeffrey <- function(x,arrayPoisson,alpha,beta){
  alpha=0.5 + sum(arrayPoisson)
  beta=length(arrayPoisson)
  return (dgamma(x,alpha,beta))
}

#punto 2 plot:

  
for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  curve(posterior(x,a,h1[i],h2[i]),col=i,add=flag, from=0, to=1.5,ylim=c(0,6) , xlab="x", ylab="y")
}

for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  curve(posterior(x,b,h1[i],h2[i]),col=i,add=flag, from=0, to=1.5,ylim=c(0,10) , xlab="x", ylab="y")
}
for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  curve(posterior(x,c,h1[i],h2[i]),col=i,add=flag, from=0, to=1.5,ylim=c(0,20) , xlab="x", ylab="y")
}

#punto 3

#mean
posteriorExpetations = function(arrayPoisson,alpha,beta){
  alpha=alpha + sum(arrayPoisson)
  beta=beta+length(arrayPoisson)
  return (alpha/beta)
}
#Maxima a posteriora (MAP)
computeMap = function(arrayPoisson,alpha,beta){
  alpha=alpha + sum(arrayPoisson)
  beta=beta+length(arrayPoisson)
  return ((alpha-1)/beta)
}

for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  print(paste("mean:", posteriorExpetations(a,h1[i],h2[i])))
  print(paste("median:", posteriorExpetations(a,h1[i],h2[i])))
}

#Punto 5 


  