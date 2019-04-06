
set.seed(1,2,3)


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

#punto 2:
a=rpois(20,0.5)
b=rpois(100,0.5)
c=rpois(1000,0.5)
sum(a)

posterior <- function(x,arrayPoisson,beta){
  alpha=sum(arrayPoisson)
  beta=beta+length(arrayPoisson)
  return (dgamma(x,alpha,beta))
}

#punto 2 plot:

  
for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  curve(posterior(x,a,h2[i]),col=i,add=flag, from=0, to=1.5,ylim=c(0,6) , xlab="x", ylab="y")
}

for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  curve(posterior(x,b,h2[i]),col=i,add=flag, from=0, to=1.5,ylim=c(0,10) , xlab="x", ylab="y")
}
for (i in 1:5)
{
  if (i ==1){flag=FALSE}
  else {flag=TRUE}
  curve(posterior(x,c,h2[i]),col=i,add=flag, from=0, to=1.5,ylim=c(0,20) , xlab="x", ylab="y")
}

