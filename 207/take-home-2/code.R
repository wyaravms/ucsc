
rm(list=ls(all=TRUE))
set.seed(7)

# data 
data=read.table("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\spring-2018\\ams-207\\take-home-2\\data\\volcano.csv", head=TRUE,sep=",")

data=data[-nrow(data),]

plot(density(data[,5]))

plot(density(log(data[,5])), main="")

volcano = log(data[,5])
n = length(volcano)
mean(volcano)

hist(volcano, xlab="Log of interevent time",main="Kilauea volcano data",breaks = 20,freq=FALSE,cex.lab=1.4,cex.main = 1.9,cex.axis=1.7)
lines(density(volcano),col="red")
#number of degrees of freedom
v=5

# gibbs sampling

# model 1
## functions to sample each parameter from its full conditionals

### update lambda
fn.update.lambda = function(data,v, mu, sigma)
{
  lambda = rgamma(1, ((v+1)/2), (v/2)+(((data - mu)^2)/(2*sigma)))
  return(lambda)
}

## udpate mu.i
fn.update.mu.i = function(data,lambda,mu,sigma,tau)
{
  va.mu.i = (tau*sigma)/((lambda*tau) +sigma)
  mean.mu.i = ((data*lambda*tau) + (mu*sigma))/((lambda*tau) + sigma) 
  
  mu.i = rnorm(1, mean.mu.i, sqrt(va.mu.i))
  return(mu.i)
}

## udpate mu
fn.update.mu = function(mu.i,n,tau,m,s)
{
  va.mu = (tau*s)/(tau + (n*s))
  mean.mu = ((s*sum(mu.i)) + (m*tau))/(tau+(n*s)) 
  
  mu = rnorm(1, mean.mu, sqrt(va.mu))
  return(mu)
}

## update sigma
fn.update.sigma = function(data,n,mu.i,lambda,a,b)
{
  sigma = 1/rgamma(1, ((n/2) + a), (b + sum(((data-mu.i)^2)*lambda)/2))  
  return(sigma)
}

## update tau
fn.update.tau = function(n,mu.i,mu,c,d)
{
  tau = 1/rgamma(1, ((n/2) + c), (sum((mu.i-mu)^2)/2 + d))
  return(tau)
}

## initial values
par.sam = NULL
par.sam$lambda = rep(0, n)
par.sam$mu.i = rep(0, n)
par.sam$mu.i[1] = 5
par.sam$mu = 5
par.sam$sigma = 0.5
par.sam$tau = 0.5
par.sam$m = 5
par.sam$s = 8
par.sam$a = 3
par.sam$b = 3
par.sam$c = 3
par.sam$d = 3

## variables for the MCMC
ns = 41000

## save simulated para
MCMC.sam.M1 = NULL
MCMC.sam.M1$lambda = array(NA, dim=c(n, ns))
MCMC.sam.M1$mu.i = array(NA, dim=c(n, ns))
MCMC.sam.M1$mu = rep(NA, ns)
MCMC.sam.M1$sigma = rep(NA, ns)
MCMC.sam.M1$tau = rep(NA, ns)

## MCMC modeling

for(i.iter in 1:ns)
{
  if((i.iter%%1000)==0)
  {
    print(paste("i.iter=", i.iter))
    print(date())
  }
  
  for (j in 1:n){
    
    ## udpate lambda
    par.sam$lambda[j] = fn.update.lambda(volcano[j],v,par.sam$mu.i[j],par.sam$sigma)
    
    ## udpate mu
    par.sam$mu.i[j] = fn.update.mu.i(volcano[j],par.sam$lambda[j], par.sam$mu, par.sam$sigma,par.sam$tau)
    
  }
  
  ## update mu
  par.sam$mu = fn.update.mu(par.sam$mu.i,n,par.sam$tau,par.sam$m,par.sam$s)
  
  ## udpate tau
  par.sam$tau = fn.update.tau(n,par.sam$mu.i,par.sam$mu,par.sam$c,par.sam$d)
  
  ## udpate sigma
  par.sam$sigma = fn.update.sigma(volcano,n,par.sam$mu.i,par.sam$lambda,par.sam$a,par.sam$b)
  
  ## save cur.sam
  MCMC.sam.M1$lambda[,i.iter] = par.sam$lambda
  MCMC.sam.M1$mu.i[,i.iter] = par.sam$mu.i
  MCMC.sam.M1$mu[i.iter] = par.sam$mu
  MCMC.sam.M1$sigma[i.iter] = par.sam$sigma
  MCMC.sam.M1$tau[i.iter] = par.sam$tau
}

n.bur = 1000
thin = 4

lambda.m1 = MCMC.sam.M1$lambda[,seq(n.bur+1, ns,by=thin)]
mu.i.m1 = MCMC.sam.M1$mu.i[,seq(n.bur+1, ns,by=thin)]
mu.m1 = MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)]
sigma.m1 = MCMC.sam.M1$sigma[seq(n.bur+1, ns,by=thin)]
tau.m1 = MCMC.sam.M1$tau[seq(n.bur+1, ns,by=thin)]

# plots of the chain
par(mfrow=c(2,3))
plot.ts(mu.m1,cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(mu)))
plot.ts(sigma.m1,cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma)))
plot.ts(tau.m1,cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(tau)))
dev.off()

# some acf plots
par(mfrow=c(1,3))
acf(mu.m1,cex.lab=1.3,cex.main = 3,cex.axis=1.3, col=1, xlab="",main="")
acf(sigma.m1,cex.lab=1.3,cex.main = 3,cex.axis=1.3, col=1, xlab="",main="")
acf(tau.m1,cex.lab=1.3,cex.main = 3,cex.axis=1.3, col=1, xlab="",main="")
dev.off()

# densities
par(mfrow=c(3,1))
plot(density(mu.m1),cex.lab=1.5,cex.main = 3,cex.axis=1.7, col=1, xlab="", main=expression(paste(mu)))
abline(v=mean(mu.m1), lty=2, col=2)
abline(v=quantile(mu.m1, probs=c(0.025, 0.975)), lty=2, col=2)
plot(density(sigma.m1),cex.lab=1.5,cex.main = 3,cex.axis=1.7, col=1, xlab="", main=expression(paste(sigma)))
abline(v=quantile(sigma.m1, probs=c(0.025, 0.975)), lty=2, col=2)
abline(v=mean(sigma.m1), lty=2, col=2)
plot(density(tau.m1),cex.lab=1.5,cex.main = 3,cex.axis=1.7, col=1, xlab="",main=expression(paste(tau)))
abline(v=quantile(tau.m1, probs=c(0.025, 0.975)), lty=2, col=2)
abline(v=mean(tau.m1), lty=2, col=2)
dev.off()


par(mfrow=c(2,2))
# densities for the lambda's
plot(density(lambda.m1[1,]),ylim=c(0,1.7),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col="grey", xlab="",main=expression(paste(lambda[i])))
for(i in 2:length(volcano)){
  lines(density(lambda.m1[i,]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col="grey", xlab="",main=expression(paste(lambda[i])))
}
lines(density(lambda.m1[9,]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=2, xlab="",main=expression(paste(lambda[i])))
lines(density(lambda.m1[15,]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=4, xlab="",main=expression(paste(lambda[i])))
legend("topright",lwd=1, cex=1.2, box.lty = 0, legend=c(expression(paste(lambda[9])),expression(paste(lambda[15]))),bty='n',col=c(2,4))
# 9, 15

# densities for the mu's
plot(density(mu.i.m1[1,]),ylim=c(0,0.8),xlim=c(1,10),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col="grey", xlab="",main=expression(paste(mu[i])))
for(i in 2:length(volcano)){
  lines(density(mu.i.m1[i,]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col="grey", xlab="",main=expression(paste(lambda[i])))
}
lines(density(mu.i.m1[9,]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=2, xlab="",main=expression(paste(lambda[i])))
lines(density(mu.i.m1[15,]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=4, xlab="",main=expression(paste(lambda[i])))
legend("topright",lwd=1, cex=1.2, box.lty = 0, legend=c(expression(paste(mu[9])),expression(paste(mu[15]))),bty='n',col=c(2,4))
# 9, 15


apply(lambda.m1, 1, mean)
apply(mu.i.m1, 1, mean)
mean(mu.m1)
mean(sqrt(sigma.m1))
mean(sqrt(tau.m1))


rm(list=setdiff(ls(), "MCMC.sam.M1"))
set.seed(7)
require(invgamma)

# data 
data=read.table("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\spring-2018\\ams-207\\take-home-2\\volcano.csv", head=TRUE,sep=",")

data=data[-nrow(data),]

plot(density(data[,5]))

plot(density(log(data[,5])), main="")

volcano = log(data[,5])
n = length(volcano)
mean(volcano)

hist(volcano)

#number of degrees of freedom
v = 5

# gibbs sampling
# model 2
## functions to sample each parameter from its full conditionals

## udpate mu
fn.update.mu = function(data,n,sigma,m,s)
{
  va.mu = (sigma*s)/(sigma + (n*s))
  mean.mu = (s*sum(data) + m*sigma)/(sigma+(n*s)) 
  
  mu = rnorm(1, mean.mu, sqrt(va.mu))
  return(mu)
}

## update sigma
fn.update.sigma = function(data,n,mu,a,b)
{
  sigma = 1/rgamma(1, ((n/2) + a), ((sum((data-mu)^2)/2) + b))  
  return(sigma)
}

## initial values
par.sam = NULL
par.sam$mu = 5
par.sam$sigma = 1.5
par.sam$m = 6
par.sam$s = 8
par.sam$a = 3
par.sam$b = 3

## variables for the MCMC
ns = 41000

## save simulated para
MCMC.sam.M2 = NULL
MCMC.sam.M2$mu = rep(NA, ns)
MCMC.sam.M2$sigma = rep(NA, ns)

## MCMC modeling

for(i.iter in 1:ns)
{
  if((i.iter%%1000)==0)
  {
    print(paste("i.iter=", i.iter))
    print(date())
  }

  ## update mu
  par.sam$mu = fn.update.mu(volcano,n,par.sam$sigma,par.sam$m,par.sam$s)
  
  ## udpate sigma
  par.sam$sigma = fn.update.sigma(volcano,n,par.sam$mu,par.sam$a,par.sam$b)
  
  ## save cur.sam
  MCMC.sam.M2$mu[i.iter] = par.sam$mu
  MCMC.sam.M2$sigma[i.iter] = par.sam$sigma
}

n.bur = 1000
thin = 4

mu.m2 = MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)]
sigma.m2 = MCMC.sam.M2$sigma[seq(n.bur+1, ns,by=thin)]

# plots of the chain
par(mfrow=c(2,2))
plot.ts(mu.m2,cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(mu)))
plot.ts(sigma.m2,cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma)))
dev.off()

# some acf plots
par(mfrow=c(1,2))
acf(mu.m2,cex.lab=1.3,cex.main = 3,cex.axis=1.3, col=1, ylab="",main="")
acf(sigma.m2,cex.lab=1.3,cex.main = 3,cex.axis=1.3, col=1, ylab="",main="")
dev.off()

# densities
par(mfrow=c(2,1))
plot(density(mu.m2),cex.lab=1.5,cex.main = 2,cex.axis=1.6, col=1, xlab="", main=expression(paste(mu)))
abline(v=quantile(mu.m2, probs=c(0.025, 0.975)), lty=2, col=2)
abline(v=mean(mu.m2), lty=2, col=2)
plot(density(sigma.m2),cex.lab=1.5,cex.main = 2,cex.axis=1.6, col=1, xlab="", main=expression(paste(sigma)))
abline(v=quantile(sigma.m2, probs=c(0.025, 0.975)), lty=2, col=2)
abline(v=mean(sigma.m2), lty=2, col=2)
dev.off()

mean(mu.m2)
mean(sqrt(sigma.m2))

