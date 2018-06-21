
rm(list=ls(all=TRUE))
library(dplyr)

setwd("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\winter-2018\\ams-206b\\final")

load("Machine1.RData")

y =  Machine
y = (as.data.frame(y))

# some summaries for sample mean and sample variance
aggregate(y, by=list(y$machine), FUN=mean, na.rm=TRUE)
aggregate(y, by=list(y$machine), FUN=var, na.rm=TRUE)

# boxplot graph
boxplot(measurements ~ machine, main="Quality control measurements", xlab="Machines")

n = 6
mm = nrow(y)

# number of each measurement
m1=9; m2=8; m3=7; m4=6;
m5=6; m6=8;

m = c(m1,m2,m3,m4,m5,m6)

# Question 2
# (c)

#Gibbs Sampling

## functions to sample each parameter from its full conditionals

### update theta.i
fn.update.theta = function(y, m, mu, sig2, tau2)
{
  theta=NULL
  for(j in 1:length(m))
  {
    #m.al = c(0,m)
    theta[j] = rnorm(1, (sum(y$measurements[y$machine==j])/sig2 + mu/tau2)/((m[j]/sig2)+(1/tau2)), sqrt(1/((m[j]/sig2)+(1/tau2))))
      #(sum(y[(sum(m.al[1:j],1)):(sum(m.al[1:(j+1)])),2])/sig2 + mu/tau2)
  }
  return(theta)
}

## update sig2
fn.update.sig2 = function(mm, v0, s0, y.theta.sq)
{
  sig2 = 1.0/rgamma(1, (mm/2)+(v0/2), (y.theta.sq)/2 + (s0^2/2))  ## mean a/b
  return(sig2)
  
}

## udpate mu
fn.update.mu = function(theta.s, tau2, mu0, omega2, n)
{
  va = 1/(n/tau2 + 1/omega2)
  mme = va*((theta.s/tau2) + (mu0/omega2))
  
  mu = rnorm(1, mme, sqrt(va))
  return(mu)
  
}

## update tau2
fn.update.tau2 = function(n, at, y.theta.mu.sq, bt)
{
  tau2 = 1.0/rgamma(1, (n/2)+at, y.theta.mu.sq/2 + bt)
  return(tau2)
  
}

# hyperparamenters generating higher variance on the distribution of the priors,
# being so a vague prior

## hyperparameter
hyper = NULL
hyper$v0 = 6
hyper$s0 = 3

hyper$mu0 = 90
hyper$omega2 = 100

hyper$at = 3
hyper$bt = 3

## initial values
par.sam = NULL
par.sam$theta = rep(mean(y$measurements), n)
par.sam$sig2 = 10^3
par.sam$mu = 90
par.sam$tau = 10^3

## variables for the MCMC
ns = 40000

## save simulated para
MCMC.sam.M1 = NULL
MCMC.sam.M1$theta = array(NA, dim=c(n, ns))
MCMC.sam.M1$sig2 = rep(NA, ns)
MCMC.sam.M1$mu = rep(NA, ns)
MCMC.sam.M1$tau = rep(NA, ns)

## MCMC modeling

for(i.iter in 1:ns)
{
  if((i.iter%%1000)==0)
  {
    print(paste("i.iter=", i.iter))
    print(date())
  }
  
  ## udpate theta
  par.sam$theta = fn.update.theta(y, m, par.sam$mu, par.sam$sig2, par.sam$tau)
  
  #theta.vector = function(par.theta,m){
    for(i in 1:length(m)){
      if(i==1){
        par.sam$theta.list=rep(par.sam$theta[i],m[i])
        } else{
      p.the = rep(par.sam$theta[i],m[i])
      par.sam$theta.list=c(par.sam$theta.list,p.the)}
    }
    #return()}

  y$theta.list = par.sam$theta.list
  #sum.dif = sapply(y$machine, function(machine) { y$measurements[y$machine==machine]-par.sam$theta } )
  ## update sig2
  par.sam$sig2 = fn.update.sig2(mm, hyper$v0, hyper$s0, (sum((y$measurements-y$theta.list)^2)))
  
  ## udpate mu
  par.sam$mu = fn.update.mu(sum(par.sam$theta), par.sam$tau, hyper$mu0, hyper$omega2, n)
  
  ## update tau2
  par.sam$tau = fn.update.tau2(n, hyper$at,  sum((par.sam$theta - par.sam$mu)^2), hyper$bt)
  
  ## save cur.sam
  MCMC.sam.M1$theta[,i.iter] = par.sam$theta
  MCMC.sam.M1$sig2[i.iter] = par.sam$sig2
  MCMC.sam.M1$mu[i.iter] = par.sam$mu
  MCMC.sam.M1$tau[i.iter] = par.sam$tau
    
}

n.bur = 1000
thin = 3

#posterior mean for the parameters
apply(MCMC.sam.M1$theta[,seq(n.bur+1, ns,by=thin)],1,mean)

mean(MCMC.sam.M1$sig2[seq(n.bur+1, ns,by=thin)])
mean(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)])
mean(MCMC.sam.M1$tau[seq(n.bur+1, ns,by=thin)])

# credibility intervals for the parameters
apply(MCMC.sam.M1$theta[,seq(n.bur+1, ns,by=thin)], 1, quantile,  probs=c(0.025, 0.975))

quantile(MCMC.sam.M1$sig2[seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC.sam.M1$tau[seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975))

# Question 3

rm(list=setdiff(ls(), "MCMC.sam.M1"))
library(dplyr)

setwd("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\winter-2018\\ams-206b\\final")

load("Machine1.RData")

y= Machine
y=(as.data.frame(y))

n = 6
mm = nrow(y)

# number of each measurement
m1=9; m2=8; m3=7; m4=6;
m5=6; m6=8;

m = c(m1,m2,m3,m4,m5,m6)


# (c)

#Gibbs Sampling

## functions to sample each parameter from its full conditionals

### update theta.i
fn.update.theta = function(y, m, mu, sig2, tau2)
{
  theta=NULL
  for(j in 1:length(m))
  {
    theta[j] = rnorm(1, (sum(y$measurements[y$machine==j])/sig2[j] + mu/tau2)/((m[j]/sig2[j])+(1/tau2)), sqrt(1/((m[j]/sig2[j])+(1/tau2))))
  }
  return(theta)
}

## update sig2
fn.update.sig2 = function(m, v0, s0, y.theta.sq.s)
{
  sigma=NULL
  for(j in 1:length(m))
  {
    sigma[j] = 1.0/rgamma(1, (m[j]/2)+(v0/2), (y.theta.sq.s[j])/2 + (s0^2/2)) 
  }
  return(sigma)
  
}

## update s02
fn.update.s02 = function(n, as, sig.sum, bs)
{
  s02 = 1.0/rgamma(1, as, (1/(sig.sum*2)) + bs)
  return(s02)
  
}

## udpate mu
fn.update.mu = function(theta.s, tau2, mu0, omega2, n)
{
  va = 1/(n/tau2 + 1/omega2)
  mme = va*((theta.s/tau2) + (mu0/omega2))
  
  mu = rnorm(1, mme, sqrt(va))
  return(mu)
  
}

## update tau2
fn.update.tau2 = function(n, at, y.theta.mu.sq, bt)
{
  tau2 = 1.0/rgamma(1, (n/2)+at, y.theta.mu.sq/2 + bt)
  return(tau2)
  
}

# hyperparamenters generating higher variance on the distribution of the priors,
# being so a vague prior

## hyperparameter
hyper = NULL
hyper$v0 = 6
hyper$s0 = 3

hyper$mu0 = 90
hyper$omega2 = 100

hyper$at = 3
hyper$bt = 3

hyper$as = 3
hyper$bs = 3

## initial values
par.sam = NULL
par.sam$theta = rep(mean(y$measurements), n)
par.sam$sig2 = rep(var(y$measurements), n)
par.sam$s02 = 10^3
par.sam$mu = 90
par.sam$tau = 10^3

## variables for the MCMC
ns = 40000

## save simulated para
MCMC.sam.M2 = NULL
MCMC.sam.M2$theta = array(NA, dim=c(n, ns))
MCMC.sam.M2$sig2 = array(NA, dim=c(n, ns))
MCMC.sam.M2$s02 = rep(NA, ns)
MCMC.sam.M2$mu = rep(NA, ns)
MCMC.sam.M2$tau = rep(NA, ns)

## MCMC modeling

for(i.iter in 1:ns)
{
  if((i.iter%%1000)==0)
  {
    print(paste("i.iter=", i.iter))
    print(date())
  }
  
  ## udpate theta
  par.sam$theta = fn.update.theta(y, m, par.sam$mu, par.sam$sig2, par.sam$tau)
  
  ## creating a vector of theta to evaluate sigma^2
  for(ii in 1:length(m)){
    if(ii==1){
      par.sam$theta.list=rep(par.sam$theta[ii],m[ii])
    } else{
      p.the = rep(par.sam$theta[ii],m[ii])
      par.sam$theta.list=c(par.sam$theta.list,p.the)}
  }
  
  y$theta.list = par.sam$theta.list
  y.theta.dif=NULL
  for(jj in 1:length(m)){
    y.theta.dif[jj] = sum((y$measurements[y$machine==jj] - y$theta.list[y$machine==jj])^2)
  }
  
  par.sam$sig2 = fn.update.sig2(m, hyper$v0, hyper$s0, y.theta.dif)
  
  par.sam$s02 = fn.update.s02(n, hyper$as, sum(par.sam$sig2), hyper$bs)
  
  ## udpate mu
  par.sam$mu = fn.update.mu(sum(par.sam$theta), par.sam$tau, hyper$mu0, hyper$omega2, n)
  
  ## update tau2
  par.sam$tau = fn.update.tau2(n, hyper$at, sum((par.sam$theta - par.sam$mu)^2), hyper$bt)
  
  ## save cur.sam
  MCMC.sam.M2$theta[,i.iter] = par.sam$theta
  MCMC.sam.M2$sig2[,i.iter] = par.sam$sig2
  MCMC.sam.M2$mu[i.iter] = par.sam$mu
  MCMC.sam.M2$tau[i.iter] = par.sam$tau
  MCMC.sam.M2$s02[i.iter] = par.sam$s02
  
}

n.bur = 1000
thin = 3

#posterior mean for the parameters
apply(MCMC.sam.M2$sig2[,seq(n.bur+1, ns,by=thin)],1,mean)
apply(MCMC.sam.M2$theta[,seq(n.bur+1, ns,by=thin)],1,mean)

mean(MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)])
mean(MCMC.sam.M2$tau[seq(n.bur+1, ns,by=thin)])
mean(MCMC.sam.M2$s02[seq(n.bur+1, ns,by=thin)])

# credibility intervals for the parameters
apply(MCMC.sam.M2$theta[,seq(n.bur+1, ns,by=thin)], 1, quantile,  probs=c(0.025, 0.975))
apply(MCMC.sam.M2$sig2[,seq(n.bur+1, ns,by=thin)], 1, quantile,  probs=c(0.025, 0.975))

quantile(MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC.sam.M2$tau[seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC.sam.M2$s02[seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975))

# Question 4

# Computing deviance information criteria
ndic = length(MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)])
l = vector("list", 6)
for(i in 1:6){
  theta.m1 = matrix(rep(MCMC.sam.M1$theta[i,seq(n.bur+1, ns,by=thin)],m[i]), nrow = ndic, ncol = m[i])
  l[[i]] = theta.m1 
}
theta.m1 = do.call("cbind", l)
sigma.m1 = (MCMC.sam.M1$sig2[seq(n.bur+1, ns,by=thin)])

dev.m1 = matrix(0, ndic, nrow(y))
for (i in 1:nrow(y)){
  dev.m1[,i] = dnorm(y$measurements[i], theta.m1[,i], sqrt(sigma.m1), log = TRUE)
}

l = vector("list", 6)
for(i in 1:6){
  theta.m2 = matrix(rep(MCMC.sam.M2$theta[i,seq(n.bur+1, ns,by=thin)],m[i]), nrow = ndic, ncol = m[i])
  l[[i]] = theta.m2
}
theta.m2 = do.call("cbind", l)

l = vector("list", 6)
for(i in 1:6){
  sigma.m2 = matrix(rep(MCMC.sam.M2$sig2[i,seq(n.bur+1, ns,by=thin)],m[i]), nrow = ndic, ncol = m[i])
  l[[i]] = sigma.m2 
}
sigma.m2 = do.call("cbind", l)

dev.m2 = matrix(0, ndic, nrow(y))
for (i in 1:nrow(y)){
  dev.m2[,i] = dnorm(y$measurements[i], theta.m2[,i], sqrt(sigma.m2[,i]), log = TRUE)
}

dev.m1 = -2*apply(dev.m1, 1, sum)
m1.DIC = mean(dev.m1) + var(dev.m1)/2

dev.m2 = -2*apply(dev.m2, 1, sum)
m2.DIC = mean(dev.m2) + var(dev.m2)/2

Gfs = c(mean(dev.m1),mean(dev.m2))
Gfs
pDs = c(var(dev.m1)/2, var(dev.m2)/2)
pDs

DICs = c(m1.DIC,m2.DIC)
DICs


# Question 5

par(mfrow=c(1,2))
### posterior distribution of the mean of sixth machine
hist(MCMC.sam.M1$theta[6,seq(n.bur+1, ns,by=thin)], freq=FALSE, main="Model 1", xlab=expression(paste(theta[6])))
abline(v=mean(MCMC.sam.M1$theta[6,seq(n.bur+1, ns,by=thin)]), col=2, lty=1)
abline(v=quantile(MCMC.sam.M1$theta[6,seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975)), col=2, lty=2)

hist(MCMC.sam.M2$theta[6,seq(n.bur+1, ns,by=thin)], freq=FALSE, main="Model 2", xlab=expression(paste(theta[6])))
abline(v=mean(MCMC.sam.M2$theta[6,seq(n.bur+1, ns,by=thin)]), col=2, lty=1)
abline(v=quantile(MCMC.sam.M2$theta[6,seq(n.bur+1, ns,by=thin)], probs=c(0.025, 0.975)), col=2, lty=2)

### posterior predictive distribution for another quality of the sixth machine
par(mfrow=c(1,2))
y.pred.1 = (rnorm(5000, MCMC.sam.M1$theta[6,seq(n.bur+1, ns,by=thin)], sqrt(MCMC.sam.M1$sig2[seq(n.bur+1, ns,by=thin)])))
hist(y.pred.1, main="Model 1", freq=FALSE, xlab=expression(paste(y[6][j])))
abline(v=mean(y.pred.1), col=2, lty=1)
abline(v=quantile(y.pred.1, probs=c(0.025, 0.975)), col=2, lty=2)

y.pred.2 = (rnorm(5000, MCMC.sam.M2$theta[6,seq(n.bur+1, ns,by=thin)], sqrt(MCMC.sam.M2$sig2[6,seq(n.bur+1, ns,by=thin)])))
hist(y.pred.2, main="Model 2", freq=FALSE, xlab=expression(paste(y[6][j])))
abline(v=mean(y.pred.2), col=2, lty=1)
abline(v=quantile(y.pred.2, probs=c(0.025, 0.975)), col=2, lty=2)

par(mfrow=c(1,2))
### posterior predictive distribution for a seventh machine
theta.pred.7th.1 = (rnorm(5000, MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)], sqrt(MCMC.sam.M1$tau[seq(n.bur+1, ns,by=thin)])))
hist(theta.pred.7th.1, main="Model 1", freq=FALSE, xlab=expression(paste(theta[7])))
abline(v=mean(theta.pred.7th.1), col=2, lty=1)
abline(v=quantile(theta.pred.7th.1, probs=c(0.025, 0.975)), col=2, lty=2)

theta.pred.7th.2 = (rnorm(5000, MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)], sqrt(MCMC.sam.M2$tau[seq(n.bur+1, ns,by=thin)])))
hist(theta.pred.7th.2, main="Model 2", freq=FALSE, xlab=expression(paste(theta[7])))
abline(v=mean(theta.pred.7th.2), col=2, lty=1)
abline(v=quantile(theta.pred.7th.2, probs=c(0.025, 0.975)), col=2, lty=2)

# prediction for a new quality of the sixth machine
theta.6th.M1 = MCMC.sam.M1$theta[6,seq(n.bur+1, ns,by=thin)]
sigma.6th.M1 = MCMC.sam.M1$sig2[seq(n.bur+1, ns,by=thin)]

theta.6th.M2 = MCMC.sam.M2$theta[6,seq(n.bur+1, ns,by=thin)]
sigma.6th.M2= MCMC.sam.M2$sig2[6,seq(n.bur+1, ns,by=thin)]

nsp = length(theta.6th.M1)
y.pred.1 = matrix(0,nsp,m6)
y.pred.2 = matrix(0,nsp,m6)

for (i in 1:nsp){
  y.pred.1[i,] = (rnorm(rep(1,m6), theta.6th.M1[i], sqrt(sigma.6th.M1[i])))
  y.pred.2[i,] = (rnorm(rep(1,m6), theta.6th.M2[i], sqrt(sigma.6th.M2[i])))
}

y.6th.M1 = (apply(y.pred.1,2,mean))
y.6th.M2 = (apply(y.pred.2,2,mean))

# prediction for a  seventh machine
theta.pred.7th.1 = (rnorm(2750, MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)], sqrt(MCMC.sam.M1$tau[seq(n.bur+1, ns,by=thin)])))
theta.pred.7th.2 = (rnorm(2750, MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)], sqrt(MCMC.sam.M2$tau[seq(n.bur+1, ns,by=thin)])))

sigma.7th.M1= MCMC.sam.M1$sig2[seq(n.bur+1, ns,by=thin)]

s02.7th.M2 = MCMC.sam.M2$s02[seq(n.bur+1, ns,by=thin)]
sigma.7th.M2 = 1/rgamma(length(s02.7th.M2),1/2,s02.7th/2)

nsp = length(theta.pred.7th.1)
y.pred.7th.1 = matrix(0,nsp,m6)
y.pred.7th.2 = matrix(0,nsp,m6)

for (i in 1:nsp){
  y.pred.7th.1[i,] = (rnorm(rep(1,m6), theta.pred.7th.1[i], sqrt(sigma.7th.M1[i])))
  y.pred.7th.2[i,] = (rnorm(rep(1,m6), theta.pred.7th.2[i], sqrt(sigma.7th.M2[i])))
}

y.7th.M1 = (apply(y.pred.7th.1,2,mean))
y.7th.M2 = (apply(y.pred.7th.2,2,mean))

plot(y.7th.M1, pch=16,ylim=c(min(y.7th.M2),max(y.7th.M2)), ylab="", main="Prediciton for a seventh machine")
points(y.7th.M2,pch=16, col=2)
legend(1,max(y.7th.M2),pch=16,col=c(1,2),legend=c("M1","M2"))
