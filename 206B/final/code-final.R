
rm(list=ls(all=TRUE))
library(dplyr)

setwd("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\winter-2018\\ams-206b\\final")

load("Machine1.RData")

y<- Machine
y<-(as.data.frame(y))

# some summaries for sample mean and sample variance
aggregate(y, by=list(y$machine), FUN=mean, na.rm=TRUE)
aggregate(y, by=list(y$machine), FUN=var, na.rm=TRUE)

# boxplot graph
boxplot(measurements ~ machine, main="Quality control measurements", xlab="Machines")

n <- 6
mm <- nrow(y)

# number of each measurement
m1=9; m2=8; m3=7; m4=6;
m5=6; m6=8;

m <- c(m1,m2,m3,m4,m5,m6)

# Question 2
# (c)

#Gibbs Sampling

## functions to sample each parameter from its full conditionals

### update theta_i
fn_update_theta <- function(y, m, mu, sig2, tau2)
{
  theta<-NULL
  for(j in 1:length(m))
  {
    #m_al <- c(0,m)
    theta[j] <- rnorm(1, (sum(y$measurements[y$machine==j])/sig2 + mu/tau2)/((m[j]/sig2)+(1/tau2)), sqrt(1/((m[j]/sig2)+(1/tau2))))
      #(sum(y[(sum(m_al[1:j],1)):(sum(m_al[1:(j+1)])),2])/sig2 + mu/tau2)
  }
  return(theta)
}

## update sig2
fn_update_sig2 <- function(mm, v0, s0, y_theta_sq)
{
  sig2 <- 1.0/rgamma(1, (mm/2)+(v0/2), (y_theta_sq)/2 + (s0^2/2))  ## mean a/b
  return(sig2)
  
}

## udpate mu
fn_update_mu <- function(theta_s, tau2, mu0, omega2, n)
{
  va <- 1/(n/tau2 + 1/omega2)
  mme <- va*((theta_s/tau2) + (mu0/omega2))
  
  mu <- rnorm(1, mme, sqrt(va))
  return(mu)
  
}

## update tau2
fn_update_tau2 <- function(n, at, y_theta_mu_sq, bt)
{
  tau2 <- 1.0/rgamma(1, (n/2)+at, y_theta_mu_sq/2 + bt)
  return(tau2)
  
}

# hyperparamenters generating higher variance on the distribution of the priors,
# being so a vague prior

## hyperparameter
hyper <- NULL
hyper$v0 <- 6
hyper$s0 <- 3

hyper$mu0 <- 90
hyper$omega2 <- 100

hyper$at <- 3
hyper$bt <- 3

## initial values
par_sam <- NULL
par_sam$theta <- rep(mean(y$measurements), n)
par_sam$sig2 <- 10^3
par_sam$mu <- 90
par_sam$tau <- 10^3

## variables for the MCMC
ns <- 40000

## save simulated para
MCMC_sam_M1 <- NULL
MCMC_sam_M1$theta <- array(NA, dim=c(n, ns))
MCMC_sam_M1$sig2 <- rep(NA, ns)
MCMC_sam_M1$mu <- rep(NA, ns)
MCMC_sam_M1$tau <- rep(NA, ns)

## MCMC modeling

for(i_iter in 1:ns)
{
  if((i_iter%%1000)==0)
  {
    print(paste("i.iter=", i_iter))
    print(date())
  }
  
  ## udpate theta
  par_sam$theta <- fn_update_theta(y, m, par_sam$mu, par_sam$sig2, par_sam$tau)
  
  #theta_vector <- function(par_theta,m){
    for(i in 1:length(m)){
      if(i==1){
        par_sam$theta_list<-rep(par_sam$theta[i],m[i])
        } else{
      p_the <- rep(par_sam$theta[i],m[i])
      par_sam$theta_list<-c(par_sam$theta_list,p_the)}
    }
    #return()}

  y$theta_list <- par_sam$theta_list
  #sum_dif <- sapply(y$machine, function(machine) { y$measurements[y$machine==machine]-par_sam$theta } )
  ## update sig2
  par_sam$sig2 <- fn_update_sig2(mm, hyper$v0, hyper$s0, (sum((y$measurements-y$theta_list)^2)))
  
  ## udpate mu
  par_sam$mu <- fn_update_mu(sum(par_sam$theta), par_sam$tau, hyper$mu0, hyper$omega2, n)
  
  ## update tau2
  par_sam$tau <- fn_update_tau2(n, hyper$at,  sum((par_sam$theta - par_sam$mu)^2), hyper$bt)
  
  ## save cur_sam
  MCMC_sam_M1$theta[,i_iter] <- par_sam$theta
  MCMC_sam_M1$sig2[i_iter] <- par_sam$sig2
  MCMC_sam_M1$mu[i_iter] <- par_sam$mu
  MCMC_sam_M1$tau[i_iter] <- par_sam$tau
    
}

n_bur <- 1000
thin <- 3

#posterior mean for the parameters
apply(MCMC_sam_M1$theta[,seq(n_bur+1, ns,by=thin)],1,mean)

mean(MCMC_sam_M1$sig2[seq(n_bur+1, ns,by=thin)])
mean(MCMC_sam_M1$mu[seq(n_bur+1, ns,by=thin)])
mean(MCMC_sam_M1$tau[seq(n_bur+1, ns,by=thin)])

# credibility intervals for the parameters
apply(MCMC_sam_M1$theta[,seq(n_bur+1, ns,by=thin)], 1, quantile,  probs=c(0.025, 0.975))

quantile(MCMC_sam_M1$sig2[seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC_sam_M1$mu[seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC_sam_M1$tau[seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975))

# Question 3

rm(list=setdiff(ls(), "MCMC_sam_M1"))
library(dplyr)

setwd("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\winter-2018\\ams-206b\\final")

load("Machine1.RData")

y<- Machine
y<-(as.data.frame(y))

n <- 6
mm <- nrow(y)

# number of each measurement
m1=9; m2=8; m3=7; m4=6;
m5=6; m6=8;

m <- c(m1,m2,m3,m4,m5,m6)


# (c)

#Gibbs Sampling

## functions to sample each parameter from its full conditionals

### update theta_i
fn_update_theta <- function(y, m, mu, sig2, tau2)
{
  theta<-NULL
  for(j in 1:length(m))
  {
    theta[j] <- rnorm(1, (sum(y$measurements[y$machine==j])/sig2[j] + mu/tau2)/((m[j]/sig2[j])+(1/tau2)), sqrt(1/((m[j]/sig2[j])+(1/tau2))))
  }
  return(theta)
}

## update sig2
fn_update_sig2 <- function(m, v0, s0, y_theta_sq_s)
{
  sigma<-NULL
  for(j in 1:length(m))
  {
    sigma[j] <- 1.0/rgamma(1, (m[j]/2)+(v0/2), (y_theta_sq_s[j])/2 + (s0^2/2)) 
  }
  return(sigma)
  
}

## update s02
fn_update_s02 <- function(n, as, sig_sum, bs)
{
  s02 <- 1.0/rgamma(1, as, (1/(sig_sum*2)) + bs)
  return(s02)
  
}

## udpate mu
fn_update_mu <- function(theta_s, tau2, mu0, omega2, n)
{
  va <- 1/(n/tau2 + 1/omega2)
  mme <- va*((theta_s/tau2) + (mu0/omega2))
  
  mu <- rnorm(1, mme, sqrt(va))
  return(mu)
  
}

## update tau2
fn_update_tau2 <- function(n, at, y_theta_mu_sq, bt)
{
  tau2 <- 1.0/rgamma(1, (n/2)+at, y_theta_mu_sq/2 + bt)
  return(tau2)
  
}

# hyperparamenters generating higher variance on the distribution of the priors,
# being so a vague prior

## hyperparameter
hyper <- NULL
hyper$v0 <- 6
hyper$s0 <- 3

hyper$mu0 <- 90
hyper$omega2 <- 100

hyper$at <- 3
hyper$bt <- 3

hyper$as <- 3
hyper$bs <- 3

## initial values
par_sam <- NULL
par_sam$theta <- rep(mean(y$measurements), n)
par_sam$sig2 <- rep(var(y$measurements), n)
par_sam$s02 <- 10^3
par_sam$mu <- 90
par_sam$tau <- 10^3

## variables for the MCMC
ns <- 40000

## save simulated para
MCMC_sam_M2 <- NULL
MCMC_sam_M2$theta <- array(NA, dim=c(n, ns))
MCMC_sam_M2$sig2 <- array(NA, dim=c(n, ns))
MCMC_sam_M2$s02 <- rep(NA, ns)
MCMC_sam_M2$mu <- rep(NA, ns)
MCMC_sam_M2$tau <- rep(NA, ns)

## MCMC modeling

for(i_iter in 1:ns)
{
  if((i_iter%%1000)==0)
  {
    print(paste("i.iter=", i_iter))
    print(date())
  }
  
  ## udpate theta
  par_sam$theta <- fn_update_theta(y, m, par_sam$mu, par_sam$sig2, par_sam$tau)
  
  ## creating a vector of theta to evaluate sigma^2
  for(ii in 1:length(m)){
    if(ii==1){
      par_sam$theta_list<-rep(par_sam$theta[ii],m[ii])
    } else{
      p_the <- rep(par_sam$theta[ii],m[ii])
      par_sam$theta_list<-c(par_sam$theta_list,p_the)}
  }
  
  y$theta_list <- par_sam$theta_list
  y_theta_dif<-NULL
  for(jj in 1:length(m)){
    y_theta_dif[jj] <- sum((y$measurements[y$machine==jj] - y$theta_list[y$machine==jj])^2)
  }
  
  par_sam$sig2 <- fn_update_sig2(m, hyper$v0, hyper$s0, y_theta_dif)
  
  par_sam$s02 <- fn_update_s02(n, hyper$as, sum(par_sam$sig2), hyper$bs)
  
  ## udpate mu
  par_sam$mu <- fn_update_mu(sum(par_sam$theta), par_sam$tau, hyper$mu0, hyper$omega2, n)
  
  ## update tau2
  par_sam$tau <- fn_update_tau2(n, hyper$at, sum((par_sam$theta - par_sam$mu)^2), hyper$bt)
  
  ## save cur_sam
  MCMC_sam_M2$theta[,i_iter] <- par_sam$theta
  MCMC_sam_M2$sig2[,i_iter] <- par_sam$sig2
  MCMC_sam_M2$mu[i_iter] <- par_sam$mu
  MCMC_sam_M2$tau[i_iter] <- par_sam$tau
  MCMC_sam_M2$s02[i_iter] <- par_sam$s02
  
}

n_bur <- 1000
thin <- 3

#posterior mean for the parameters
apply(MCMC_sam_M2$sig2[,seq(n_bur+1, ns,by=thin)],1,mean)
apply(MCMC_sam_M2$theta[,seq(n_bur+1, ns,by=thin)],1,mean)

mean(MCMC_sam_M2$mu[seq(n_bur+1, ns,by=thin)])
mean(MCMC_sam_M2$tau[seq(n_bur+1, ns,by=thin)])
mean(MCMC_sam_M2$s02[seq(n_bur+1, ns,by=thin)])

# credibility intervals for the parameters
apply(MCMC_sam_M2$theta[,seq(n_bur+1, ns,by=thin)], 1, quantile,  probs=c(0.025, 0.975))
apply(MCMC_sam_M2$sig2[,seq(n_bur+1, ns,by=thin)], 1, quantile,  probs=c(0.025, 0.975))

quantile(MCMC_sam_M2$mu[seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC_sam_M2$tau[seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975))
quantile(MCMC_sam_M2$s02[seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975))

# Question 4

# Computing deviance information criteria
ndic = length(MCMC_sam_M2$mu[seq(n_bur+1, ns,by=thin)])
l <- vector("list", 6)
for(i in 1:6){
  theta_m1 = matrix(rep(MCMC_sam_M1$theta[i,seq(n_bur+1, ns,by=thin)],m[i]), nrow = ndic, ncol = m[i])
  l[[i]] = theta_m1 
}
theta_m1 = do.call("cbind", l)
sigma_m1 = (MCMC_sam_M1$sig2[seq(n_bur+1, ns,by=thin)])

dev.m1 = matrix(0, ndic, nrow(y))
for (i in 1:nrow(y)){
  dev.m1[,i] = dnorm(y$measurements[i], theta_m1[,i], sqrt(sigma_m1), log = TRUE)
}

l <- vector("list", 6)
for(i in 1:6){
  theta_m2 = matrix(rep(MCMC_sam_M2$theta[i,seq(n_bur+1, ns,by=thin)],m[i]), nrow = ndic, ncol = m[i])
  l[[i]] = theta_m2
}
theta_m2 = do.call("cbind", l)

l <- vector("list", 6)
for(i in 1:6){
  sigma_m2 = matrix(rep(MCMC_sam_M2$sig2[i,seq(n_bur+1, ns,by=thin)],m[i]), nrow = ndic, ncol = m[i])
  l[[i]] = sigma_m2 
}
sigma_m2 = do.call("cbind", l)

dev.m2 = matrix(0, ndic, nrow(y))
for (i in 1:nrow(y)){
  dev.m2[,i] = dnorm(y$measurements[i], theta_m2[,i], sqrt(sigma_m2[,i]), log = TRUE)
}

dev.m1 = -2*apply(dev.m1, 1, sum)
m1_DIC = mean(dev.m1) + var(dev.m1)/2

dev.m2 = -2*apply(dev.m2, 1, sum)
m2_DIC = mean(dev.m2) + var(dev.m2)/2

Gfs = c(mean(dev.m1),mean(dev.m2))
Gfs
pDs = c(var(dev.m1)/2, var(dev.m2)/2)
pDs

DICs = c(m1_DIC,m2_DIC)
DICs


# Question 5

par(mfrow=c(1,2))
### posterior distribution of the mean of sixth machine
hist(MCMC_sam_M1$theta[6,seq(n_bur+1, ns,by=thin)], freq=FALSE, main="Model 1", xlab=expression(paste(theta[6])))
abline(v=mean(MCMC_sam_M1$theta[6,seq(n_bur+1, ns,by=thin)]), col=2, lty=1)
abline(v=quantile(MCMC_sam_M1$theta[6,seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975)), col=2, lty=2)

hist(MCMC_sam_M2$theta[6,seq(n_bur+1, ns,by=thin)], freq=FALSE, main="Model 2", xlab=expression(paste(theta[6])))
abline(v=mean(MCMC_sam_M2$theta[6,seq(n_bur+1, ns,by=thin)]), col=2, lty=1)
abline(v=quantile(MCMC_sam_M2$theta[6,seq(n_bur+1, ns,by=thin)], probs=c(0.025, 0.975)), col=2, lty=2)

### posterior predictive distribution for another quality of the sixth machine
par(mfrow=c(1,2))
y_pred_1 <- (rnorm(5000, MCMC_sam_M1$theta[6,seq(n_bur+1, ns,by=thin)], sqrt(MCMC_sam_M1$sig2[seq(n_bur+1, ns,by=thin)])))
hist(y_pred_1, main="Model 1", freq=FALSE, xlab=expression(paste(y[6][j])))
abline(v=mean(y_pred_1), col=2, lty=1)
abline(v=quantile(y_pred_1, probs=c(0.025, 0.975)), col=2, lty=2)

y_pred_2 <- (rnorm(5000, MCMC_sam_M2$theta[6,seq(n_bur+1, ns,by=thin)], sqrt(MCMC_sam_M2$sig2[6,seq(n_bur+1, ns,by=thin)])))
hist(y_pred_2, main="Model 2", freq=FALSE, xlab=expression(paste(y[6][j])))
abline(v=mean(y_pred_2), col=2, lty=1)
abline(v=quantile(y_pred_2, probs=c(0.025, 0.975)), col=2, lty=2)

par(mfrow=c(1,2))
### posterior predictive distribution for a seventh machine
theta_pred_7th_1 <- (rnorm(5000, MCMC_sam_M1$mu[seq(n_bur+1, ns,by=thin)], sqrt(MCMC_sam_M1$tau[seq(n_bur+1, ns,by=thin)])))
hist(theta_pred_7th_1, main="Model 1", freq=FALSE, xlab=expression(paste(theta[7])))
abline(v=mean(theta_pred_7th_1), col=2, lty=1)
abline(v=quantile(theta_pred_7th_1, probs=c(0.025, 0.975)), col=2, lty=2)

theta_pred_7th_2 <- (rnorm(5000, MCMC_sam_M2$mu[seq(n_bur+1, ns,by=thin)], sqrt(MCMC_sam_M2$tau[seq(n_bur+1, ns,by=thin)])))
hist(theta_pred_7th_2, main="Model 2", freq=FALSE, xlab=expression(paste(theta[7])))
abline(v=mean(theta_pred_7th_2), col=2, lty=1)
abline(v=quantile(theta_pred_7th_2, probs=c(0.025, 0.975)), col=2, lty=2)

# prediction for a new quality of the sixth machine
theta_6th_M1 = MCMC_sam_M1$theta[6,seq(n_bur+1, ns,by=thin)]
sigma_6th_M1 = MCMC_sam_M1$sig2[seq(n_bur+1, ns,by=thin)]

theta_6th_M2 = MCMC_sam_M2$theta[6,seq(n_bur+1, ns,by=thin)]
sigma_6th_M2= MCMC_sam_M2$sig2[6,seq(n_bur+1, ns,by=thin)]

nsp = length(theta_6th_M1)
y_pred_1 = matrix(0,nsp,m6)
y_pred_2 = matrix(0,nsp,m6)

for (i in 1:nsp){
  y_pred_1[i,] = (rnorm(rep(1,m6), theta_6th_M1[i], sqrt(sigma_6th_M1[i])))
  y_pred_2[i,] = (rnorm(rep(1,m6), theta_6th_M2[i], sqrt(sigma_6th_M2[i])))
}

y_6th_M1 = (apply(y_pred_1,2,mean))
y_6th_M2 = (apply(y_pred_2,2,mean))

# prediction for a  seventh machine
theta_pred_7th_1 = (rnorm(2750, MCMC_sam_M1$mu[seq(n_bur+1, ns,by=thin)], sqrt(MCMC_sam_M1$tau[seq(n_bur+1, ns,by=thin)])))
theta_pred_7th_2 = (rnorm(2750, MCMC_sam_M2$mu[seq(n_bur+1, ns,by=thin)], sqrt(MCMC_sam_M2$tau[seq(n_bur+1, ns,by=thin)])))

sigma_7th_M1= MCMC_sam_M1$sig2[seq(n_bur+1, ns,by=thin)]

s02_7th_M2 = MCMC_sam_M2$s02[seq(n_bur+1, ns,by=thin)]
sigma_7th_M2 = 1/rgamma(length(s02_7th_M2),1/2,s02_7th/2)

nsp = length(theta_pred_7th_1)
y_pred_7th_1 = matrix(0,nsp,m6)
y_pred_7th_2 = matrix(0,nsp,m6)

for (i in 1:nsp){
  y_pred_7th_1[i,] = (rnorm(rep(1,m6), theta_pred_7th_1[i], sqrt(sigma_7th_M1[i])))
  y_pred_7th_2[i,] = (rnorm(rep(1,m6), theta_pred_7th_2[i], sqrt(sigma_7th_M2[i])))
}

y_7th_M1 = (apply(y_pred_7th_1,2,mean))
y_7th_M2 = (apply(y_pred_7th_2,2,mean))

plot(y_7th_M1, pch=16,ylim=c(min(y_7th_M2),max(y_7th_M2)), ylab="", main="Prediciton for a seventh machine")
points(y_7th_M2,pch=16, col=2)
legend(1,max(y_7th_M2),pch=16,col=c(1,2),legend=c("M1","M2"))
