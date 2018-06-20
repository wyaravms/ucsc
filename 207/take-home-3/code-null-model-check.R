
rm(list=ls(all=TRUE))
set.seed(7)
library(mvtnorm)

sst = read.csv("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\spring-2018\\ams-207\\take-home-3\\sst.csv")
attach(sst)
Type = factor(Type)

nd = N
n = nrow(sst)
y = temp
X = cbind(rep(1,n),lat,lon,Type)
ncov = ncol(X)

# exploratory plots
pairs(sst[,-5], pch=16, panel=panel.smooth, cex.axis=1.5, cex.lab=1.5)

# fitting a null model (M2) to compare with the main model (M1)
fit.null = lm(temp ~ 1, data=sst, x=TRUE, y=TRUE)
X.null = fit.null$x
ncov = ncol(X.null)
summary(fit.null)

# gibbs steps
nmc = 41000

# hyperparameters
alpha = 3; beta = 3
a = 2; b =2

post.mu.i.null = array(NA, dim=c(nmc, n))
post.sigma.i.null = array(NA, dim=c(nmc, n))
post.beta.null = array(NA, dim=c(nmc, ncov))
post.sigma.null = array(NA, nmc)
post.tau.null = array(NA, nmc)

# initial values
post.tau.null[1] = 1.5
post.mu.i.null[1,] = y
post.sigma.i.null[1,] = 1.5
post.sigma.null[1] = 1.5

for (i in 2:nmc){
  
  m.beta.null = solve(t(X.null)%*%X.null)%*%t(X.null)%*%post.mu.i.null[i-1,]
  var.beta.null = post.tau.null[i-1]*solve(t(X.null)%*%X.null)
  
  post.beta.null[i,] = rmvnorm(1, m.beta.null, var.beta.null)
  
  a.tau.null = (n)/2 - 1 
  b.tau.null = sum((post.mu.i.null[i-1,] - X.null%*%post.beta.null[i,])^2)/2
  
  sigma.hat.null = sum((y-X.null%*%m.beta.null)^2)/(n-qr(X.null)$rank)
  
  #post.tau[i] = 1/rgamma(1, ((n-qr(X)$rank)/2), (n-qr(X)$rank)*sigma.hat/2)
  post.tau.null[i] = 1/rgamma(1, a.tau.null, b.tau.null)
  
  m.mui.null = ((y*post.tau.null[i]*nd) + (X.null%*%post.beta.null[i,]*post.sigma.i.null[i-1,]))/(post.tau.null[i]*nd + post.sigma.i.null[i-1,])
  var.mui.null = (post.sigma.i.null[i-1,]*post.tau.null[i])/((post.tau.null[i]*nd)+(post.sigma.i.null[i-1,]))
  
  post.mu.i.null[i,] = rnorm(n,as.vector(m.mui.null),sqrt(var.mui.null))
  
  a.sigmai.null = 1/2 + (alpha + 1)
  b.sigmai.null = ((((y - post.mu.i.null[i,])^2)*nd)/2) + alpha*post.sigma.null[i-1]
  
  post.sigma.i.null[i,] = 1/rgamma(n,a.sigmai.null,b.sigmai.null)
  
  a.sigma.null = n*(alpha + 1) + a
  b.sigma.null = (alpha*sum(1/(post.sigma.i.null[i,]))) + b
  
  post.sigma.null[i] = rgamma(1, a.sigma.null, b.sigma.null)
  
  cat(i, "/", nmc, "\r")
}

n.bur = 1000; thin = 4

post.beta.null = post.beta.null[seq(n.bur+1, nmc,by=thin),]
post.tau.null = post.tau.null[seq(n.bur+1, nmc,by=thin)]
post.mu.i.null = post.mu.i.null[seq(n.bur+1, nmc,by=thin),]
post.sigma.i.null = post.sigma.i.null[seq(n.bur+1, nmc,by=thin),]
post.sigma.null = post.sigma.null[seq(n.bur+1, nmc,by=thin)]

nmcb = length(post.sigma.null)

par(mfrow=c(2,3))
plot.ts(post.beta.null)
plot.ts(post.tau.null)
plot.ts(post.mu.i.null[,1])
plot.ts(post.sigma.i.null[,1])
plot.ts(post.sigma.null)

# Posterior predictions
pred.y.m2 = matrix(0, n, nmcb)
for (i in 1:nmcb)
  pred.y.m2[,i] = rnorm(n, post.mu.i.null[i,], sqrt(post.sigma.i.null[i,]/N))

#Gelfand and Ghosh
G1 = sum((apply(pred.y, 1, mean)-y)^2)
P1 = sum(apply(pred.y,1,var))
D1G = G1 + P1

G2 = sum((apply(pred.y.m2, 1,mean)-y)^2)
P2 = sum(apply(pred.y.m2,1,var))
D2G = G2 + P2

GGs = c(D1G,D2G)
GGs

# DIC (deviance information criteria)
m1 = matrix(0, nmcb, length(y)) # posterior distribution from the main model
for (i in 1:length(y))
  m1[,i] = dnorm(y[i], post.mu.i[,i], sqrt(post.sigma.i[,i]/N[i]), log = TRUE)

m2 = matrix(0, nmcb, length(y))
for (i in 1:length(y))
  m2[,i] = dnorm(y[i], post.mu.i.null[,i], sqrt(post.sigma.i.null[,i]/N[i]), log = TRUE)

m1.m = -2*apply(m1, 1, sum)
m1.DIC = mean(m1.m) + var(m1.m)/2

m2.m = -2*apply(m2, 1, sum)
m2.DIC = mean(m2.m) + var(m2.m)/2

pDs = c(var(m1.m)/2, var(m2.m)/2)
pDs

DICs = c(m1.DIC,m2.DIC)
DICs

# DIC (deviance information criteria)
mean.mu.i.m1 = apply(post.mu.i,2,mean)
mean.sigma.i.m1 = apply(post.sigma.i,2,mean)

mean.mu.i.m2 = apply(post.mu.i.null,2,mean)
mean.sigma.i.m2 = apply(post.sigma.i.null,2,mean)

# deviance statistics
mm1 = matrix(0, 1, n)
for (i in 1:n)
  mm1[,i] = dnorm(y[i], mean.mu.i.m1[i], sqrt(mean.sigma.i.m1[i]/N[i]), log = TRUE)

mm2 = matrix(0, 1, n)
for (i in 1:n)
  mm2[,i] = dnorm(y[i], mean.mu.i.m2[i], sqrt(mean.sigma.i.m2[i]/N[i]), log = TRUE)

mm1 = matrix(0, 1, n)
mm1 = dnorm(y,mean.mu.i.m1, sqrt(mean.sigma.i.m1/N), log = TRUE)

mm2 = matrix(0, 1, n)
mm2 = dnorm(y,mean.mu.i.m2, sqrt(mean.sigma.i.m2/N), log = TRUE)

mm1 = sum(mm1)
mm2 = sum(mm2)

m1 = matrix(0, nmcb, n)
for (i in 1:nmcb)
  m1[i,] = dnorm(y, post.mu.i[i,], sqrt(post.sigma.i[i,]/N), log = TRUE)

m2 = matrix(0, nmcb, n)
for (i in 1:nmcb)
  m2[i,] = dnorm(y, post.mu.i.null[i,], sqrt(post.sigma.i.null[i,]/N), log = TRUE)

pD1 = 2*(mm1 - mean(m1)) 
pD2 = 2*(mm2 - mean(m2))

DIC1 = -2*mm1 + 2*pD1  
DIC2 = -2*mm2 + 2*pD2 

pDs = c(pD1,pD2)
DICs = c(DIC1,DIC2)

DICs

# DIC (deviance information criteria) from the slides
likelihood = function(y, mu.dic, sigma.dic){
  val = sum(dnorm(y, mu.dic, sigma.dic, log = TRUE))
  return(val)
}

hlp = NULL
for(t in 1:nmcb){
  hlp[t] = likelihood(y, post.mu.i[i,], sqrt(post.sigma.i[i,]/N))
}
mu.mean.dic = apply(post.mu.i, 2, mean)
sigma.mean.dic = apply(post.sigma.i, 2, mean)

lph = likelihood(y, mu.mean.dic, sqrt(sigma.mean.dic/N))
pdic=2*(lph-mean(hlp))
DIC.1=-2*lph+2*pdic

hlp.2 = NULL
for(t in 1:nmcb){
  hlp.2[t] = likelihood(y, post.mu.i.null[i,], sqrt(post.sigma.i.null[i,]/N))
}
mu.mean.dic.2 = apply(post.mu.i.null, 2, mean)
sigma.mean.dic.2 = apply(post.sigma.i.null, 2, mean)

lph.2 = likelihood(y, mu.mean.dic.2, sqrt(sigma.mean.dic.2/N))
pdic.2 = 2*(lph.2-mean(hlp.2))
DIC.2 = -2*lph.2+2*pdic.2

c(DIC.1,DIC.2)
