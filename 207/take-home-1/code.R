
# AMS 207 - take home 1

rm(list=ls(all=TRUE))
set.seed(7)
require(lme4)
require(invgamma)

data(Dyestuff)
Dyestuff$Batch = as.numeric(Dyestuff$Batch)

N = 6
n = 5
Nn = N*n

dyes = (as.data.frame(Dyestuff))
aggregate(dyes, by=list(dyes$Batch), FUN=mean)
aggregate(dyes, by=list(dyes$Batch), FUN=var)

par(mfrow=c(1,1))
boxplot(split(dyes$Yield ,dyes$Batch), ylab="Yield (in grams)",xlab="Batch")
points(dyes)


# gibbs sampling

## functions to sample each parameter from its full conditionals

### update b.i
fn.update.b = function(y, N, n, mu, sig2y, sig2b)
{
  b<-NULL
  for(j in 1:N)
  {
    b[j] = rnorm(1, ((sum(dyes$Yield[dyes$Batch==j] - mu)/sig2y)*(1/((n/sig2y)+(1/sig2b)))), sqrt(1/((n/sig2y)+(1/sig2b))))
  }
  return(b)
}

## udpate mu
fn.update.mu = function(y.b, sig2y,Nn)
{
  va = sig2y/Nn
  mme = va*(y.b/sig2y)
  
  mu = rnorm(1, mme, sqrt(va))
  return(mu)
}

## update sig2y
fn.update.sig2y = function(Nn, y.mu.b.sq)
{
  sigy = rinvgamma(1, shape=((Nn-2)/2.0), scale=1/(y.mu.b.sq/2.0))  ## mean a/b
  return(sigy)
}

## update sig2b
fn.update.sig2b = function(N, y.b.sq)
{
  sig2b = rinvgamma(1, shape=((N-2)/2.0), scale=1/(y.b.sq/2.0))
  return(sig2b)
}

## initial values
par.sam = NULL
par.sam$b = rep(0, N)
par.sam$mu = 1000
par.sam$sig2y = 10^2
par.sam$sig2b = 10^2

## variables for the MCMC
ns = 41000

## save simulated para
MCMC.sam.M1 = NULL
MCMC.sam.M1$b = array(NA, dim=c(N, ns))
MCMC.sam.M1$mu = rep(NA, ns)
MCMC.sam.M1$sig2y = rep(NA, ns)
MCMC.sam.M1$sig2b = rep(NA, ns)

## MCMC modeling

for(i.iter in 1:ns)
{
  if((i.iter%%1000)==0)
  {
    print(paste("i.iter=", i.iter))
    print(date())
  }
  
  ## udpate theta
  par.sam$b = fn.update.b(dyes, N, n, par.sam$mu, par.sam$sig2y, par.sam$sig2b)
  
  for(i in 1:N){
    if(i==1){
      par.sam$b.list = rep(par.sam$b[i],n)
    } else{
      p.the = rep(par.sam$b[i],n)
      par.sam$b.list = c(par.sam$b.list,p.the)}
  }
  
  dyes$b.list = par.sam$b.list
  
  par.sam$mu = fn.update.mu((sum(dyes$Yield - dyes$b.list)), par.sam$sig2y, Nn)
  
  ## update sigma y
  par.sam$sig2y = fn.update.sig2y(Nn, sum((dyes$Yield - par.sam$mu - dyes$b.list)^2))
  
  ## udpate sigma b
  par.sam$sig2b = fn.update.sig2b(N, sum(par.sam$b^2))
  
  ## save cur.sam
  MCMC.sam.M1$b[,i.iter] = par.sam$b
  MCMC.sam.M1$mu[i.iter] = par.sam$mu
  MCMC.sam.M1$sig2y[i.iter] = par.sam$sig2y
  MCMC.sam.M1$sig2b[i.iter] = par.sam$sig2b
}

n.bur = 1000
thin = 10

#posterior mean for the parameters
mean(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)])
mean(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])
mean(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])
mean(sqrt(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)]))
mean(sqrt(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)]))

apply(MCMC.sam.M1$b[,seq(n.bur+1, ns,by=thin)], 1, mean)

#credibility intervals
quantile(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)], probs=c(0.025, .5, 0.975))
quantile(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)], probs=c(0.025, .5, 0.975))
quantile(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)], probs=c(0.025, .5, 0.975))
quantile(sqrt(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)]), probs=c(0.025, .5, 0.975))
quantile(sqrt(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)]), probs=c(0.025, .5, 0.975))

apply(MCMC.sam.M1$b[,seq(n.bur+1, ns,by=thin)], 1, quantile, probs=c(.025,.5, .975), na.rm=TRUE)

# the crude densities for mu and sigma2
par(mfrow=c(1,1))
plot(density(MCMC.sam.M1$b[1,seq(n.bur+1, ns,by=thin)]),cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, xlab="", main="")
lines(density(MCMC.sam.M1$b[2,seq(n.bur+1, ns,by=thin)]), col=2)
lines(density(MCMC.sam.M1$b[3,seq(n.bur+1, ns,by=thin)]), col=3)
lines(density(MCMC.sam.M1$b[4,seq(n.bur+1, ns,by=thin)]), col=4)
lines(density(MCMC.sam.M1$b[5,seq(n.bur+1, ns,by=thin)]), col=5)
lines(density(MCMC.sam.M1$b[6,seq(n.bur+1, ns,by=thin)]), col=6)

legend(-260,0.011, lty = 1, cex=1.2,box.lty = 0, lwd = 1, legend=c(expression(b[1]), expression(b[2]), expression(b[3]), expression(b[4]), expression(b[5]), expression(b[6])), col=1:6)

plot.ts(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)])
plot.ts(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])
plot.ts(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])

acf(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)])
acf(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])
acf(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])

par(mfrow=c(3,1))
plot(density(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)]), cex.lab=1.5,cex.main = 3,cex.axis=2,main = expression(paste(mu)), xlab="")
plot(density(sqrt(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])),cex.lab=1.5,cex.main = 3,cex.axis=2, xlab="", main = expression(paste(sigma[y])))
plot(density(sqrt(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])),xlim=c(0,300), cex.lab=1.5,cex.main = 3,cex.axis=2,xlab="", main = expression(paste(sigma[b])))

par(mfrow=c(3,1))
plot(density(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)]), cex.lab=1.5,cex.main = 3,cex.axis=2,main = expression(paste(mu)), xlab="")
plot(density(log(sqrt(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)]))), cex.lab=1.5,cex.main =3,cex.axis=2, xlab="", main = expression("log(" * sigma[y] * ")"))
plot(density(log(sqrt(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)]))), cex.lab=1.5,cex.main = 3,cex.axis=2, xlab="", main = expression("log(" * sigma[b] * ")"))


# graphs for densities for the mu and transformation log(sigma)
# just to compare with the results given by the normal approximation
# for the transformed parameters
par(mfrow=c(1,3))
plot.ts(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)])
plot.ts(log(sqrt(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])))
plot.ts(log(sqrt(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])))

par(mfrow=c(1,3))
plot(density(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)]), main = expression(paste(mu)))
plot(density(log(sqrt(MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)]))), main = expression(paste(sigma[y])))
plot(density(log(sqrt(MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)]))), main = expression(paste(sigma[b])))

# part 4
rm(list=ls(all=TRUE))
library(LearnBayes)
library(mvtnorm)
data=matrix(c(1545,1440,1440,1520,1580,
              1540,1555,1490,1560,1495,
              1595,1550,1605,1510,1560,
              1445,1440,1595,1465,1545,
              1595,1630,1515,1635,1625,
              1520,1455,1450,1480,1445),c(6,5),byrow=TRUE)

log.post.var.comp=function(theta,y)
{
  mu = theta[1]; sigma.y = exp(theta[2]); sigma.b = exp(theta[3])
  Y=apply(y,1,mean); n=dim(y)[2]
  S=apply(y,1,var)*(n-1)
  loglike=sum(dnorm(Y,mu,sqrt(sigma.y^2/n+sigma.b^2),log=TRUE)+
                dgamma(S,shape=(n-1)/2,rate=1/(2*sigma.y^2),log=TRUE))
  return(loglike+theta[2]+theta[3])
}

logf = function(theta, data) {
  if (is.matrix(theta) == TRUE) {
    val = matrix(0, c(dim(theta)[1], 1))
    for (j in 1:dim(theta)[1]) val[j] = log.post.var.comp(theta[j,], data)
  }
  else val = log.post.var.comp(theta, data)
  return(val)
}

limits=c(-5,10,-5,10)
ng=100
z0 = mean1[1]
x0 = seq(limits[1], limits[2], len = ng)
y0 = seq(limits[3], limits[4], len = ng)
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den = logf(cbind(z0,X[1:n2], Y[1:n2]), data)
den = den - max(den)
den = matrix(den, c(ng, ng))
contour(x0, y0, den,levels = seq(-8.9, 0, by = 2.3),xlim=c(3.1,4.8),ylim=c(-4,6.5),cex.lab=1.5,cex.main = 3,cex.axis=1.8,lwd = 1,
        xlab=expression(sigma[y]),ylab=expression(sigma[b]),add=TRUE,cex.axis=1.2, cex.lab=1.2)


limits=c(1200,1700,-5,10)
ng=100
x0 = seq(limits[1], limits[2], len = ng)
y0 = seq(limits[3], limits[4], len = ng)
z0 = mean1[3]
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den11 = logf(cbind(X[1:n2], Y[1:n2],z0), data)
den11 = den11 - max(den11)
den11 = matrix(den11, c(ng, ng))
contour(x0, y0, den11, levels = seq(-8.9, 0, by = 2.3), xlim=c(1400,1650),ylim=c(2.5,5.5),cex.lab=1.5,cex.main = 3,cex.axis=1.8,lwd = 1,
        xlab=expression(mu),ylab=expression(sigma[y]),cex.axis=1.2,cex.lab=1.2, col=1)

limits=c(1200,1900,-5,10)
ng=100
x0 = seq(limits[1], limits[2], len = ng)
z0 = mean1[2]
y0 = seq(limits[3], limits[4], len = ng)
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den22 = logf(cbind(X[1:n2],z0,Y[1:n2]), data)
den22 = den22 - max(den22)
den22 = matrix(den22, c(ng, ng))
contour(x0, y0, den22, levels = seq(-8.9, 0, by = 2.3), xlim=c(1250,1800),ylim=c(-1.5,7),cex.lab=1.5,cex.main = 3,cex.axis=1.8,lwd = 1,
        xlab=expression(mu),ylab=expression(sigma[b]),add=TRUE,cex.axis=1.2, cex.lab=1.2, col=1)
legend(1250,1.8,lty = 1, cex=1, box.lty = 0, lwd = 1, legend=c("Exact","Normal"), col=c(1,2))

start = c(1500, 3, 3)
theta.mle = optim(par=start, log.post.var.comp, y=data, method="BFGS",hessian = TRUE, control = list(fnscale = -1))
theta.mle$par
theta.mle$hessian

fisher = solve(-theta.mle$hessian)
fisher

mean1 = c(theta.mle$par[1],theta.mle$par[2],theta.mle$par[3])
Var1  = fisher
tpar=list(m=mean1,var=Var1,df=2)

logf.n = function(theta, tpar) {
  if (is.matrix(theta) == TRUE) {
    val = matrix(0, c(dim(theta)[1], 1))
    for (j in 1:dim(theta)[1]) val[j] = dmvnorm(theta[j,], mean = tpar$m, sigma = tpar$var,log=TRUE)
  } else val = dmvnorm(theta, mean = tpar$m, sigma = tpar$var,log=TRUE)
  return(val)
}

limits=c(-5,10,-5,10)
ng=100
z0 = mean1[1]
x0 = seq(limits[1], limits[2], len = ng)
y0 = seq(limits[3], limits[4], len = ng)
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den1 = logf.n(cbind(z0, X[1:n2], Y[1:n2]), tpar)
den1 = den1 - max(den1)
den1 = matrix(den1, c(ng, ng))
contour(x0, y0, den1, levels = seq(-8.9, 0, by = 2.3), xlim=c(3.1,4.8),ylim=c(1,6.5),lwd = 1,
        xlab=expression(sigma[y]),ylab=expression(sigma[b]),add = TRUE,cex.axis=1.2, cex.lab=1.2, col=2)
legend(3.1,0,lty = 1, cex=1, box.lty = 0, lwd = 1, legend=c("Exact","Normal"), col=c(1,2))


limits=c(1200,1700,-5,10)
ng=100
x0 = seq(limits[1], limits[2], len = ng)
y0 = seq(limits[3], limits[4], len = ng)
z0 = mean1[3]
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den2 = logf.n(cbind(X[1:n2], Y[1:n2],z0), tpar)
den2 = den2 - max(den2)
den2 = matrix(den2, c(ng, ng))
contour(x0, y0, den2, levels = seq(-8.9, 0, by = 2.3), xlim=c(1400,1650),ylim=c(1.5,6),lwd = 1,
        xlab=expression(mu),ylab="",cex.axis=1.2, add=TRUE,cex.lab=1.2, col=2)

legend(1402,3.3,lty = 1, cex=1, box.lty = 0, lwd = 1, legend=c("Exact","Normal"), col=c(1,2))

limits=c(1200,1700,-5,10)
ng=100
x0 = seq(limits[1], limits[2], len = ng)
z0 = mean1[2]
y0 = seq(limits[3], limits[4], len = ng)
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den3 = logf.n(cbind(X[1:n2],z0,Y[1:n2]), tpar)
den3 = den3 - max(den3)
den3 = matrix(den3, c(ng, ng))
contour(x0, y0, den3, levels = seq(-19.9, 0, by = 6.3), xlim=c(1450,1600),ylim=c(1.5,6),lwd = 1,
        xlab=expression(mu),ylab=expression(sigma[b]),add=TRUE,cex.axis=1.2, cex.lab=1.2, col=2)

legend(1400,5.8, lty = 1, cex=1.2, legend=c(expression(sigma[y]), expression(sigma[b])), col=1:2)


# reject sampling
log.rej = function(theta,y,tpar){
  l.rej = log.post.var.comp(theta,y) - dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df,log=TRUE) 
  return(l.rej)
}

tpar=list(m=mean1,var=Var1,df=3)

start = c(1500, 3, 3)
theta.rej = optim(par=start, log.rej, y=data, tpar=tpar,method="BFGS",hessian = TRUE, control = list(fnscale = -1))
theta.rej$par

dmax=log.rej(theta.rej$par,y=data,tpar=tpar)
dmax

# Accept/reject method
set.seed(7)
M = 55000
d = length(tpar$m)
theta = rmt(M,mean = c(tpar$m), S = tpar$var, df = tpar$df)
lf = matrix(0, c(dim(theta)[1], 1))

for (j in 1:dim(theta)[1]) lf[j] = log.post.var.comp(theta[j,],y=data)

lg = dmt(theta,  mean = c(tpar$m), S = tpar$var, df = tpar$df, log = TRUE)
if (d == 1) {
  prob = exp(c(lf) - lg - dmax)
  draws1=(theta[runif(M) < prob])
} else {
  prob = exp(lf - lg - dmax)
  draws1=(theta[runif(M) < prob, ])
}

nd = nrow(draws1)
accept.rate = nd/M
accept.rate  

plot(draws1[,1],draws1[,2],
     xlab=expression(mu),ylab=expression(sigma[y]),xlim=c(1450,1600),ylim=c(1.5,6),pch=1)

plot(draws1[,1],draws1[,3],
     xlab=expression(mu),ylab=expression(sigma[b]),xlim=c(1350,1800),ylim=c(-1.5,7),pch=1)

plot(draws1[,2],draws1[,3],
     xlab=expression(sigma[y]),ylab=expression(sigma[b]),xlim=c(3.1,4.8),ylim=c(-4,7),pch=1)


# SIR method
set.seed(7)
sir.sexp = function (logf, tpar, n, infor) 
{
  k = length(tpar$m)
  theta = rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
  lf = matrix(0, c(dim(theta)[1], 1))
  for (j in 1:dim(theta)[1]) lf[j] = logf(theta[j, ], infor)
  lp = dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
           log = TRUE)
  md = max(lf - lp)
  wt = exp(lf - lp - md)
  probs = wt/sum(wt)
  indices = sample(1:n, size = n, prob = probs, replace = TRUE)
  if (k > 1) 
    theta = theta[indices, ]
  else theta = theta[indices]
  return(theta)
}

theta=sir.sexp(log.post.var.comp,tpar,nd,data)

plot(theta[,1],theta[,2],
     xlab=expression(mu),ylab=expression(sigma[y]),xlim=c(1450,1600),ylim=c(1.5,6),pch=1)

plot(theta[,1],theta[,3],
     xlab=expression(mu),ylab=expression(sigma[b]),xlim=c(1450,1600),ylim=c(-1.5,8),pch=1)

plot(theta[,2],theta[,3],
     xlab=expression(sigma[y]),ylab=expression(sigma[b]),xlim=c(3.1,4.8),ylim=c(-4,7),pch=1)

par(mfrow=c(1,3))
plot(density(draws1[,1]),axes = FALSE, main="", xlab="", col=4)
lines(density(theta[,1]), col=2)
lines(density(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)]), col=3)
legend("topright",lty = 1, cex=1,box.lty = 0, lwd = 1, legend=c(expression(mu[Gibbs]),expression(mu[RS]), expression(mu[SIR])), col=c(3,4,2))
axis(1)
axis(2)

plot(density(draws1[,2]),axes = FALSE, main="", xlab="", col=4)
lines(density(theta[,2]), col=2)
lines(density(log(sqrt((MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])))), col=3)
legend("topleft",lty = 1, cex=1,box.lty = 0, lwd = 1, legend=c(expression(sigma[y][Gibbs]),expression(sigma[y][RS]), expression(sigma[y][SIR])), col=c(3,4,2))
axis(1)
axis(2)

plot(density(draws1[,3]),axes = FALSE, main="", xlab="", col=4)
lines(density(theta[,3]), col=2)
lines(density(log(sqrt((MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])))), col=3)
legend("topleft",lty = 1, cex=1, box.lty = 0, lwd = 1, legend=c(expression(sigma[b][Gibbs]),expression(sigma[b][RS]), expression(sigma[b][SIR])), col=c(3,4,2))
axis(1)
axis(2)

#par(mfrow=c(1,3))
plot(density(draws1[,1]), main="", xlab="", cex.lab=1.5,cex.main = 3,cex.axis=2,col=4)
lines(density(theta[,1]), col=2)
lines(density(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)]), col=3)
legend(1300,0.0174,lty = 1, cex=1.3,box.lty = 0, lwd = 1, legend=c(expression(mu[Gibbs]),expression(mu[RS]), expression(mu[SIR])), col=c(3,4,2))

plot(density(draws1[,2]), main="", xlab="", cex.lab=1.5,cex.main = 3,cex.axis=2,col=4)
lines(density(theta[,2]), col=2)
lines(density(log(sqrt((MCMC.sam.M1$sig2y[seq(n.bur+1, ns,by=thin)])))), col=3)
legend(3.4,2.7,lty = 1, cex=1.3,box.lty = 0, lwd = 1, legend=c(expression(sigma[y][Gibbs]),expression(sigma[y][RS]), expression(sigma[y][SIR])), col=c(3,4,2))

plot(density(draws1[,3]), main="",xlim=c(0,6),cex.lab=1.5,cex.main = 3,cex.axis=2, xlab="", col=4)
lines(density(theta[,3]), col=2)
lines(density(log(sqrt((MCMC.sam.M1$sig2b[seq(n.bur+1, ns,by=thin)])))), col=3)
legend(0.1,0.78,lty = 1, cex=1.3, box.lty = 0, lwd = 1, legend=c(expression(sigma[b][Gibbs]),expression(sigma[b][RS]), expression(sigma[b][SIR])), col=c(3,4,2))
