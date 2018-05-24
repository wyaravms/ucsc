
# HW 2 - AMS 207
rm(list=ls(all=TRUE))
set.seed(7)

# Question 8

# number of deaths, number of occurrences 
# Annual deaths from horse kicks in the Prussian army (1875-1894)

# data
n = 20
data = c(3, 5, 7, 9, 10, 18, 6, 14, 11, 9, 5, 11, 15, 6, 11, 17, 12, 15, 8, 4)

mean(data)
var(data)

# plot for data
plot(data, type="l", main="Deaths from horse kicks in Prussian Army")
plot(density(data))
hist(data, freq=F, ylim=c(0,0.12), xlim=c(0,18),main="Deaths from horse kicks in Prussian Army")
lines(density(data, bw = "SJ"))

### Model 1 - Poisson ###

# mle estimator
loglike_pois = function(par, data){
  x = data
  lambda = par[1]
  sum(dpois(x, lambda, log = TRUE))
}

mle_pois = optim(par=c(0.3), loglike_pois, data=data, method="BFGS", hessian = TRUE, control = list(fnscale = -1))
mle_pois$par

#plot logLike:
p.seq <- seq(0, 50, 0.001)
plot(p.seq, sapply(p.seq, loglike_pois, data=data), type="l", cex.axis=1.6, cex.lab=1.4, ylab="Log-Likelihood", xlab=expression(paste(lambda)))

# like
like_pois = function(data,lambda){prod(dpois(data, lambda, log = FALSE))}
# plot Like
p.seq <- seq(0, 20, 0.001)
plot(p.seq, sapply(p.seq, like_pois,data=data), type="l", cex.axis=1.6, cex.lab=1.4, ylab="Likelihood", xlab=expression(paste(lambda)))

# plot for the mass function
plot(data, dpois(data, mle_pois$par), ylim=c(0,0.14), type="h", lwd=1, col="blue", cex.axis=1.6, cex.lab=1.4,
     main="", ylab="Probability", xlab="Number of occurrences" )
abline(h=0, col="green2")
N=10^4
points(data, dbinom(data, N, 9.8/N), pch=19, col="darkgreen")
legend(13,0.145,lty = c(1,NA), pch=c(NA,19), cex=1.4, bty = "n", lwd = 1, col=c("blue","darkgreen"), legend=c("Poisson","Binomial") )

# plot posterior, prior likelihood
colors = c("orange","yellow", "purple","red","blue","green")
a = c(0.5,1,9,10,8,9)
b = c(0.01,5,1,1,2,9)
curve(dgamma(x,shape=a[4] + sum(data),rate=b[4] + n),cex.axis=1.6, cex.lab=1.4,col="red",xlab=expression(paste(lambda)),ylab="Density",lwd=1,from=0,to=18)
curve(dgamma(x,shape=a[4],rate=b[4]),col="blue",lwd=1,add=TRUE)
curve(dgamma(x, shape=sum(data)+1,rate=n), add=TRUE,type="l", col="green", xlab="Number of occurrences",ylab="Probability")
legend(0,0.58,lty = c(1), cex=1.4, bty = "n", lwd = 1, col=c("green","blue","red"), legend=c("Likelihood","Prior","Posterior") )

# prior plots
colors = c("orange","yellow", "purple","red","blue","green")
a = c(0.5,1,9,10,8,9)
b = c(0.01,1,1,1,2,9)
range_pois = seq(0,18,0.001)
plot(range_pois, dgamma(range_pois,a[1] , b[1]), ylim=c(0,1.3), col=colors[1],
     type="l", ylab="Prior Density", xlab = expression(paste(lambda)), cex.axis=1.6, cex.lab=1.4, main="")
for (i in 1:length(a)){
  lines(range_pois, dgamma(range_pois, a[i] , b[i]), 
        type="l", ylab="",  col=colors[i], lty=1)
}

legend(10.7,1.3,lty = c(1,rep(1,length(a))), cex=1.4, bty = "n", lwd = 1, 
       legend=c("Gamma(0.5,0.01)","Gamma(1,1)","Gamma(9,1)","Gamma(10,1)","Gamma(8,2)","Gamma(9,9)")
       , col=colors)

# posterior plots
dev.off()
colors = c("orange","yellow", "purple","red","blue","green")
a = c(0.5,1,9,10,8,9)
b = c(0.01,1,1,1,2,9)
range_pois = seq(0,18,0.001)
plot(range_pois, dgamma(range_pois,a[1] + sum(data), b[1] + n), col=colors[1], xlim=c(4,16), ylim=c(0,0.8),
     type="l", ylab="Posterior density", xlab = expression(paste(lambda)), cex.axis=1.6, cex.lab=1.4, main="")
for (i in 2:length(a)){
  lines(range_pois, dgamma(range_pois, a[i] + sum(data), b[i] + n), 
        type="l", ylab="", col=colors[i], lty=1)
}
legend(11.3,0.8,lty = c(1,rep(1,length(a))), cex=1.4, bty = "n", lwd = 1, 
       legend=c("Gamma(0.5,0.01)","Gamma(1,1)","Gamma(9,1)","Gamma(10,1)","Gamma(8,2)","Gamma(9,9)"), 
       col=colors)
#abline(v=mle_pois$par)


# posterior mean 
a = c(0.5,1,9,10,8,9)
b = c(0.01,1,1,1,2,9)
mean_pois=rep(0,length(a))
mode_pois=rep(0,length(a))
var_pois=rep(0,length(a))
quant_pois=matrix(0,length(a),3)

for (i in 1:length(a)){
  mean_pois[i] = (a[i] + sum(data))/ (b[i] + n)
  mode_pois[i] = (a[i] + sum(data) - 1)/ (b[i] + n)
  var_pois[i] = (a[i] + sum(data))/ (b[i] + n)^2
  quant_pois[i,] = qgamma(c(0.025, .5, 0.975), a[i] + sum(data), b[i] + n)
}

mean_pois
mode_pois
quant_pois


### Model 2 - Binomial ###

#large number of trials
N = 10^4

# mle estimator
loglike_bino <- function(par, data){
  x = data
  prob = par[1]
  sum(dbinom(x, N, prob, log = TRUE))
}

mle_bino = optim(par=c(0.01), loglike_bino, data=data, method="Brent", lower = 0, upper = 1, control = list(fnscale = -1))

mle_bino$par
loglike_bino(c(mle_bino$par),data=data)

dev.off()
# plot for the mass function
plot(data, dbinom(data, N, mle_bino$par), type="h", lwd=1, col="blue", ylim=c(0,0.14),cex.axis=1.6, cex.lab=1.4,
     main="", ylab="Probability", xlab="Possible Values")
abline(h=0, col="green2")
points(data, dbinom(data, N, 9.8/N), pch=19, col="darkgreen")

legend(13,0.145,lty = c(1,NA), pch=c(NA,19), cex=1.4, bty = "n", lwd = 1, col=c("blue","darkgreen"), legend=c("Poisson","Binomial") )

#plot logLike:
p.seq = seq(0.01, 0.99, 0.001)
plot(p.seq, sapply(p.seq, loglike_bino, data=data), type="l", cex.axis=1.6, cex.lab=1.4, ylab="Log-Likelihood",xlab="Probability Sucess")
#optimum:
optimize(loglike_bino, data=data,lower=0, upper=1, maximum=TRUE)

# like
like_bino = function(data,prob,N){prod(dbinom(data, N, prob, log = FALSE))}
# plot Like
p.seq = seq(0.00001,0.002,10^-5)
plot(p.seq, sapply(p.seq, like_bino, data=data, N=N), type="l", cex.axis=1.6, cex.lab=1.4, ylab="Likelihood", xlab="Probability Sucess")


# likelihood posterior prior
colors = c("black","green","orange","purple","red","yellow","blue")
alpha=c(0.5,1,3,2,5,10,60)
beta=c(0.5,1,3,4,24,10,40)
curve(dbeta(x, alpha[2] + sum(data), N*n - sum(data) + beta[2]),cex.axis=1.6, cex.lab=1.4,col="red",xlab="prob",ylab="",lwd=1,from=0,to=0.002)
curve(dbeta(x, alpha[2], beta[2]),col="blue",lwd=1,add=TRUE)
curve(dbeta(x, sum(data), N*n - sum(data)), add=TRUE,type="l", col="green", xlab="Number of occurrences",ylab="Probability")
legend(0,5700,lty = c(1), cex=1.4, bty = "n", lwd = 1, col=c("green","blue","red"), legend=c("Likelihood","Prior","Posterior") )


# prior plots
colors = c("black","green","orange","purple","red","yellow","blue")
alpha=c(0.5,1,3,2,5,10,60)
beta=c(0.5,1,3,4,24,10,40)
range_prior = seq(0.0001,1,10^-2)

plot(range_prior, dbeta(range_prior, alpha[1], beta[1]), col=colors[1],ylim=c(0,12), xlim=c(0,1),
     type="l", ylab="", xlab = "prob", cex.axis=1.6, cex.lab=1.4, main="")
for (i in 2:length(alpha)){
  lines(range_prior, dbeta(range_prior, alpha[i], beta[i]), 
        type="l", ylab="",  col=colors[i], lty=1)
}
legend(0.65,12.9,lty = c(1,rep(1,length(a))), cex=1.4, bty = "n", lwd = 1, 
       legend=c("Beta(0.5,0.5)","Beta(1,1)","Beta(3,3)","Beta(2,4)","Beta(5,24)","Beta(10,10)","Beta(60,40)")
       , col=colors)

# posterior plots
dev.off()
colors = c("black","green","orange","purple","red","yellow","blue" )
alpha=c(0.5,1,3,2,5,10,60)
beta=c(0.5,1,3,4,24,10,40)
range_post = seq(0.00001,0.003,10^-5)
plot(range_post, dbeta(range_post, alpha[1] + sum(data), N*n - sum(data) + beta[1]), col=colors[1],
     type="l", ylab="", xlab = "prob", cex.axis=1.6, cex.lab=1.4, main="",xlim=c(0.0006,0.002))
for (i in 2:length(alpha)){
  lines(range_post, dbeta(range_post, alpha[i] + sum(data), N*n - sum(data) + beta[i]), 
        type="l", ylab="", col=colors[i], lty=1)
}
legend(0.0015,5900,lty = c(1,rep(1,length(a))), cex=1.4, bty = "n", lwd = 1, 
       legend=c("Beta(0.5,0.5)","Beta(1,1)","Beta(3,3)","Beta(2,4)","Beta(5,24)","Beta(10,10)","Beta(60,40)"), 
       col=colors)
#abline(v=mle_bino$par)


# posterior mean 
alpha=c(0.5,1,3,2,5,10,60)
beta=c(0.5,1,3,4,24,10,40)
mean_bino=rep(0,length(alpha))
mode_bino=rep(0,length(alpha))
quant_bino=matrix(0,length(alpha),3)

for (i in 1:length(alpha)){
  mean_bino[i] = (alpha[i] + sum(data))/ (alpha[i] + beta[i] + N*n)
  mode_bino[i] = (alpha[i] + sum(data) - 1)/ ( alpha[i] + beta[i] + N*n - 2)
  quant_bino[i,] = qbeta(c(0.025, .5, 0.975), alpha[i] + sum(data), N*n - sum(data) + beta[i])
}

mean_bino
mode_bino
quant_bino


# Part 3

# Information Criteria
a = c(0.5,1,9,10,8,9)
b = c(0.01,1,1,1,2,9)

alpha=c(0.1,1,3,2,5,10,40)
beta=c(0.1,1,2,4,14,20,60)

i=4;j=2

# data
N = 10^4
n = 20
data = c(3, 5, 7, 9, 10, 18, 6, 14, 11, 9, 5, 11, 15, 6, 11, 17, 12, 15, 8, 4)

ns = 10000

# likelihood
fm1  = function(lambda){sum(dpois(data,lambda,log=TRUE))}
fm2  = function(prob){sum(dbinom(data,N,prob,log=TRUE))}

# maximum likelihood estimation
lambda.mle = optimize(fm1,lower=0, upper=15, maximum=TRUE)$m
prob.mle = optimize(fm2,lower=0, upper=1, maximum=TRUE)$m

# AIC and BIC
k = 1
AIC1 = -2*fm1(lambda.mle)+k
AIC2 = -2*fm2(prob.mle)+k
BIC1 = -2*fm1(lambda.mle)+k*log(n)
BIC2 = -2*fm2(prob.mle)+k*log(n)
AICs = c(AIC1,AIC2)
BICs = c(BIC1,BIC2)

AICs
BICs

# posterior distribution
lambda.post = rgamma(ns, a[i] + sum(data), b[i] + n)
prob.post = rbeta(ns, alpha[j] + sum(data), N*n - sum(data) + beta[j])

# DIC (deviance information criteria)
means1 = mean(lambda.post)
means2 = mean(prob.post)
mu11   = fm1(means1)
mu22   = fm2(means2)

pred1 = NULL; pred2 = NULL
for(t in 1:ns){
  pred1[t] = fm1(lambda.post[t])
  pred2[t] = fm2(prob.post[t])
}

pD1 = 2*(mu11 - mean(pred1)) 
pD2 = 2*(mu22 - mean(pred2))

DIC1 = -2*mu11 + 2*pD1  
DIC2 = -2*mu22 + 2*pD2 

pDs = c(pD1,pD2)
DICs = c(DIC1,DIC2)

DICs

#Gelfand and Ghosh
a = c(0.5,1,9,10,8,9)
b = c(0.01,1,1,1,2,9)

alpha=c(0.1,1,3,2,5,10,40)
beta=c(0.1,1,2,4,14,20,60)

i=4;j=2
nsp = 10000

# posterior distribution
lambda.post = rgamma(nsp,a[i] + sum(data), b[i] + n)
prob.post = rbeta(nsp, alpha[i] + sum(data), N*n - sum(data) + beta[i])

# predicted values
pred_values_1 = matrix(0,nsp,n)
pred_values_2 = matrix(0,nsp,n)

for(i in 1:length(lambda.post))
{
  #obtain one prediction for death
  pred_values_1[i,]=rpois(rep(1,n),rep(lambda.post[i],n))
  pred_values_2[i,]=rbinom(rep(1,n),N,rep(prob.post[i],n))
}

#G term of gelfand and ghosh
G1 = sum((apply(pred_values_1, 2,mean)-data)^2)
P1 = sum(apply(pred_values_1,2,var))
D1G = G1 + P1

G2 = sum((apply(pred_values_2, 2,mean)-data)^2)
P2 = sum(apply(pred_values_2,2,var))
D2G = G2 + P2

GGs = c(D1G,D2G)
GGs

## ALL 4 criterion
cbind(AICs, BICs, DICs,GGs)

# Bayes Factor
rm(list=ls(all=TRUE))
set.seed(7)
require(rmutil)

a = c(0.5,1,9,10,8,9);b = c(0.01,1,1,1,2,9)
alpha=c(0.5,1,3,2,5,10,40);beta=c(0.5,1,2,4,14,20,60)
i=2;j=2;N=10^4;ns = 500000;n=20

# data
data = c(3, 5, 7, 9, 10, 18, 6, 14, 11, 9, 5, 11, 15, 6, 11, 17, 12, 15, 8, 4)

# density of the parameters
f1  = function(lambda){sum(dpois(data,lambda,log=TRUE))}
f2  = function(prob){sum(dbinom(data,N,prob,log=TRUE))}

# generate ns parameters from the prior
BFprior1 = rgamma(ns, a[i], b[i])
BFprior2 = rbeta(ns, alpha[j], beta[j])

# using the density with the prior informations
# density of the parameters
BFden1 = matrix(0, ns, length(data))
BFden2 = matrix(0, ns, length(data))
for (i in 1:length(data)){
  BFden1[,i] = dpois(data[i],BFprior1,log=TRUE)
  BFden2[,i] = dbinom(data[i],N,BFprior2,log=TRUE)
}

BF1 = apply(BFden1, 1, sum)
BF2 = apply(BFden2, 1, sum)

BF12 = mean(exp(BF1)) / mean(exp(BF2))
BF12
BF21 = mean(exp(BF2)) / mean(exp(BF1))
BF21

cbind(BF12,BF21)

# using marginals 
poissongammaden = function(data,n,a,b) {
  a*log(b) - ((a+sum(data))*(log(n+b))) - sum(lfactorial(data)) + (lgamma(a+sum(data))) - (lgamma(a))}

binomialbetaden = function(data,n,N,alpha,beta) {
  sum(log(choose(N, data))) + lgamma(alpha+beta) + lgamma(alpha+sum(data)) + lgamma(n*N - sum(data) +beta) -
    lgamma(alpha) - lgamma(beta) - lgamma(n*N + alpha +beta)}

PG = poissongammaden(data,n, a[2], b[2])
BB = binomialbetaden(data,n,N,alpha[2],beta[2])

exp(PG-BB)
exp(BB-PG)

#Using the integration of the marginal
poisson_gamma_pdf = function(theta, par)
  # theta: parameter to be integrate
  # par: a list consists of data, alpha and beta
{
  #lambda = par$lambda
  a = par$a
  b = par$b
  data = par$data
  sum(dpois(data, lambda = theta, log=TRUE)) * dgamma(theta, shape = a, rate = b, log=TRUE)
}

bino_beta_pdf = function(theta, par)
  # theta: parameter to be integrate
  # par: a list consists of data, N, alpha and beta
{
  N = par$N
  #p = par$p
  data = par$data
  alpha = par$alpha
  beta = par$beta
  sum(dbinom(data, size=N, prob=theta, log=TRUE)) * dbeta(theta, shape1=alpha, shape2=beta, log=T)
}

bayes.factor =function(a, b, alpha, beta, N, data)
{
  m1=integrate(poisson_gamma_pdf, 0, 100, par=list(a=a, b=b, data=data))
  m2=integrate(bino_beta_pdf, 0, 1, par=list(N=N, data = data, alpha=alpha, beta=beta))
  m1$value / m2$value
}

bayes.factor(1, 1, 1, 1, N, data)
