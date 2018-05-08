
rm(list=ls(all=TRUE))
library(LearnBayes)
library(mvtnorm)

# pump failured
n = 8

# sufficient statistic
y1 = 23721
sumy = 15962989

# transformed parameters ~ 
log(sumy)
log(y1)

logpost = function(theta,postpar) {
  theta1 = theta[1]
  theta2 = theta[2]
  sumy = postpar$sumy
  y1 = postpar$y1
  n = postpar$n
  logpostr = -n*theta1 - ((1/exp(theta1)*(sumy - n*(y1 - exp(theta2))))) + theta1 + theta2
  return(logpostr)
}

postpar=list(sumy=sumy,y1=y1,n=n)

logf = function(theta, postpar) {
  if (is.matrix(theta) == TRUE) {
    val = matrix(0, c(dim(theta)[1], 1))
    for (j in 1:dim(theta)[1]) val[j] = logpost(theta[j,], postpar)
  }
  else val = logpost(theta, postpar)
  return(val)
}

limits = c(12,18,0,25)
ng = 100
x0 = seq(limits[1], limits[2], len = ng)
y0 = seq(limits[3], limits[4], len = ng)
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den = logf(cbind(X[1:n2], Y[1:n2]), postpar)
den = den - max(den)
den = matrix(den, c(ng, ng))
contour(x0, y0, den,levels = seq(-8.9, 0, by = 2.3), xlim=c(12.5,18),ylim=c(0,19),lwd = 1,
        xlab=expression(theta[1]),ylab=expression(theta[2]),cex.axis=1.2, cex.lab=1.2)

# estimating the mode of the posterior
logpost = function(theta,postpar) {
  theta1 = theta[1]
  theta2 = theta[2]
  sumy = postpar$sumy
  y1 = postpar$y1
  n = postpar$n
  logpostr = -n*theta1 - ((1/exp(theta1)*(sumy - n*(y1 - exp(theta2))))) + theta1 + theta2
  return(logpostr)
}

postpar=list(sumy=sumy,y1=y1,n=n)

start=c(25,25)
theta_mle = optim(par=start, logpost, postpar=postpar, method="BFGS",hessian = TRUE, control = list(fnscale = -1))
theta_mle$par
theta_mle$hessian

fisher = solve(-theta_mle$hessian)
fisher

mean1 = c(theta_mle$par[1],theta_mle$par[2])
Var1  = fisher
tpar=list(m=mean1,var=Var1,df=1)

logf_n = function(theta, tpar) {
  if (is.matrix(theta) == TRUE) {
    val = matrix(0, c(dim(theta)[1], 1))
    for (j in 1:dim(theta)[1]) val[j] = dmvnorm(theta[j,], mean = tpar$m, sigma = tpar$var,log=TRUE)
  } else val = dmvnorm(theta, mean = tpar$m, sigma = tpar$var,log=TRUE)
  return(val)
}

limits = c(12,18,0,25)
ng = 100
x0 = seq(limits[1], limits[2], len = ng)
y0 = seq(limits[3], limits[4], len = ng)
X = outer(x0, rep(1, ng))
Y = outer(rep(1, ng), y0)
n2 = ng^2
den1 = logf_n(cbind(X[1:n2], Y[1:n2]), tpar)
den1 = den1 - max(den1)
den1 = matrix(den1, c(ng, ng))
contour(x0, y0, den1, levels = seq(-8.9, 0, by = 2.3), xlim=c(12.5,18),ylim=c(0,19),lwd = 1,
        xlab=expression(theta[1]),ylab=expression(theta[2]),cex.axis=1.2, cex.lab=1.2, add=TRUE, col=2)

# reject sampling
log_rej = function(theta,postpar,tpar){
  l_rej = logpost(theta,postpar) - dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df,log=TRUE) 
  return(l_rej)
}

tpar = list(m=mean1,var=Var1,df=3)

start = c(15,15)
theta_rej = optim(par=start,log_rej,postpar=postpar,tpar=tpar,method="BFGS",hessian = TRUE, control = list(fnscale = -1))
theta_rej$par

dmax = log_rej(theta_rej$par,postpar=postpar,tpar=tpar)
dmax

# Accept/reject method
M = 15000
d = length(tpar$m)
theta = rmt(M,mean = c(tpar$m), S = tpar$var, df = tpar$df)
lf = matrix(0, c(dim(theta)[1], 1))

for (j in 1:dim(theta)[1]) lf[j] = logpost(theta[j,],postpar=postpar)

lg = dmt(theta,  mean = c(tpar$m), S = tpar$var, df = tpar$df, log = TRUE)
if (d == 1) {
  prob = exp(c(lf) - lg - dmax)
  draws1=(theta[runif(n) < prob])
} else {
  prob = exp(lf - lg - dmax)
  draws1=(theta[runif(n) < prob, ])
}

nd = nrow(draws1)
accept_rate = nd/M
accept_rate  

plot(draws1[,1],draws1[,2],
     xlab=expression(theta[1]),xlim=c(12.5,18),ylim=c(0,19),ylab=expression(theta[2]),pch=1)
contour(x0, y0, den, levels = seq(-8.9, 0, by = 2.3), xlim=c(12.5,18),ylim=c(0,19),lwd = 1,
        xlab=expression(theta[1]),ylab=expression(theta[2]),cex.axis=1.2, cex.lab=1.2, add=TRUE)

# SIR method
set.seed(7)
sir_sexp = function (logf, tpar, n, infor) 
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

theta = sir_sexp(logpost,tpar,nd,postpar)

plot(theta[,1],theta[,2],xlim=c(12.5,18),ylim=c(0,19),
     xlab=expression(theta[1]),ylab=expression(theta[2]),pch=1)
contour(x0, y0, den, levels = seq(-8.9, 0, by = 2.3), xlim=c(13,17),ylim=c(0,19),lwd = 1,
        xlab=expression(theta[1]),ylab=expression(theta[2]),cex.axis=1.2, cex.lab=1.2, add=TRUE)

plot(c(0,0,0),S$summary,type="b",lwd=3,xlim=c(-1,21),
        ylim=c(5,11), xlab="Observation removed",ylab="log K")
for (i in 1:20)
 lines(c(i,i,i),S$summary.obs[i,],type="b")

# density function
post = function(theta1,theta2,sumy,y1,n) {
  logpostr = -n*theta1 - ((1/exp(theta1)*(sumy - n*(y1 - exp(theta2))))) + theta1 + theta2
  return(exp(logpostr- max(logpostr)))} 


I = integrate(post,10,18,theta2=12,sumy=sumy,y1=y1,n=n)
par(mfrow=c(1,1))

#using normal distribution as proposal distribution
curve(post(x,theta2=12,sumy=sumy,y1=y1,n=n)/I$value,from=10,to=18,
      ylab="Density", xlab=expression(theta[1]),lwd=3, main="Densities",ylim=c(0,1))
curve(dnorm(x,mean1[1],sqrt(0.1666)),add=TRUE)

curve(dnorm(x,15,0.5),add=TRUE)
legend("topright",legend=c("Exact","Normal"),lwd=c(3,1))
curve(post(x,theta2=12,sumy=sumy,y1=y1,n=n)/I$value/dnorm(x,15,0.5),from=10,to=18, ylab="Weight",xlab="log trans",main="Weight = g/p")

# using t distribution as proposal
curve(post(x,theta2=12,sumy=sumy,y1=y1,n=n)/I$value,from=13,to=18,
      ylab="Density", xlab=expression(theta[1]),lwd=3, main="Densities",ylim=c(0,1))

curve(dmt(x,15,0.3,3),add=TRUE)
legend("topright",legend=c("Exact","Normal"),lwd=c(3,1))

curve(post(x,theta2=12,sumy=sumy,y1=y1,n=n)/I$value/dmt(x,15,0.5,2),from=10,to=18, ylab="Weight",xlab="log trans",main="Weight = g/p")

# importance sampling
post = function(theta,postpar) {
  sumy = postpar$sumy
  y1 = postpar$y1
  n = postpar$n
  theta1 = theta[1]
  theta2 = theta[2]
  logpostr = -n*theta1 - ((1/exp(theta1)*(sumy - n*(y1 - exp(theta2))))) + theta1 + theta2
  return(logpostr)} 

tpar = list(m=mean1,var=2*Var1,df=2)
postpar = list(sumy=sumy,y1=y1,n=n)

imporsampl = function (logf, tpar, h, n, data) 
{
  theta = rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
  lf = matrix(0, c(dim(theta)[1], 1))
  for (j in 1:dim(theta)[1]) lf[j] = logf(theta[j, ], data)
  H = lf
  for (j in 1:dim(theta)[1]) H[j] = h(theta[j, ])
  lp = dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
           log = TRUE)
  md = max(lf - lp)
  wt = exp(lf - lp - md)
  est = sum(wt * H)/sum(wt)
  SEest = sqrt(sum((H - est)^2 * wt^2))/sum(wt)
  return(list(est = est, se = SEest, theta = theta, wt = wt))
}

# mean
myfunc = function(theta) return(theta[2])

s = imporsampl(post,tpar,myfunc,10000,postpar)
cbind(s$est,s$se)
mu = s$est
mu

hist(s$wt,freq=FALSE)

# second moment
myfunc = function(theta) return(theta[1]^2)

s = imporsampl(post,tpar,myfunc,10000,postpar)
cbind(s$est,s$se)
mu2 = s$est
mu2

var_theta = mu2 - mu^2
var_theta

# variance
myfunc = function(theta) return((theta[2]-mu)^2)

s=imporsampl(post,tpar,myfunc,10000,postpar)
cbind(s$est,s$se)
var2 = s$est
var2

mycontour(post,c(10,20,0,25),postpar,tpar)
points(s$theta[,1],s$theta[,2])

# laplace approximation

myfunc = function(theta) return(log(theta))

post_laplace = function(theta,postpar,myfunc) {
  sumy = postpar$sumy
  y1 = postpar$y1
  n = postpar$n
  theta1 = theta[1]
  theta2 = theta[2]
  logpostr = -n*theta1 - ((1/exp(theta1)*(sumy - n*(y1 - exp(theta2))))) + theta1 + theta2 - log(myfunc(theta1)) - log(myfunc(theta2))
  return(logpostr)
} 

start=c(15,13)
theta_lap=optim(par=start,post_laplace,postpar=postpar, myfunc=myfunc,method="BFGS",hessian = TRUE, control = list(fnscale = -1))
theta_lap$par
theta_lap$hessian

mean1_lap = c(theta_lap$par[1],theta_lap$par[2])
mean1_lap

fisher_lap = solve(-theta_lap$hessian)
fisher_lap

mean1 = c(theta_mle$par[1],theta_mle$par[2])
Var1  = fisher

numer = (mean1_lap * sqrt(det(fisher_lap))*exp(-logpost(mean1_lap,postpar)))
denom = (sqrt(det(fisher))*exp(-logpost(mean1,postpar)))

expec_lap = numer/denom
expec_lap

library(psych)

m=outer(solve(fisher),1/mean1,FUN="+") 
A=solve(m[1,,])
A=sqrt(det(A))
B=sqrt(det(fisher))

# part 7
b = function(theta) return(exp(theta[,1]))
m = function(theta, smin) return(smin - exp(theta[,2]))

t0 = 10^6
R = function(t0,theta,smin) return(exp(-(t0 - m(theta,smin))/b(theta)))
hist(R(t0,theta,postpar$y1))
plot(density(R(t0,theta,postpar$y1)), main="")
abline(v=mean(R(t0,theta,postpar$y1)), col=2, lty=2)
abline(v=quantile(R(t0,theta,postpar$y1), probs = c(0.025, 0.975)), col=2, lty=2)

mean(R(t0,theta,postpar$y1))
var(R(t0,theta,postpar$y1))
sd(R(t0,theta,postpar$y1))
quantile(R(t0,theta,postpar$y1), probs = c(0.025, 0.975))
