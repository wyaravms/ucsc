
# HW 3 - AMS 207

rm(list=ls(all=TRUE))
setwd("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\spring-2018\\ams-207\\hw3")
set.seed(7)

library(LearnBayes)
library(mvtnorm)

data(birthweight)
attach(birthweight)
gender = factor(gender)

#birthweight = birthweight[-19,]
n = nrow(birthweight)

par(mfrow=c(1,2))
boxplot(age ~ gender, ylab="Age", xlab="Gender",  cex.axis=1.5, cex.lab=1.5)
boxplot(weight ~ gender, ylab="Weight",xlab="Gender", cex.axis=1.5, cex.lab=1.5)
dev.off()

plot(density(weight), axes=FALSE, main="Weight", cex.axis=1.5, cex.lab=1.5)
axis(1, at=c(2000, 2500, 3000, 3500, 3900), cex.axis=1.6)
axis(2, cex.axis=1.6)

# scatter plots
ageM = split(age, gender)$"0"
weightM = split(weight, gender)$"0"
ageF = split(age, gender)$"1"
weightF = split(weight, gender)$"1"
varM = cbind(weightM, ageM)
varF = cbind(weightF, ageF)
pairs(varM, pch=16, panel=panel.smooth, labels=c("Weight.M","Age.M"), cex.axis=1.5, cex.lab=1.5)
pairs(varF, pch=16, panel=panel.smooth, labels=c("Weight.F","Age.F"), cex.axis=1.5, cex.lab=1.5)

# fitting the model using lm
fit = lm(birthweight$weight ~ birthweight$age + birthweight$gender, qr=TRUE, x=TRUE, y=TRUE)
ls(fit)

nlm.coef = function(fit){
  X = model.matrix(fit)
  y = fit$y
  
  #rank 
  rank = qr(fit$x)$rank
  
  # extracting the decomposition Q and R
  Q = qr.Q(qr(X))
  R = qr.R(qr(X))
  
  # computing beta.hat best way
  beta.qr = backsolve(R,(t(Q)%*%y))
  
  # residual sums of squares
  rss = t(y-X%*%beta.qr)%*%(y-X%*%beta.qr)
  # estimate of residuals
  sigma = rss/(nrow(X) - qr(X)$rank) 
  
  # inverse of R
  Rinv = backsolve(r = R, x = diag(ncol(R)))
  
  # identity matrix to confirm
  Rind = round(R%*%Rinv, 3)
  
  # full variance-covariance solved system
  vcov = Rinv%*%t(Rinv)*as.vector(sigma)^2
  
  list(X=X, y=y, beta.hat=beta.qr, sigma=sigma, Rinv=Rinv, Rind=Rind, vcov=vcov, rank=rank)
}

# parameters used in the posterior
beta.hat = nlm.coef(fit)
X = beta.hat$X
y = beta.hat$y
p = beta.hat$rank
beta.mean = beta.hat$beta.hat
n.beta = nrow(beta.mean)
R = beta.hat$Rinv
sigma.hat = beta.hat$sigma

nmc = 10000

post.sigma = rep(NA, nmc)
post.beta = array(NA, dim=c(nmc, n.beta))

post.sigma = sqrt(1/rgamma(nmc, ((n-p)/2), (n-p)/2*sigma.hat))

# slow
for(i in 1:nmc)
  post.beta[i,] = rmvnorm(1, mean = t(beta.mean), sigma = tcrossprod(R)*(post.sigma[i])^2)

# much faster
post.beta = rmnorm(nmc, mean = rep(0,p), varcov = tcrossprod(R))
post.beta = array(1, c(nmc, 1))%*%t(beta.mean) + array(post.sigma, c(nmc, p))*post.beta

# another way
R = qr.R(qr(X))
z = rmnorm(nmc,rep(0,n.beta),diag(ncol(R))) 

beta = backsolve(R,t(z))

post.beta = array(NA, dim=c(nmc, n.beta))
post.beta = post.sigma*t(beta) + array(1, c(nmc, 1))%*%t(beta.mean)

apply(post.beta,2,mean)
apply(post.beta,2,sd)
apply(post.beta,2,quantile, probs=c(0.025,0.5, 0.975))
mean(post.sigma)
sd(post.sigma)
quantile(post.sigma, probs=c(0.025, 0.5,0.975))

write.table(post.beta,file="beta.txt")
write.table(post.sigma,file="sigma.txt")

# densities
par(mfrow=c(2,2))
for(i in 1:3){
  j = i-1
  plot(density(post.beta[,i]), axes=FALSE, cex.axis=1.5, cex.lab=1.5, main=substitute(paste(beta[j]),list(j=j)))
  axis(1, cex.axis=1.5)
  axis(2, cex.axis=1.5)
  abline(v=quantile(post.beta[,i], probs=c(0.025, 0.975)), lty=2, col=4)
  abline(v=mean(post.beta[,i]), lty=2, col=4)
}
plot(density(post.sigma), axes=FALSE,  cex.axis=1.5, cex.lab=1.5,main=substitute(paste(sigma)))
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
abline(v=quantile(post.sigma, probs=c(0.025, 0.975)), lty=2, col=4)
abline(v=mean(post.sigma), lty=2, col=4)
dev.off()

# histograms
par(mfrow=c(2,2))
for(i in 1:3){
  j = i-1
  hist(post.beta[,i], axes=FALSE, cex.axis=1.5, cex.lab=1.5,xlab=substitute(paste(beta[j]),list(j=j)), main="", freq=FALSE)
  axis(1, cex.axis=1.5)
  axis(2, cex.axis=1.5)
  abline(v=quantile(post.beta[,i], probs=c(0.025, 0.975)), lty=2, col=4)
  abline(v=mean(post.beta[,i]), lty=2, col=4)
}
hist(post.sigma, axes=FALSE, cex.axis=1.5, cex.lab=1.5,xlab=substitute(paste(sigma)), main="",freq=FALSE)
axis(1, cex.axis=1.5)
axis(2, cex.axis=1.5)
abline(v=quantile(post.sigma, probs=c(0.025, 0.975)), lty=2, col=4)
abline(v=mean(post.sigma), lty=2, col=4)
dev.off()

# Posterior predictions
pred.y = matrix(0, n, nmc)
for (i in 1:nmc)
  pred.y[,i] = rnorm(n, X%*% post.beta[i,], 1*post.sigma[i])

# summary of observed versus predicted
pred.mean = apply(pred.y, 1, mean)
pred.ci=apply(pred.y,1,quantile,c(.05,.95))
plot(y, pred.mean, ylim = range(pred.ci), cex.axis=1.5, cex.lab=1.5,pch = 20, ylab="Predicted", xlab="Observed")
abline(0, 1)
segments(x0 = y, y0 = pred.ci[1,], x1 = y, y1 = pred.ci[2,], col = 'blue')
points(y, pred.mean, pch = 20, col = 'black')

# summary of each predictive distribution by the 5th and 95th quantiles
pred.sum=apply(pred.y,1,quantile,c(.05,.95))
par(mfrow=c(1,1))
ind=1:length(y)
matplot(rbind(ind,ind),pred.sum,type="l",cex.axis=1.5, cex.lab=1.5,lty=1,col=4,xlab="index",ylab="Weight")
points(ind,y,pch=20)
out=(y>pred.sum[2,])
text(ind[out], y[out], label=y[out], pos = 4)

# possible outlier
hist(pred.y[19,], freq=FALSE, main="",cex.axis=1.5, cex.lab=1.5, xlab=expression(paste('y'[19])))
abline(v=y[19], col=4)
pv.y19 = length(which((pred.y[19,]>(y[19])) == TRUE))
pv.y19 = round(pv.y19/nmc,4)
legend("topleft",pch=NA, cex=1.5, box.lty = 0, legend=bquote("p =" ~ .(pv.y19)),bty='n',col=c(2))

# 36 and 40 weeks combined with Female and Male
cov1 = c(1,36,0)
cov2 = c(1,36,1)
cov3 = c(1,40,0)
cov4 = c(1,40,1)
X1 = rbind(cov1, cov2, cov3, cov4)

# expected response for the covariates
pred.cov = array(0, c(nmc, nrow(X1)))
for (i in 1:nrow(X1))
  pred.cov[,i] = t(X1[i,] %*% t(post.beta))

plot(density(pred.cov[,1]),xlim=c(2200,max(density(pred.cov[,1])$x)), 
     ylim=c(min(density(pred.cov[,1])$y),max(density(pred.cov[,3])$y)),axes=FALSE,
     cex.lab=1.5,cex.main = 1.5, col="red",main="")
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5) 
lines(density(pred.cov[,2]),cex.lab=1.5,cex.main = 1.5, col="blue",main="Female")
legend("topleft",lwd=1, lty=c(1,1), cex=1.3, box.lty = 0, 
       legend=c(expression('F-36 weeks'), expression('M-36 weeks')),bty='n',col=c(4,2))

plot(density(pred.cov[,3]),xlim=c(2500,max(density(pred.cov[,3])$x)),axes=FALSE,cex.lab=1.5,cex.main = 1.5, col="red", lty=1,main="")
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5) 
lines(density(pred.cov[,4]),cex.lab=1.5,cex.main = 1.5, col="blue",lty=1,main="Female")
legend("topleft",lwd=1, lty=c(1,1), cex=1.3, box.lty = 0, 
       legend=c(expression('F-40 weeks'), expression('M-40 weeks')),bty='n',col=c(4,2))

apply(pred.cov,2,quantile, probs=c(0.025,0.5, 0.975))

# Posterior predictions
pred.cov = array(0, c(nmc, nrow(X1)))
for (i in 1:nrow(X1))
  pred.cov[,i] = t(X1[i,] %*% t(post.beta)) + rnorm(nmc)*post.sigma 

plot(density(pred.cov[,1]),xlim=c(500,max(density(pred.cov[,4])$x)), 
     ylim=c(min(density(pred.cov[,1])$y),max(density(pred.cov[,3])$y)),axes=FALSE,
     cex.lab=1.5,cex.axis=1.5,cex.main = 1.5, col="red",main="")
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5) 
lines(density(pred.cov[,2]),cex.lab=1.5,cex.main = 1.5, col="blue",main="Female")

lines(density(pred.cov[,3]),cex.lab=1.5,cex.main = 1.5, col="red", lty=2,main="40 weeks")
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5) 
lines(density(pred.cov[,4]),cex.lab=1.5,cex.main = 1.5, col="blue",lty=2,main="Female")
legend("topleft",lwd=1, lty=c(1,1,2,2), cex=1.5, box.lty = 0, legend=c(expression('F-36 weeks'), expression('M-36 weeks'),expression('F-40 weeks'), expression('M-40 weeks')),bty='n',col=c(4,2,4,2))

apply(pred.cov,2,quantile, probs=c(0.025,0.5, 0.975))

# using the distribution
pred.cov = matrix(0, nmc, 4)
for (i in 1:4)
  pred.cov[,i] = rnorm(nmc, X1[i,]%*%t(post.beta), 1*post.sigma)

# histograms
par(mfrow=c(2,2))
hist(pred.cov[,1],cex.lab=1.5,cex.axis=1.5,cex.main = 1.5, freq=FALSE, col="#92C5DE", main="36 weeks - Female", xlab = "Weight")
hist(pred.cov[,2],cex.lab=1.5,cex.axis=1.5,cex.main = 1.5, freq=FALSE, col="#92C5DE", main="36 weeks - Male", xlab = "Weight")
hist(pred.cov[,3],cex.lab=1.5,cex.axis=1.5,cex.main = 1.5, freq=FALSE, col="#92C5DE", main="40 weeks - Female", xlab = "Weight")
hist(pred.cov[,4],cex.lab=1.5,cex.axis=1.5,cex.main = 1.5, freq=FALSE, col="#92C5DE", main="40 weeks - Male", xlab = "Weight")
dev.off()

# residuals

# hat matrix
hat.matrix = function(fit.qr) {
  # Q factor
  Q = qr.qy(fit.qr, diag(1, nrow = nrow(fit.qr$qr), ncol = fit.qr$rank))
  # QQ'
  tcrossprod(Q)
}

hat.m = hat.matrix(fit$qr)

outlier.prob = function (fit,beta.hat, hat.m, post.sigma, k) 
{
  e.hat = as.vector(beta.hat$y-beta.hat$X%*%beta.hat$beta.hat)
  h = diag(hat.m)
  prob = 0 * e.hat
  for (i in 1:length(prob)) {
    z1 = (k - e.hat[i]/sqrt(post.sigma))/sqrt(h[i])
    z2 = (-k - e.hat[i]/sqrt(post.sigma))/sqrt(h[i])
    prob[i] = mean(1 - pnorm(z1) + pnorm(z2))
  }
  list(prob=prob, e.hat=e.hat, h=h)
}

res.sum = outlier.prob(fit, beta.hat, hat.m, post.sigma, 3)
prob.outlier = res.sum$prob

par(mfrow=c(1,1))
plot(y,prob.outlier, ylab="Probability of outliers", cex.axis=1.5, cex.lab=1.5,cex.main = 1.5)
out = (prob.outlier > 2*pnorm(-3))
text(y[out], prob.outlier[out], label=which(out == TRUE), pos = 4)

# fitted values
pred.fitted = array(0, c(nmc, nrow(X)))
for (i in 1:nrow(X))
  pred.fitted[,i] = t(X[i,] %*% t(post.beta))

# residuals against the fitted values  
pred.mean.fitted = apply(pred.fitted, 2, mean)
plot(pred.mean.fitted,res.sum$e.hat, ylim=c(-300,510), cex.axis=1.5, cex.lab=1.5,cex.main = 1.5,ylab="Residuals", xlab="Fitted Values")
abline(h=0, lty=2, col=4)
identify(pred.mean.fitted,res.sum$e.hat)

r = res.sum$e.hat/(as.numeric(sqrt(beta.hat$sigma))*sqrt(1-res.sum$h))
#r = res.sum$e.hat/(as.numeric(sqrt(beta.hat$sigma)))
plot(pred.mean.fitted,r, ylim=c(-2,3), cex.axis=1.4, cex.lab=1.4,cex.main = 1.5,ylab="Standardized Residuals", xlab="Fitted Values")
abline(h=c(-2,0,2), lty=2, col=4)
identify(pred.mean.fitted,r)

e.bayes = array(0,c(nmc,n))
for(i in 1:nmc){
  e.bayes[i,] = (beta.hat$y-beta.hat$X%*%post.beta[i,])
}

cor.pred = array(0,c(nmc))
for(j in 1:nmc){
  cor.pred[j] = cor(e.bayes[j,],pred.fitted[j,])
}

plot(density(cor.pred), cex.axis=1.5, cex.lab=1.5,cex.main = 1.4, main = "Correlation")
abline(v=0, lty=2, col=4)

# correlation between the vectors of residuals
i = seq(1,24,2)
j = seq(2,24,2)
cor(res.sum$e.hat[i],res.sum$e.hat[j])

