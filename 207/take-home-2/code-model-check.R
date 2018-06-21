
mc = length(mu.m1)
n = length(volcano)

# Model 1
# posterior predictive distribution for each observation
pred.volcano.m1 = matrix(0, mc, n)
for (j in 1:n){
  pred.volcano.m1[,j] = rnorm(mc, mu.i.m1[j,], sqrt(sigma.m1/lambda.m1[j,]))
}

# posterior predictive for a new observation
pred.mu.m1 = rnorm(mc, mu.m1, sqrt(tau.m1))
pred.lambda.m1 = rgamma(mc, v/2, v/2)
pred.volcano.new.m1 = rnorm(mc, pred.mu.m1, sqrt(sigma.m1/pred.mu.m1))

# Model 2
# posterior predictive distribution for each observation
pred.volcano.m2 = matrix(0, mc, n)
for (j in 1:n){
  pred.volcano.m2[,j] = rnorm(mc, mu.m2, sqrt(sigma.m2))
}

# posterior predictive for a new observation
pred.volcano.new.m2 = rnorm(mc, mu.m2, sqrt(sigma.m2))

# graph histogram and the future data
hist(volcano, xlab="Log of interevent time",breaks=20, ylim=c(0,0.4),main="",freq=FALSE,cex.lab=1.4,cex.main = 1.9,cex.axis=1.7)
lines(density(pred.volcano.new.m1),col="red")
lines(density(pred.volcano.new.m2),col="blue")
legend(1,0.4,lwd=1, cex=1.6, box.lty = 0, legend=c("M1", "M2"),bty='n',col=c(2,4))

# graph of the predictive distribution
plot(density(pred.volcano.m1[,1]),ylim=c(0,0.5),xlim=c(-5,15),cex.lab=1.5,cex.main = 1.5,cex.axis=1.8, col="grey", xlab="",main="M1")
for(i in 2:length(volcano)){
  lines(density(pred.volcano.m1[,i]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col="grey", xlab="",main=expression(paste(lambda[i])))
}
lines(density(pred.volcano.m1[,9]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=2, xlab="",main=expression(paste(lambda[i])))
lines(density(pred.volcano.m1[,15]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=4, xlab="",main=expression(paste(lambda[i])))
#lines(density(volcano),col="red")
#legend("topleft",lwd=1, cex=1.6, box.lty = 0, legend=c(expression('log T'), expression('log T'[rep])),bty='n',col=c(2,"grey"))
legend("left",lwd=1, cex=1.1, box.lty = 0, legend=c(expression('log T'[9]),expression('log T'[15])),bty='n',col=c(2,4))


# graph of the predictive distribution
plot(density(pred.volcano.m2[,1]),ylim=c(0,0.43),xlim=c(-5,15),cex.lab=1.5,cex.main = 1.5,cex.axis=1.8, col="grey", xlab="",main="M2")
for(i in 2:length(volcano)){
  lines(density(pred.volcano.m2[,i]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col="grey", xlab="",main=expression(paste(lambda[i])))
}
lines(density(pred.volcano.m1[,9]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=2, xlab="",main=expression(paste(lambda[i])))
lines(density(pred.volcano.m1[,15]),cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=4, xlab="",main=expression(paste(lambda[i])))
#lines(density(volcano),col="red")
#legend("topleft",lwd=1, cex=1.6, box.lty = 0, legend=c(expression('log T'), expression('log T'[rep])),bty='n',col=c(2,"grey"))
legend("topleft",lwd=1, cex=1.1, box.lty = 0, legend=c(expression('log T'[9]),expression('log T'[15])),bty='n',col=c(2,4))


# DIC (deviance information criteria)
mean.lambda.m1 = apply(MCMC.sam.M1$lambda[,seq(n.bur+1, ns,by=thin)],1,mean)
mean.mu.i.m1 = apply(MCMC.sam.M1$mu.i[,seq(n.bur+1, ns,by=thin)],1,mean)
mean.mu.m1 = mean(MCMC.sam.M1$mu[seq(n.bur+1, ns,by=thin)])
mean.sigma.m1 = mean(MCMC.sam.M1$sigma[seq(n.bur+1, ns,by=thin)])
mean.tau.m1 = mean(MCMC.sam.M1$tau[seq(n.bur+1, ns,by=thin)])

mean.mu.m2 = mean(MCMC.sam.M2$mu[seq(n.bur+1, ns,by=thin)])
mean.sigma.m2 = mean(MCMC.sam.M2$sigma[seq(n.bur+1, ns,by=thin)])

# deviance statistics
mm1 = matrix(0, 1, length(volcano))
for (i in 1:length(volcano)){
  mm1[,i] = dnorm(volcano[i], mean.mu.i.m1[i], sqrt(mean.sigma.m1/mean.lambda.m1[i]), log = TRUE)
}

mm2 = matrix(0, 1, length(volcano))
for (i in 1:length(volcano)){
  mm2[,i] = dnorm(volcano[i], mean.mu.m2, sqrt(mean.sigma.m2), log = TRUE)
}

mm1 = sum(mm1)
mm2 = sum(mm2)

m1 = matrix(0, length(mu.m1), length(volcano))
for (i in 1:length(volcano)){
  m1[,i] = dnorm(volcano[i], mu.i.m1[i,], sqrt(sigma.m1/lambda.m1[i,]), log = TRUE)
}

m2 = matrix(0, length(mu.m2), length(volcano))
for (i in 1:length(volcano)){
  m2[,i] = dnorm(volcano[i], mu.m2, sqrt(sigma.m2), log = TRUE)
}

pD1 = 2*(mm1 - mean(m1)) 
pD2 = 2*(mm2 - mean(m2))

DIC1 = -2*mm1 + 2*pD1  
DIC2 = -2*mm2 + 2*pD2 

pDs = c(pD1,pD2)
DICs = c(DIC1,DIC2)

pDs
DICs

# another way
m1 = matrix(0, length(mu.m1))
for (i in 1:length(mu.m1)){
  m1[i] = -2*sum(dnorm(volcano, mu.i.m1[,i], sqrt(sigma.m1[i]/lambda.m1[,i]), log = TRUE))
}

m2 = matrix(0, length(mu.m2))
for (i in 1:length(mu.m2)){
  m2[i] = -2*sum(dnorm(volcano, mu.m2[i], sqrt(sigma.m2[i]), log = TRUE))
}

m1.dic = mean(m1) + var(m1)/2
m2.dic = mean(m2) + var(m2)/2

m1 = matrix(0, length(mu.m1), length(volcano))
for (i in 1:length(volcano)){
  m1[,i] = dnorm(volcano[i], mu.i.m1[i,], sqrt(sigma.m1/lambda.m1[i,]), log = TRUE)
}

m2 = matrix(0, length(mu.m2), length(volcano))
for (i in 1:length(volcano)){
  m2[,i] = dnorm(volcano[i], mu.m2, sqrt(sigma.m2), log = TRUE)
}

m1.2v = -2*apply(m1, 1, sum)
m1.DIC = mean(m1.2v) + var(m1.2v)/2

m2.2v = -2*apply(m2, 1, sum)
m2.DIC = mean(m2.2v) + var(m2.2v)/2

pDs.2v = c(var(m1.2v)/2, var(m2.2v)/2)
pDs.2v

DICs.2v = c(m1.DIC,m2.DIC)
DICs.2v


#Gelfand and Ghosh

#G term of gelfand and ghosh
G1 = sum((apply(pred.volcano.m1, 2,mean)-volcano)^2)
P1 = sum(apply(pred.volcano.m1,2,var))
D1G = G1 + P1

G2 = sum((apply(pred.volcano.m2, 2,mean)-volcano)^2)
P2 = sum(apply(pred.volcano.m2,2,var))
D2G = G2 + P2

GGs = c(D1G,D2G)
GGs

# Bayes factor approximation
M = 1000000 

bf.tau.m1 = 1/rgamma(M, par.sam$c, par.sam$d)
bf.mu.m1 = rnorm(M, par.sam$m, sqrt(par.sam$s))
prior.mean.m1 = rnorm(M, bf.mu.m1, sqrt(bf.tau.m1))

bf.sigma.m1 = 1/rgamma(M, par.sam$a, par.sam$b)
bf.lambda.m1 = matrix(0, M, length(volcano))
for (i in 1:length(volcano)) {
  bf.lambda.m1[,i] = rgamma(M, v/2, v/2)
}

prior.sd.m1 = matrix(0, M, length(volcano))
for (i in 1:length(volcano)) {
  prior.sd.m1[,i] = sqrt(bf.sigma.m1 /bf.lambda.m1[,i])
}

bf.m1 = matrix(0, M, length(volcano))

for (i in 1:length(volcano))
  bf.m1[,i] = dnorm(volcano[i], prior.mean.m1, prior.sd.m1[,i], log = TRUE)

bf.m1 = apply(bf.m1, 1, sum)

prior.mean.m2 = rnorm(M, par.sam$m, sqrt(par.sam$s))
prior.sd.m2 = sqrt(1/rgamma(M, par.sam$a, par.sam$b))

bf.m2 = matrix(0, M, length(volcano))
for (i in 1:length(volcano))
  bf.m2[,i] = dnorm(volcano[i], prior.mean.m2, prior.sd.m2, log = TRUE)
bf.m2 = apply(bf.m2, 1, sum)

BF12 = mean(exp(bf.m1)) / mean(exp(bf.m2))
BF12

BF21 = mean(exp(bf.m2)) / mean(exp(bf.m1))
BF21

# predictive first
# model 1
# posterior predictive distribution for each observation
prior.volcano.m1 = matrix(0, mc, n)
for (j in 1:n){
  prior.volcano.m1[,j] = rnorm(mc, prior.mean.m1[j], prior.sd.m1[j])
}

# Model 2
# posterior predictive distribution for each observation
prior.volcano.m2 = matrix(0, mc, n)
for (j in 1:n){
  prior.volcano.m2[,j] = rnorm(mc, prior.mean.m2[j], prior.sd.m2[j])
}

max.prior.m1 = (apply(prior.volcano.m1, 1,mean))
max.pred.m1 = (apply(pred.volcano.m1, 1,mean))
plot(max.prior.m1[1:200],max.pred.m1[1:200], pch=18, xlim=c(0,7),ylim=c(0,7))
abline(0,1)

plot(volcano, apply(pred.volcano.m1, 2, mean),ylim=c(0,7))

# model 1
# predictive test statistics
par(mfrow=c(2,2))
max.pred.m1 = (apply(pred.volcano.m1, 1,max))
hist(max.pred.m1, xlab="Maximum", main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=max(volcano), lty=2,lwd=2, col=2)   
p.max.pred.m1 = length(which((max.pred.m1>(max(volcano))) == TRUE))
p.max.pred.m1 = round(p.max.pred.m1/mc,3)
legend(13,4500,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.max.pred.m1)),bty='n',col=c(2))

min.pred.m1 = (apply(pred.volcano.m1, 1,min))
hist(min.pred.m1, xlab="Minimum",main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=min(volcano), lty=2, lwd=2,col=2)    
p.min.pred.m1 = length(which((min.pred.m1>(min(volcano))) == TRUE))
p.min.pred.m1 = round(p.min.pred.m1/mc,3)
legend(-20,6900,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.min.pred.m1)),bty='n',col=c(2))

mean.pred.m1 = (apply(pred.volcano.m1, 1,mean))
hist(mean.pred.m1, xlab="Mean", main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=mean(volcano), lty=2, lwd=2,col=2)   
p.mean.pred.m1 = length(which((mean.pred.m1>(mean(volcano))) == TRUE))
p.mean.pred.m1 = round(p.mean.pred.m1/mc,3)
legend(5.5,2200,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.mean.pred.m1)),bty='n',col=c(2))

#med.pred.m1 = (apply(pred.volcano.m1, 1,median))
#hist(med.pred.m1)
#abline(v=median(volcano), lty=2, col=2)               
#p.med.pred.m1 = length(which((med.pred.m1>(median(volcano))) == TRUE))
#p.med.pred.m1 = p.med.pred.m1/mc
#legend(5.7,1100,lty = 2, cex=1, box.lty = 0, lwd = 2, legend=bquote("p =" ~ .(p.med.pred.m1)),bty='n',col=c(2))

sd.pred.m1 = (apply(pred.volcano.m1, 1,sd))
hist(sd.pred.m1, xlab="Standard Deviation", main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=sd(volcano), lty=2,lwd=2, col=2)               
p.sd.pred.m1 = length(which((sd.pred.m1>(sd(volcano))) == TRUE))
p.sd.pred.m1 = round(p.sd.pred.m1/mc,3)
legend(2,5000,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.sd.pred.m1)),bty='n',col=c(2))


# model 2
# predictive test statistics
par(mfrow=c(2,2))
max.pred.m2 = (apply(pred.volcano.m2, 1,max))
hist(max.pred.m2, xlab="Maximum", main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=max(volcano), lty=2,lwd=2, col=2)   
p.max.pred.m2 = length(which((max.pred.m2>(max(volcano))) == TRUE))
p.max.pred.m2 = round(p.max.pred.m2/mc,3)
legend(8.7,3400,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.max.pred.m2)),bty='n',col=c(2))

min.pred.m2 = (apply(pred.volcano.m2, 1,min))
hist(min.pred.m2, xlab="Minimum",main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=min(volcano), lty=2, lwd=2,col=2)    
p.min.pred.m2 = length(which((min.pred.m2>(min(volcano))) == TRUE))
p.min.pred.m2 = round(p.min.pred.m2/mc,3)
legend(-2.5,3400,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.min.pred.m2)),bty='n',col=c(2))

mean.pred.m2 = (apply(pred.volcano.m2, 1,mean))
hist(mean.pred.m2, xlab="Mean", main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=mean(volcano), lty=2, lwd=2,col=2)   
p.mean.pred.m2 = length(which((mean.pred.m2>(mean(volcano))) == TRUE))
p.mean.pred.m2 = round(p.mean.pred.m2/mc,3)
legend(4.3,2000,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.mean.pred.m2)),bty='n',col=c(2))

#med.pred.m2 = (apply(pred.volcano.m2, 1,median))
#hist(med.pred.m2)
#abline(v=median(volcano), lty=2, col=2)               
#p.med.pred.m2 = length(which((med.pred.m2>(median(volcano))) == TRUE))
#p.med.pred.m2 = p.med.pred.m2/mc
#legend(5.7,1100,lty = 2, cex=1, box.lty = 0, lwd = 2, legend=bquote("p =" ~ .(p.med.pred.m2)),bty='n',col=c(2))

sd.pred.m2 = (apply(pred.volcano.m2, 1,sd))
hist(sd.pred.m2, xlab="Standard Deviation", main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=sd(volcano), lty=2,lwd=2, col=2)               
p.sd.pred.m2 = length(which((sd.pred.m2>(sd(volcano))) == TRUE))
p.sd.pred.m2 = round(p.sd.pred.m2/mc,3)
legend(1.3,3000,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.sd.pred.m2)),bty='n',col=c(2))
