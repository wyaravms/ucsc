
rm(list=ls(all=TRUE))
#set.seed(3)
library(mvtnorm)

sst = read.csv("C:\\Users\\WYARAVMS\\Google Drive\\phd-ucsc\\spring-2018\\ams-207\\take-home-3\\data\\sst.csv")
attach(sst)
Type = factor(Type)

nd = N
n = nrow(sst)
y = temp
X = cbind(rep(1,n),lat,lon,Type)
ncov = ncol(X)

hist(y,axes=FALSE, freq=FALSE, breaks = 20,main="Sea Surface Temperature", cex.lab=1.4,cex.main = 1.9,cex.axis=1.7)
lines(density(y), col=4,main="Temperature", cex.axis=1.5, cex.lab=1.5)
axis(1, cex.axis=1.6)
axis(2, cex.axis=1.6)

# exploratory plots
pairs(sst[,c(-1,-5)], pch=16, panel=panel.smooth, cex.axis=1.5, cex.lab=1.5)

fit = lm(temp ~ lat + lon + Type, data=sst, x=TRUE, y=TRUE)

X = fit$x
ncov = ncol(X)
summary(fit)

# gibbs function

sst.gibbs = function(X, y, nd, n, ncov, nmc, n.bur, thin){
  
  # hyperparameters
  alpha = 3; beta = 3
  a = 2; b =2
  
  post.mu.i = array(NA, dim=c(nmc, n))
  post.sigma.i = array(NA, dim=c(nmc, n))
  post.beta = array(NA, dim=c(nmc, ncov))
  post.sigma = array(NA, nmc)
  post.tau = array(NA, nmc)
  
  # initial values
  post.tau[1] = 1.5
  post.mu.i[1,] = y
  post.sigma.i[1,] = 1.5
  post.sigma[1] = 1.5
  
  for (i in 2:nmc){
    
    m.beta = solve(t(X)%*%X)%*%t(X)%*%post.mu.i[i-1,]
    var.beta = post.tau[i-1]*solve(t(X)%*%X)
    
    post.beta[i,] = rmvnorm(1, m.beta, var.beta)
    
    a.tau = (n/2) - 1 
    #b.tau = sum((post.mu.i[i-1,] - X%*%post.beta[i,])^2)/2
    b.tau = (post.mu.i[i-1,] - (X%*%post.beta[i,]))
    
    #sigma.hat = sum((y-X%*%m.beta)^2)/(n-qr(X)$rank)
    
    #post.tau[i] = 1/rgamma(1, ((n-qr(X)$rank)/2), (n-qr(X)$rank)*sigma.hat/2)
    post.tau[i] = 1/rgamma(1, a.tau, (t(b.tau)%*%b.tau)/2)
    
    m.mui = ((y*post.tau[i]*nd) + (X%*%post.beta[i,]*post.sigma.i[i-1,]))/(post.tau[i]*nd + post.sigma.i[i-1,])
    var.mui = (post.sigma.i[i-1,]*post.tau[i])/((post.tau[i]*nd)+(post.sigma.i[i-1,]))
    
    post.mu.i[i,] = rnorm(n,as.vector(m.mui),sqrt(var.mui))
    
    a.sigmai = 1/2 + (alpha + 1)
    b.sigmai = ((((y - post.mu.i[i,])^2)*nd)/2) + alpha*post.sigma[i-1]
    
    post.sigma.i[i,] = 1/rgamma(n,a.sigmai,b.sigmai)
    
    a.sigma = n*(alpha + 1) + a
    b.sigma = (alpha*sum(1/(post.sigma.i[i,]))) + b
    
    post.sigma[i] = rgamma(1, a.sigma, b.sigma)
    
    cat(i, "/", nmc, "\r")
  }
  
  n.bur = n.bur; thin = thin
  
  post.beta = post.beta[seq(n.bur+1, nmc,by=thin),]
  post.tau = post.tau[seq(n.bur+1, nmc,by=thin)]
  post.mu.i = post.mu.i[seq(n.bur+1, nmc,by=thin),]
  post.sigma.i = post.sigma.i[seq(n.bur+1, nmc,by=thin),]
  post.sigma = post.sigma[seq(n.bur+1, nmc,by=thin)]
  
  list(post.beta=post.beta, post.tau=post.tau, post.mu.i=post.mu.i, post.sigma.i=post.sigma.i, post.sigma=post.sigma)
}

posts = sst.gibbs(X = X, y = y, nd=nd, n=n, ncov=ncov, nmc=41000, n.bur=1000, thin=4)

nmcb = length(posts$post.sigma)
post.beta = posts$post.beta
post.tau = posts$post.tau
post.mu.i = posts$post.mu.i
post.sigma.i = posts$post.sigma.i
post.sigma = posts$post.sigma

# plotting the posterior distributions
# the chains, densities and correlations plots
par(mfrow=c(2,2))
plot.ts(post.tau, cex.lab=1.5, cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(tau^2)))
plot.ts(post.sigma, cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma^2)))
acf(post.tau, cex.lab=1.5, cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(tau^2)))
acf(post.sigma, cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma^2)))

par(mfrow=c(2,3))
for(i in 1:3)
  plot.ts(post.beta[,i], cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(beta[i]),list(i=i)), ylab="")
for(i in 1:3)
  acf(post.beta[,i], cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(beta[i]),list(i=i)), ylab="")

par(mfrow=c(2,3))
for(i in 4:6)
  plot.ts(post.beta[,i],  cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(beta[i]),list(i=i)), ylab="")
for(i in 4:6)
  acf(post.beta[,i],  cex.lab=1.5, cex.main = 2,cex.axis=1.8,main=substitute(paste(beta[i]),list(i=i)), ylab="")

par(mfrow=c(3,3))
for(i in 10:18)
  plot.ts(post.mu.i[,i],  cex.lab=1.5, cex.main = 2,cex.axis=1.8,main=substitute(paste(mu[i]),list(i=i)), ylab="")

par(mfrow=c(3,3))
for(i in 10:18)
  acf(post.mu.i[,i], cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(mu[i]),list(i=i)), ylab="")

# densities
par(mfrow=c(3,2))
for(i in 1:6){
  plot(density(post.beta[,i]), main=substitute(paste(beta[i]),list(i=i)),cex.lab=1.5, cex.main = 2,cex.axis=1.8, ylab="")
  abline(v=mean(post.beta[,i]), lty=1, col=2)
  abline(v=quantile(post.beta[,i], probs=c(0.025, 0.975)), lty=2, col=2)
}

par(mfrow=c(2,1))
plot(density(post.tau), cex.lab=1.5, cex.main = 2,cex.axis=1.8, col=1, ylab="",main=expression(paste(tau^2)))
abline(v=mean(post.tau), lty=1, col=2)
abline(v=quantile(post.tau, probs=c(0.025, 0.975)), lty=2, col=2)
plot(density(post.sigma), cex.lab=1.5,cex.main = 2,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma^2)))
abline(v=mean(post.sigma), lty=1, col=2)
abline(v=quantile(post.sigma, probs=c(0.025, 0.975)), lty=2, col=2)

dev.off()

# plot for sigma
# summary of each posterior distribution by the 5th and 95th quantiles
post.sum.sigma=apply(post.sigma.i,2,quantile,c(.05,.95))
post.sum.mean.sigma=apply(post.sigma.i,2,mean)
par(mfrow=c(1,1))
ind=1:length(y)
matplot(rbind(ind,ind),sqrt(post.sum.sigma),type="l",cex.axis=1.5, cex.lab=1.5,lty=1,col=4,xlab="index",ylab=expression(paste(sigma[i])))
points(ind,sqrt(post.sum.mean.sigma),pch=20)

# plot for mu
# summary of each posterior distribution by the 5th and 95th quantiles
post.sum.mu=apply(post.mu.i,2,quantile,c(.05,.95))
post.sum.mean.mu=apply(post.mu.i,2,mean)
par(mfrow=c(1,1))
ind=1:length(y)
matplot(rbind(ind,ind),post.sum.mu,type="l",cex.axis=1.5, cex.lab=1.5,lty=1,col=4,xlab="index",ylab=expression(paste(mu[i])))
points(ind,post.sum.mean.mu,pch=20)

# Posterior predictions for each observed observation
pred.y = matrix(0, n, nmcb)
for (i in 1:nmcb)
  pred.y[,i] = rnorm(n, post.mu.i[i,], sqrt(post.sigma.i[i,]/N))

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
matplot(rbind(ind,ind),pred.sum,type="l",cex.axis=1.5, cex.lab=1.5,lty=1,col=4,xlab="index",ylab="Temperature")
points(ind,y,pch=20)
out=(y<pred.sum[1,])
#text(ind[out], y[out], label=y[out], pos = 4)

# plot of the replicated & test statistics

# predictive test statistics
par(mfrow=c(2,2))
pred.max = (apply(pred.y, 2,max))
hist(pred.max, xlab="Maximum", freq=FALSE, breaks=20, main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=max(y), lty=2,lwd=2, col=2)   
p.pred.max = length(which((pred.max>(max(y))) == TRUE))
p.pred.max = round(p.pred.max/nmcb,5)
legend("topright",pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.pred.max)),bty='n',col=c(2))

pred.min = (apply(pred.y, 2,min))
hist(pred.min, xlab="Minimum", freq=FALSE, breaks=20, main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=min(y), lty=2,lwd=2, col=2)   
p.pred.min = length(which((pred.min>(min(y))) == TRUE))
p.pred.min = round(p.pred.min/nmcb,5)
legend(6,0.5,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.pred.min)),bty='n',col=c(2))

pred.mean = (apply(pred.y, 2,mean))
hist(pred.mean, xlab="Mean", freq=FALSE, breaks=20, main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=mean(y), lty=2,lwd=2, col=2)   
p.pred.mean = length(which((pred.mean>(mean(y))) == TRUE))
p.pred.mean = round(p.pred.mean/nmcb,5)
legend(18.5,3.2,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.pred.mean)),bty='n',col=c(2))

pred.sd = (apply(pred.y, 2, sd))
hist(pred.sd, xlab="Standard Deviation", freq=FALSE, breaks=20, main="",cex.lab=1.4,cex.main = 3,cex.axis=1.8)
abline(v=sd(y), lty=2,lwd=2, col=2)   
p.pred.sd = length(which((pred.sd>(sd(y))) == TRUE))
p.pred.sd = round(p.pred.sd/nmcb,5)
legend(1.5,4.5,pch=NA, cex=1.6, box.lty = 0, legend=bquote("p =" ~ .(p.pred.sd)),bty='n',col=c(2))

# plot of grids and predictions for new observations
x.grid = seq(min(lon), max(lon), by=0.1)
y.grid = seq(min(lat), max(lat), by=0.1) 

lat.long = expand.grid(x.grid, y.grid)

device.m = matrix(data = c(0,1,0,0,0,0,1,0,0,0,0,1), nrow = 4, ncol =3)

data = merge(lat.long,device.m, by = NULL)
data = cbind(rep(1,nrow(lat.long)),data)
n.l.l = nrow(data)

# using only 5000 sample from the posterior distributions
n.pred = 5000
pred.sigma.i = matrix(NA, n.l.l, n.pred)
for (i in 1:n.l.l)
  pred.sigma.i[i,] = 1/rgamma(n.pred,alpha+1,alpha*post.sigma[1:n.pred])

pred.beta = cbind(post.beta[1:n.pred,1], post.beta[1:n.pred,3], post.beta[1:n.pred,2], post.beta[1:n.pred,4], post.beta[1:n.pred,5],post.beta[1:n.pred,6])

mean.mu.i = as.matrix(data)%*%t(as.matrix(pred.beta))

pred.mu.i = matrix(0, n.l.l, n.pred)
for (i in 1:n.l.l)
  pred.mu.i[i,] = rnorm(n.pred, mean.mu.i[i,] , sqrt(post.tau[1:n.pred]))

pred.y.new = matrix(0, n.l.l, n.pred)
for (i in 1:n.l.l)
  pred.y.new[i,] = rnorm(n.pred, pred.mu.i[i,] , sqrt(pred.sigma.i[i,]))

pred.y.new.mean = apply(pred.y.new, 1, mean)
pred.y.new.sd = apply(pred.y.new, 1, sd)

# building the new design matrix
l = list(as.matrix(pred.y.new.mean), as.matrix(data[,-1]))
data.new = do.call(cbind,l)
colnames(data.new) = c("y", "lon", "lat", "dev2", "dev3", "dev4" )
data.new = data.frame(data.new)

data.new.dev1 = data.new[data.new$dev2 == "0" & data.new$dev3 == "0" & data.new$dev4 == "0",]
data.new.dev2 = data.new[data.new$dev2 == "1" & data.new$dev3 == "0" & data.new$dev4 == "0",]
data.new.dev3 = data.new[data.new$dev2 == "0" & data.new$dev3 == "1" & data.new$dev4 == "0",]
data.new.dev4 = data.new[data.new$dev2 == "0" & data.new$dev3 == "0" & data.new$dev4 == "1",]

l.sd = list(as.matrix(pred.y.new.sd), as.matrix(data[,-1]))
data.new.sd = do.call(cbind,l.sd)
colnames(data.new.sd) = c("y", "lon", "lat", "dev1", "dev2", "dev3" )
data.new.sd = data.frame(data.new.sd)

data.new.dev1.sd = data.new.sd[data.new.sd$dev1 == "0" & data.new.sd$dev2 == "0" & data.new.sd$dev3 == "0",]
data.new.dev2.sd = data.new.sd[data.new.sd$dev1 == "1" & data.new.sd$dev2 == "0" & data.new.sd$dev3 == "0",]
data.new.dev3.sd = data.new.sd[data.new.sd$dev1 == "0" & data.new.sd$dev2 == "1" & data.new.sd$dev3 == "0",]
data.new.dev4.sd = data.new.sd[data.new.sd$dev1 == "0" & data.new.sd$dev2 == "0" & data.new.sd$dev3 == "1",]

# plotting the grids and maps
library(ggplot2)
library(ggmap)
library(MASS)


map = get_map(location = c(lon = 29,lat = 33), zoom = 6, source = "google")
ggmap(map)

# map of the observed data
ggmap(map) +
  geom_point(aes(x = lon, y = lat, group = Type,
                 colour = Type),  data = sst, size=3)+
  geom_point(shape=22) +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

# map of the whole set of predictions
ggmap(map) +
  geom_tile(data = data.new, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: bucket", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

#maps for the means of the predictions of the new observations of each device
ggmap(map) +
  geom_tile(data = data.new.dev1, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: bucket", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

ggmap(map) +
  geom_tile(data = data.new.dev2, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: d.buoy", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

ggmap(map) +
  geom_tile(data = data.new.dev3, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: eri", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

ggmap(map) +
  geom_tile(data = data.new.dev4, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: f.buoy", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

# maps for the sd
ggmap(map) +
  geom_tile(data = data.new.dev1.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: bucket", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

ggmap(map) +
  geom_tile(data = data.new.dev2.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: d.buoy", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

ggmap(map) +
  geom_tile(data = data.new.dev3.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: eri", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

ggmap(map) +
  geom_tile(data = data.new.dev4.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: f.buoy", size=16)+
  theme(legend.text = element_text(size = 16), legend.title = element_text(size=16))

# plot of the grids
ggplot(data = data.new.dev1, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev1, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: bucket")

ggplot(data = data.new.dev2, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev2, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: d.buoy")

ggplot(data = data.new.dev3, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev3, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: eri")

ggplot(data = data.new.dev4, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev4, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "blue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: f.buoy")

# maps of the sd
ggplot(data = data.new.dev1.sd, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev1.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "steelblue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: bucket")

ggplot(data = data.new.dev2.sd, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev1.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "steelblue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: d.buoy")

ggplot(data = data.new.dev3.sd, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev1.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "steelblue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: eri")

ggplot(data = data.new.dev4.sd, aes(lon, lat, fill = y)) +
  geom_tile(data = data.new.dev1.sd, aes(lon, lat, fill = y)) +
  scale_fill_gradient(low = "white",high = "steelblue") +
  labs(x = "Longitude", y = "Latitude", fill = "temperature", title = "Temperature Device type: f.buoy")


# fitting a null model (M2) to compare with the main model (M1)
fit.null = lm(temp ~ 1, data=sst, x=TRUE, y=TRUE)
X.null = fit.null$x
ncov.null = ncol(X.null)
summary(fit.null)

posts.null = sst.gibbs(X=X.null, y=y, nd=nd, n=n, ncov=ncov.null, nmc=41000, n.bur=1000, thin=4)

nmcb = length(posts.null$post.sigma)
post.beta.null = posts.null$post.beta
post.tau.null = posts.null$post.tau
post.mu.i.null = posts.null$post.mu.i
post.sigma.i.null = posts.null$post.sigma.i
post.sigma.null = posts.null$post.sigma

# plotting the posterior distributions
# the chains, densities and correlations plots
par(mfrow=c(2,2))
plot.ts(post.tau.null, cex.lab=1.5, cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(tau^2)))
plot.ts(post.sigma.null, cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma^2)))
acf(post.tau.null, cex.lab=1.5, cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(tau^2)))
acf(post.sigma.null, cex.lab=1.5,cex.main = 3,cex.axis=1.8, col=1, ylab="",main=expression(paste(sigma^2)))

par(mfrow=c(1,2))
plot.ts(post.beta.null, cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(beta[0]),list(i=0)), ylab="")
acf(post.beta.null, cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(beta[0]),list(i=0)), ylab="")

par(mfrow=c(3,3))
for(i in 10:18)
  plot.ts(post.mu.i.null[,i],  cex.lab=1.5, cex.main = 2,cex.axis=1.8,main=substitute(paste(mu[i]),list(i=i)), ylab="")

par(mfrow=c(3,3))
for(i in 10:18)
  acf(post.mu.i.null[,i], cex.lab=1.5, cex.main = 2,cex.axis=1.8, main=substitute(paste(mu[i]),list(i=i)), ylab="")


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

# DIC (deviance information criteria) from the slides
likelihood = function(y, mu.dic, sigma.dic){
  val = sum(dnorm(y, mu.dic, sigma.dic, log = TRUE))
  return(val)
}

hlp = NULL
for(t in 1:nmcb){
  hlp[t] = likelihood(y, post.mu.i[t,], sqrt(post.sigma.i[t,]/N))
}
mu.mean.dic = apply(post.mu.i, 2, mean)
sigma.mean.dic = apply(post.sigma.i, 2, mean)

lph = likelihood(y, mu.mean.dic, sqrt(sigma.mean.dic/N))
pdic=2*(lph-mean(hlp))
DIC.1=-2*lph+2*pdic

hlp.2 = NULL
for(t in 1:nmcb){
  hlp.2[t] = likelihood(y, post.mu.i.null[t,], sqrt(post.sigma.i.null[t,]/N))
}
mu.mean.dic.2 = apply(post.mu.i.null, 2, mean)
sigma.mean.dic.2 = apply(post.sigma.i.null, 2, mean)

lph.2 = likelihood(y, mu.mean.dic.2, sqrt(sigma.mean.dic.2/N))
pdic.2 = 2*(lph.2-mean(hlp.2))
DIC.2 = -2*lph.2+2*pdic.2

c(DIC.1,DIC.2)
