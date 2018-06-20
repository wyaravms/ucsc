
rm(list=ls(all=TRUE))
set.seed(27)

# take home - final - AMS 206

y = c(82, 79, 74, 83, 80, 81, 84, 81)
x1 = c(10, 9, 9, 11, 11, 10, 10, 12)
x2 = c(15, 14, 13, 15, 14, 14, 16, 13)

X = cbind(x1, x2)

fit.1 = lm(y ~ -1 + x1 + x2) 
summary(fit.1)
p = 2
n = length(y)

XtX = t(X)%*%X
invXtX = solve(XtX)
Xty = t(X)%*%y

# Question 1

# 1.a

beta.hat = solve(t(X)%*%X)%*%t(X)%*%y
sigma.hat = sum((y-X%*%beta.hat)^2)/(n-p)

# 1.b
lambda.1 = c(1,0)
var.est.1 = sqrt((sigma.hat)%*%(t(lambda.1)%*%invXtX%*%lambda.1))

t(lambda.1)%*%beta.hat + qt(0.975,6)*var.est.1
t(lambda.1)%*%beta.hat - qt(0.975,6)*var.est.1

lambda.2 = c(2,1)
var.est.2 = sqrt((sigma.hat)%*%(t(lambda.2)%*%invXtX%*%lambda.2))

t(lambda.2)%*%beta.hat + qt(0.975,6)*var.est.2
t(lambda.2)%*%beta.hat - qt(0.975,6)*var.est.2

# 1.c
lambda.3 = c(0,1)
var.est.3 = sqrt((sigma.hat)%*%(t(lambda.3)%*%invXtX%*%lambda.3))

t(lambda.3)%*%beta.hat + qt(0.975,6)*var.est.3
t(lambda.3)%*%beta.hat - qt(0.975,6)*var.est.3

t.val.3 = (t(lambda.3)%*%beta.hat - 3) / var.est.3
pval.3 = (1-pt(abs(t.val.3),n-p))*2

# 1.d
lambda.4 = c(1,-1)
var.est.4 = sqrt((sigma.hat)%*%(t(lambda.4)%*%invXtX%*%lambda.4))

t(lambda.4)%*%beta.hat + qt(0.975,6)*var.est.4
t(lambda.4)%*%beta.hat - qt(0.975,6)*var.est.4

t.val.4 = (t(lambda.4)%*%beta.hat - 0) / var.est.4
pval.4 = (1-pt(abs(t.val.4),n-p))*2

# Question 2

# 2.a
fit.2 = lm(y ~ -1 + x1 + x2)
summary(fit.2)

lambda.5 = c(0,1)
var.est.5 = sqrt((sigma.hat^2)%*%(t(lambda.5)%*%invXtX%*%lambda.5))

t.val.5 = (t(lambda.5)%*%beta.hat - 0) / var.est.5
pval.5 = (1-pt(abs(t.val.5),n-p))*2

# 2.b
library(ellipse)
plot(ellipse(fit.2, which = c(1,2), level = 0.95), type = 'l', cex.lab=1.5,cex.main = 2,cex.axis=1.8)
points(fit.2$coefficients[1], fit.2$coefficients[2], col=2)
segments(fit.2$coefficients[1],0,fit.2$coefficients[1],fit.2$coefficients[2])
segments(0,fit.2$coefficients[2],fit.2$coefficients[1],fit.2$coefficients[2])
legend(3.5,4.5, pch=1, box.lty = 0,cex=1.6, legend=c(expression(paste(hat(beta)[1]:hat(beta)[2]))),bty='n',col=c(2))



beta1.hat <- beta.hat[1]
beta2.hat <- beta.hat[2]
se1 <- std.error[1]
se2 <- std.error[2]
alpha <- 0.05
beta1.lo <- beta1.hat + se1*qt(alpha/2, df=n-p)
beta1.hi <- beta1.hat - se1*qt(alpha/2, df=n-p)
beta2.lo <- beta2.hat + se2*qt(alpha/2, df=n-p)
beta2.hi <- beta2.hat - se2*qt(alpha/2, df=n-p)
beta1.lo.ex <- beta1.hat + se1*qt((alpha/2)/2, df=n-p)
beta1.hi.ex <- beta1.hat - se1*qt((alpha/2)/2, df=n-p)
beta2.lo.ex <- beta2.hat + se2*qt((alpha/2)/2, df=n-p)
beta2.hi.ex <- beta2.hat - se2*qt((alpha/2)/2, df=n-p)
plot(beta1.hat, beta2.hat, xlim=c(-1,5), ylim=c(1,6),
     xlab=expression(beta[1]), ylab=expression(beta[2]))
segments(x0=c(beta1.hat,beta1.lo), y0=c(beta2.lo,beta2.hat),
         x1=c(beta1.hat,beta1.hi), y1=c(beta2.hi,beta2.hat), col="grey")
polygon(x=c(beta1.lo.ex,beta1.lo.ex,beta1.hi.ex,beta1.hi.ex),
        y=c(beta2.lo.ex,beta2.hi.ex,beta2.hi.ex,beta2.lo.ex))

# 2.c

fit.3 = lm(y ~ x1 + x2)
summary(fit.3)

plot(ellipse(fit.3, which = c(1,2), level = 0.95), type = 'l')
plot(ellipse(fit.3, which = c(2,3), level = 0.95), type = 'l')


# Question 4

data = matrix(c(0.620, 1, 1,
  1.342, 1, 1,
  0.669, 1, 1,
  0.687, 1, 1,
  0.155, 1, 1,
  2.000, 1, 1,
  1.228, 1, 2,
  3.762, 1, 2,
  2.219, 1, 2,
  4.207, 1, 2,
  0.615, 1, 3,
  2.245, 1, 3,
  2.077, 1, 3,
  3.357, 1, 3,
  1.182, 2, 1,
  1.068, 2, 1,
  2.545, 2, 1,
  2.233, 2, 1, 
  2.664, 2, 1, 
  1.002, 2, 1, 
  2.506, 2, 1,
  4.285, 2, 1,
  1.696, 2, 1,
  3.080, 2, 2, 
  2.741, 2, 2,
  2.522, 2, 2,
  1.647, 2, 2,
  1.999, 2, 2,
  2.939, 2, 2,
  2.240, 2, 3,
  0.330, 2, 3,
  3.453, 2, 3,
  1.527, 2, 3,
  0.809, 2, 3,
  1.942, 2, 3), nrow = 35, ncol = 3, byrow = TRUE, dimnames=list(c(seq(1,35,1)),c("Scores","BG","Major")))

data = as.data.frame(data)  

fit.4 = lm(Scores ~ as.factor(Major)*as.factor(BG), data=data)
summary(fit.4)

anova(fit.4)

names <- c('Major','BG')
data[, names] = lapply(data[ , names], factor)

levels(data$Major)[c(1,2,3)] = c('Economics', 'Anthropology', "Sociology")
levels(data$BG)[c(1,2)] = c("Rural", "Urban")

interaction.plot(data$Major, data$BG, data$Scores, trace.lab="BG", legend=F,cex.axis=1.5, cex.lab=1.5, xlab="Major", ylab="Mean of Score", cex=1.5)
legend(2.3,1.6, lty = c(2,1), box.lty = 0,cex=1.5, legend=c("Rural","Urban"))

interaction.plot(data$BG, data$Major, data$Scores, trace.lab="Major",cex.axis=1.5, cex.lab=1.5, xlab="BG", ylab="Mean of Score", cex=1.5)

library(ggplot2)
ggplot(aes(y = Scores, x = BG, fill = Major), data = data) + 
  geom_boxplot()

fit.5 = lm(Scores ~ 1, data=data)
summary(fit.5)

anova(fit.4, fit.5)
