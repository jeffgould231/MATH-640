## Example I.8 ##
X		<- 27
n		<- 100
curve(dbeta(x, X + 1, n - X + 1), from = 0, to = 1, xlab = expression(theta),
      ylab = expression(paste('p(',theta ,'|X)')), main = 'Posterior Density',
      lwd = 2, col = 'blue')

meanTheta	<- (X+1)/(n+2) # empirical Bayesian estimator
meanTheta

varTheta	<- ((X+1)*(n-X+1))/((n+3)*(n+2)^2)
varTheta
sqrt(varTheta)

modeTheta	<- X/n
modeTheta



## Example I.R3 ##
#install.packages('mcmcplots')
library(mcmcplots)
set.seed(526)
betaSamples	<- matrix(rbeta(10000, X + 1, n - X + 1), ncol = 1)
colnames(betaSamples)	<- 'theta'

denplot(betaSamples, main = "Sampled Posterior Density")

plot(density(betaSamples), xlab = expression(theta),
     ylab = expression(paste("p(",theta,"|X)")),
     main = "Sampled Posterior Density\nwith Credible Interval",
     lwd = 2, col = "blue", xlim = c(0,0.5))

credInt		<- quantile(betaSamples, probs = c(0.025, 0.975))

densOut		<- density(betaSamples)
idStart		<- max(which(densOut$x < credInt[1])) + 1
idEnd		<- min(which(densOut$x > credInt[2])) - 1
gx			<- densOut$x[idStart:idEnd]
gy			<- densOut$y[idStart:idEnd]
px			<- rep(0, length(gy))
polygon(c(gx, rev(gx)), c(px, rev(gy)), border = FALSE, col = rgb(0, 0, 1, alpha = 0.5))
abline(h = 0)