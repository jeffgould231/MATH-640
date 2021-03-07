## Example I.R5 ##
# install.packages('MCMCpack')
library(MCMCpack)

setwd('/Users/mjm556/Dropbox/Courses/MATH 640/Spring 2021/Note Sets/Part I/Data')
hers <- read.table('hers.txt', header = TRUE)
head(hers)

y		<- hers$x
ybar	<- mean(y)
sig		<- var(y)
n		<- length(y)
B		<- 10000

set.seed(3)

rmu		<- ybar + (sqrt(sig/n))*rt(B, df = n-1)
rsig	<- rinvgamma(B, (n-1)/2, ((n-1)/2)*sig)

median(rmu)
median(rsig)

quantile(rmu, probs = c(0.025, 0.975))
quantile(rsig, probs = c(0.025, 0.975))

par(mfrow = c(1,2))
plot(density(rmu), xlab = expression(mu), col = 'blue', lwd = 2, main = 'Posterior Mean')
plot(density(rsig), xlab = expression(sigma^2), col = 'blue', lwd = 2, main = 'Posterior Variance')



## Example I.R6 ##
set.seed(3)
B		<- 10000
rsig.R6	<- rinvgamma(B, (n-1)/2, ((n-1)/2)*sig)
rmu.R6	<- rnorm(B, ybar, sqrt(rsig.R6/n))

median(rmu.R6)
median(rsig.R6)

quantile(rmu.R6, probs = c(0.025, 0.975))
quantile(rsig.R6, probs = c(0.025, 0.975))

par(mfrow = c(1,2))
plot(density(rmu.R6), xlab = expression(mu), col = 'blue', lwd = 2)
plot(density(rsig.R6), xlab = expression(sigma^2), col = 'blue', lwd = 2)



## Example I.R7 ##
# library(MCMCpack)
set.seed(2584)
B		<- 10000
rtheta 	<- rdirichlet(B, c(338 + 1, 296 + 1, 282 + 1, 493 + 1))

rtMat	<- apply(rtheta, 2, quantile, probs = c(0.5, 0.025, 0.975))
colnames(rtMat)	<- c('Nat. Front', 'Republicans', 'En Marche!', 'Other')
round(t(rtMat), 4)

plot(density(rtMat[,1]), col = 'blue', lwd = 2, main = '',
     ylab = expression(paste('p(',theta,'|',y,')')),
     xlab = expression(theta), xlim = c(0.12,0.32), ylim = c(0,20))
lines(density(rtMat[,2]), col = 'red', lwd = 2, lty = 2)
lines(density(rtMat[,3]), col = 'green', lwd = 2, lty = 3)

lines(density(rtMat[,4]), col = 'gold', lwd = 2, lty = 4)

mean(rtheta[,3] >= 0.2401)
