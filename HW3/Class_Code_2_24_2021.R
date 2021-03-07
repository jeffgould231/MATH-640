### Example I.R8 ###
library(mvtnorm)
library(MCMCpack)

hers  <- read.table('hersreg.txt', header = TRUE)

B     <- 10000
n     <- nrow(hers)
X     <- model.matrix(chtchol ~ treatment + sbp + dbp + statins, data = hers)
Y     <- hers$chtchol
K     <- ncol(X)

bhat	<- c(solve(t(X)%*%X)%*%(t(X)%*%Y))
SSY	<- t(Y - X%*%bhat)%*%(Y - X%*%bhat)
XtXi	<- solve(t(X)%*%X)
rbeta	<- matrix(0, nrow = B, ncol = K)

set.seed(90210)
rsig	<- rinvgamma(B, (n-K)/2, (1/2)*SSY)
for(i in 1:B){
  CovX		<- rsig[i]*XtXi
  rbeta[i,]	<- c(rmvnorm(1, mean = bhat, sigma = CovX))
}

rbMat	   <- apply(rbeta, 2, quantile, probs = c(0.5, 0.025, 0.975))
postp    <- apply(rbeta > 0, 2, mean)
outMat   <- cbind(t(rbMat), postp)
rownames(outMat) <- colnames(X)
colnames(outMat) <- c(rownames(rbMat), 'P(b > 0)')
round(outMat,4)

cbind(round(t(rbMat),4), apply(rbeta > 0, 2, mean))

## compare ##
lm(Y ~ X - 1)


## Example I.R9 ##
lyme	<- read.table('lymeMN.txt', header = TRUE)
head(lyme)

y		<- lyme$count
B		<- 100000

alpha	<- sum(y)
beta	<- length(y)

set.seed(217)
lambda	<- rgamma(B, alpha, beta)

## Normal Approx ##
mu		<- (sum(y)-1)/length(y)
sig2	<- (sum(y)-1)/(length(y)^2)

set.seed(217)
lamNA	<- rnorm(B, mean = mu, sd = sqrt(sig2))

quantile(lambda, probs = c(0.5, 0.025, 0.975))
quantile(lamNA, probs = c(0.5, 0.025, 0.975))

plot(density(lambda), col = 'blue', lwd = 2, main = 'Posterior Comparison', ylab = expression(paste('p(', lambda,'|y)')), xlab = expression(lambda))
lines(density(lamNA), col = 'gray', lwd = 2, lty = 3)