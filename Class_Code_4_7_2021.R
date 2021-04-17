##### Example II.R7 #####
library(MCMCpack)
hers	<- read.table('hersreg.txt', header = TRUE)
head(hers)

hp		<- hers[which(hers$treatment == 0),] # placebo
hht		<- hers[which(hers$treatment == 1),] # active treatment
B		  <- 2*10000

np		<- nrow(hp)
yp		<- hp$chtchol
mup		<- mean(yp)

nh		<- nrow(hht)
yh		<- hht$chtchol
muh		<- mean(yh)

sigmas	<- matrix(1, nrow = B, ncol = 2)
mus			<- matrix(0, nrow = B, ncol = 2)
mus[1,]		  <- c(mean(yp), mean(yh)) # raw means as starting values
sigmas[1,]	<- c(var(hp$chtchol), var(hht$chtchol)) # raw variances as starting values

set.seed(1222)
for(b in 2:B){
  ## update sigmas ##
  sigmas[b,1] <- rinvgamma(1, np/2, (1/2)*sum((yp - mus[b-1,1])^2))
  sigmas[b,2] <- rinvgamma(1, nh/2, (1/2)*sum((yh - mus[b-1,2])^2))
  
  ## update mus ##
  mus[b,1] <- rnorm(1, mup, sqrt(sigmas[b-1,1]/np))
  mus[b,2] <- rnorm(1, muh, sqrt(sigmas[b-1,2]/nh))
}

apply(mus[-c(1:(B/2)),], 2, mean)
apply(sigmas[-c(1:(B/2)),], 2, mean)

diff_pvh <- mus[-c(1:(B/2)),1] - mus[-c(1:(B/2)),2]
quantile(diff_pvh, probs = c(0.5, 0.025, 0.975))
mean(diff_pvh > 0)
mean(diff_pvh < 0)

library(mcmcplots)
diffMat <- matrix(diff_pvh, nrow = B)
colnames(diffMat) <- 'Difference'
mcmcplot1(diffMat)

library(coda)
geweke.diag(diffMat)
geweke.diag(sigmas[-c(1:(B/2)),])
geweke.diag(mus[-c(1:(B/2)),])

traplot(sigmas[-c(1:(B/2)),])



### Example II.R8 ###
library(mvtnorm)
library(MCMCpack)
library(mcmcplots)

hers	<- read.table('hersreg.txt', header = TRUE)

set.seed(50)

B	<- 2*10000
k	<- dim(hers)[2]-1
n	<- dim(hers)[1]
y	<- hers[,1]
X <- model.matrix(~ treatment + sbp + statins, data = hers)

sig		  <- rep(0,B)
betamat	<- matrix(0, nrow = B, ncol = k)
bhat		<- c(solve(t(X)%*%X)%*%(t(X)%*%y))
v			  <- solve(t(X)%*%X) 

sig[1]			<- (1/(n-k))*sum((y-X%*%bhat)^2) # MSE
betamat[1,]	<- bhat # OLS

for (b in 2:B) { 
  # Sample Beta #
  betamat[b,] <- c(rmvnorm(1, bhat, sig[b-1]*v)) 
  
  # Sample Sigma^2 #
  sig[b] <- rinvgamma(1, n/2, sum((y-X%*%betamat[b-1,])^2)/2)
}
colnames(betamat)	<- colnames(X)

betai	  <- betamat[(B/2+1):B,]
sigi		<- sig[(B/2+1):B]

round(t(apply(betai, 2, quantile, probs = c(0.5, 0.025, 0.975))), 4)
rmeanplot(betai)
geweke.diag(betai)
par(mfrow = c(1,1))
caterplot(betai, denstrip = TRUE)
caterplot(betai[,-3], denstrip = TRUE) # dropping sbp

mean(sigi)
sigMat				<- matrix(sigi, ncol = 1)
colnames(sigMat)	<- 'sigma^2'
mcmcplot1(sigMat, greek = TRUE)
geweke.diag(sigi)



##### Example II.R9 #####
laplace   <- read.table('laplace.txt', header = TRUE)
x			    <- laplace$x

logPost	  <- function(x, mu){ - sum(abs(x-mu)) }

# mug       <- seq(min(x), max(x), by = 0.001)
mug       <- seq(26.5, 27.75, by = 0.001)
plotmu    <- vector('numeric', length = length(mug))
for(i in 1:length(mug)){
  plotmu[i] <- exp(logPost(x, mug[i]))
}
plot(mug, plotmu, type = 'l', xlab = 'x', ylab = expression(f[X](x)))

B	        <- 2*10000
mu	      <- vector(length = B)
arr	      <- vector(length = B)
propVar	  <- 0.5 # 0.5

### step 1 ###
mu[1]	    <- median(x)

### step 2 ###
set.seed(8562)
for(i in 2:B){
  # muStar	<- rnorm(1, median(x), propVar)
  muStar	<- rnorm(1, mu[i-1], propVar)
  r		    <- exp(logPost(x, muStar))/exp(logPost(x, mu[i-1]))
  # r		<- exp(logPost(x, muStar) - logPost(x, mu[i-1]))
  u		    <- runif(1)
  if(u < min(r,1)){
    mu[i]	  <- muStar
    arr[i]	<- 1
  } else{
    mu[i]	  <- mu[i-1]
    arr[i]	<- 0
  }
}

median(mu[-(1:(B/2))])
hist(mu[-(1:(B/2))])
plot(mu[-(1:(B/2))], type = 'l')
acf(mu[-(1:(B/2))])

mean(arr[-(1:(B/2))])

geweke.diag(mu[-(1:(B/2))])