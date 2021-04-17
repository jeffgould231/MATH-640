##### Example II.R1 #####
set.seed(1875)

B 		<- 10000
rtheta	<- rbeta(B, 8, 14)
yobs	<- c(1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

yrep	<- matrix(0, nrow = B, ncol = 20)
Tyrep	<- vector(length = B)

for(b in 1:B){
  yrep[b,] <- rbinom(length(yobs), 1, rtheta[b])
  Tyrep[b] <- length(rle(yrep[b,])$lengths)-1
}

Ty <- 3
hist(Tyrep, xlab = expression(paste(T,'(',y^rep,',',theta,')')),
     col = 'darkgray')
abline(v=Ty, col = 'blue', lwd = 3)

mean(Tyrep)
median(Tyrep)
1-mean(Tyrep > Ty)


##### Example II.R2 #####
library(MCMCpack)
library(moments)

hers <- read.table('hers.txt', header = TRUE)
head(hers)

y		<- hers$x
ybar	<- mean(y)
sig		<- var(y)
n		<- length(y)
B		<- 10000

set.seed(12345)

rmu		<- ybar + (sqrt(sig/n))*rt(B, df = n-1)
rsig	<- rinvgamma(B, (n-1)/2, ((n-1)/2)*sig)

yrep	<- matrix(0, nrow = B, ncol = length(y))
Tyrep	<- matrix(0, nrow = B, ncol = 4)

for(b in 1:B){
  yrep[b,]	<- rnorm(length(y), rmu[b], sqrt(rsig[b]))
  Tyrep[b,1]	<- mean(yrep[b,])
  Tyrep[b,2]	<- var(yrep[b,])
  Tyrep[b,3]	<- skewness(yrep[b,])
  Tyrep[b,4]	<- kurtosis(yrep[b,])
}

Ty1	<- mean(y)
Ty2 <- var(y)
Ty3	<- skewness(y)
Ty4	<- kurtosis(y)

1-sum(Tyrep[,1] > Ty1)/B
1-sum(Tyrep[,2] > Ty2)/B
1-sum(Tyrep[,3] > Ty3)/B
sum(Tyrep[,4] > Ty4)/B

mean(Tyrep[,4])
median(Tyrep[,4])
hist(Tyrep[,4], xlim = c(min(Tyrep[,4]), 6), col = 'darkgray')
abline(v=Ty4, col = 'blue', lwd = 3)


##### Example II.R3 #####
library(mcmcplots)

rmuMat				<- matrix(rmu, ncol = 1)
colnames(rmuMat)	<- 'mu'

rmuList	<- mcmc.list(list(mcmc(rmu)))

traplot(rmuMat, greek = TRUE)
rmeanplot(rmuMat, greek = TRUE)
autplot1(rmuList)
mcmcplot1(rmuMat, greek = TRUE)


##### Example II.R4 #####
library(coda)

geweke.diag(mcmc(rmu))

rmu1	<- mcmc(rmu)

set.seed(125263)
rmu2	<- mcmc(ybar + (sqrt(sig/n))*rt(B, df = n-1))

set.seed(542)
rmu3	<- mcmc(ybar + (sqrt(sig/n))*rt(B, df = n-1))

set.seed(8686)
rmu4	<- mcmc(ybar + (sqrt(sig/n))*rt(B, df = n-1))

allChains <- mcmc.list(list(rmu1, rmu2, rmu3, rmu4)) 

gelman.diag(allChains)


##### Example II.R5 #####
# from Example II.R1 #
yobs	<- c(1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)

set.seed(1875)

B 		<- 10000
rtheta	<- rbeta(B, 8, 14)

n			<- length(yobs)
X			<- sum(yobs)
thb 		<- mean(rtheta)

binomLik	<- function(n,x,p){ choose(n, x)*(p^x)*(1-p)^(n-x)}

pDIC		<- 2*(log(binomLik(n, X, thb)) - (1/B)*sum(log(binomLik(n, X, rtheta))))
DIC			<- -2*log(binomLik(n, X, thb)) + 2*pDIC

# using built-in functions #

pDIC		<- 2*(log(dbinom(X, n, thb)) - (1/B)*sum(log(dbinom(X, n, rtheta))))
DIC			<- -2*log(dbinom(X, n, thb)) + 2*pDIC

# from Example II.R2 #
hers <- read.table('hers.txt', header = TRUE)

y		<- hers$x
ybar	<- mean(y)
sig		<- var(y)
n		<- length(y)
B		<- 10000

set.seed(12345)

rmu		<- ybar + (sqrt(sig/n))*rt(B, df = n-1)
rsig	<- rinvgamma(B, (n-1)/2, ((n-1)/2)*sig)

lppdv		<- vector(length = n)
pwaicv		<- vector(length = n)

for(i in 1:n){
  yi			<- y[i]
  liki		<- dnorm(yi, rmu, sqrt(rsig)) # built-in function for likelihood of y
  lppdv[i]	<- 1/B*sum(liki)
  pwaicv[i]	<- (1/(B-1))*sum((log(liki) - mean(log(liki)))^2)
}

lppd	<- sum(log(lppdv))
pwaic	<- (1/(B-1))*sum(pwaicv)
WAIC	<- -2*lppd +2*pwaic