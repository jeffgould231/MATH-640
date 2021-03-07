## Example I.R10 ##
library(mvtnorm)
oldFaithful	<- read.table('oldFaithful.txt', header = TRUE)

fit		<- glm(etime ~ waiting, data = oldFaithful, family = binomial(link = 'logit'))
bhat	<- coef(fit)
vbeta	<- vcov(fit)
B		<- 10000

beta	<- rmvnorm(B, mean = bhat, sigma = vbeta)

round(t(apply(beta, 2, quantile, probs = c(0.5, 0.025, 0.975))), 4)


## Example I.R11 ##
seizure	<- read.table('seizures.txt', header = TRUE)

fit		<- glm(y ~ trt + base + age, data = seizure, family = poisson(link = 'log'))
bhat	<- coef(fit)
vbeta	<- vcov(fit)
B		<- 10000

beta	<- rmvnorm(B, mean = bhat, sigma = vbeta)

round(t(apply(beta, 2, quantile, probs = c(0.5, 0.025, 0.975))), 4)