## Example I.R1 ##
set.seed(1212)
x		<- rnorm(1000)

plot(density(x), lwd = 2, lty = 2, ylim = c(0,0.4))
curve(dnorm(x), add = TRUE, col = 'blue', lwd = 2)



## Example I.R2 ##
set.seed(327)
y		<- rpois(1000,3)

plot(0:15, dpois(0:15, 3), type = 'h', ylim=c(0.,0.25))
points(0:15, dpois(0:15, 3), pch = 16)
lines(0:(length(table(y))-1), table(y)/1000, col = 'blue', type = 'h')
points(0:(length(table(y))-1), table(y)/1000, pch = 16, col ='blue')
