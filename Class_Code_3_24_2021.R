##### Example II.R6 #####
M		  <- 2.56
curve(dbeta(x, 6, 3), from = 0, to = 1, col = 'blue', ylab = 'f(x)', lwd = 2)
curve(M*dunif(x), add = TRUE, lwd = 2)

x1 <- seq(0.1, 0.99, by = 0.001)
max(dbeta(x1, 6, 3))

M		  <- 2.56
B     <- 10000
theta	<- vector(length = B)
b		  <- 1
count	<- 1

set.seed(87)
while(b < (B + 1)){
  ## step 1 ##
  tb	<- runif(1) # proposal
  U	  <- runif(1) # always U(0,1)
  
  ## step 2 ##
  r	<- dbeta(tb, 6, 3)/(M*dunif(tb)) # importance ratio scaled by 1/M
  if(U < r){
    theta[b]	<- tb
    b			    <- b + 1
  }
  count	  <- count + 1
}

b/count
plot(density(theta), ylab = 'f(x)', xlab = 'x')
curve(dbeta(x, 6, 3), add = TRUE, col = 'blue')

plot(1:B, theta, type = 'l')