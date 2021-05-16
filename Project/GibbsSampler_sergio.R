library(tidyverse)

omega <- 0.633
theta <- 0.327
xi <- 0.514


btc_candles <- read_csv("btcOpenCandles.csv") %>%
  mutate(dVt = lead(AnnualizedRV) - AnnualizedRV) 

test <- btc_candles %>%
  filter(time >= as.Date("2019-07-01"))
train <- btc_candles %>%
  filter(time < as.Date("2019-07-01"))

Vt <- train$AnnualizedRV
dVt <- train$dVt


omega_func <- function(Omega, Theta, Xi, q, r, .dVt = dVt, .Vt = Vt){
  exp(-1/2*Xi^(2)*sum(((.dVt - Theta * (Omega - .Vt))^2) / .Vt)-q*Omega)*Omega^(r-1)
}

#logL_theta <- function(Theta, Omega, Xi,  .dVt = dVt, .Vt = Vt){
# -1 / (2 * Xi^2) * sum(((.dVt - Theta * (Omega - .Vt))^2) / .Vt)
#}

burnIn <- 5000
thin <- 10
keep <- 1000
B <- burnIn + keep * thin

omega <- vector("numeric", B)
xi <- vector("numeric", B)
theta <- vector("numeric", B)

omega[1] = mean(Vt)
xi[1] = sd(Vt)^2
theta[1] = 0.20

n <- length(Vt)
a<-4
d<-1
s<-3
t<-2
tune<-0.13

omega_accept <- 0
theta_accept <- 0
set.seed(111)
for (b in 2:B) {
  
  
  ##Gibbs for xi
  phi_alpha <- (n/2)+a
  phi_beta <- (sum(((dVt - theta[b-1]*(omega[b-1] - Vt))^(2))/(Vt))-2*d)/2
  
  phi <- MCMCpack::rinvgamma(1, phi_alpha, phi_beta)
  xi[b] = sqrt(phi)
  
  U <- runif(1)
  
  #Gibbs for theta
  mu<- sum((dVt*(omega[b-1]-Vt)))/sum((omega[b-1]-Vt)^(2))
  sd<-sqrt(xi[b-1]^(2)*(sum((omega[b-1]-Vt)^(2)/Vt)^(-1)))
  
  theta_new<- rnorm(1,mean=mu,sd=sd)
  theta[b] = theta_new
  
  ## MA for omega
  omega_star <- rnorm(1, mean = mean(Vt), sd = tune*sd(Vt))
  
  r <- omega_func(omega_star, theta[b-1], xi[b-1],s,t)/omega_func(omega[b-1], theta[b-1], xi[b-1],s,t)
  
  if(U < min(r,1)){
    omega[b] = omega_star
    if(b > burnIn){
      omega_accept <- omega_accept + 1
    }
  }else{
    omega[b] = omega[b-1]
  }
  
}



HestonModel <- function(dBt, theta, omega, Vt, xi){
  return(
    theta * (Vt - omega) + xi * sqrt(Vt) * dBt
  )
}


simDVT <- function(Vt_i, theta, omega, xi, date, dBt = NULL, sims = 1000){
  if(is.null(dBt)) dBt <- rnorm(sims)
  
  dYt <- t(sapply(dBt, HestonModel, Vt_i, theta, omega, xi)) %>% as.data.frame() 
  dYt$Date = date
  return(dYt)
  
}


ParamsDf <- data.frame(theta = theta,
                       omega =omega,
                       xi = xi) %>%
  expand_grid(test%>% select(time, AnnualizedRV)) %>%
  rename(date = time,
         Vt_i = AnnualizedRV)

set.seed(111)
testSim <- pmap_dfr(ParamsDf, simDVT)

checkResults <- testSim %>%
  pivot_longer(cols = V1:V1000, values_to = "simDVt") %>%
  select(-name) %>%
  left_join(test %>% select(Date = time, dVt, AnnualizedRV) %>% mutate(Vt1 = lead(AnnualizedRV))) %>%
  mutate(Vt1_pred = AnnualizedRV + simDVt)

write_rds(checkResults, "HestonModel2Sim.rds")

checkResults %>%
  group_by(Date) %>%
  mutate(upper95 = quantile(Vt1_pred, 0.95),
         lower95 = quantile(Vt1_pred, 0.05)) %>%
  mutate(upper = max(Vt1_pred),
         lower = min(Vt1_pred),
         pointEst = mean(Vt1_pred)) %>%
  select(Date, upper95, lower95, Vt1, Vt1_pred) %>%
  distinct() %>%
  mutate(in95 = (Vt1 >=lower95 & Vt1 <= upper95)) %$%mean(in95, na.rm = T) #%$% cor(Vt1_pred, Vt1, use = "complete.obs")^2

checkResults %>%
  select(Date, Vt1_pred, Vt1) %>%
  ungroup() %>%
  group_by(Date, Vt1) %>%
  summarise(pointEst = mean(Vt1_pred)) %$% cor(Vt1, pointEst, use = "complete.obs")^2

### WAIC
keepIdx <- seq(burnIn+1, B, by = thin)
ThetaKeep <- theta[keepIdx]
OmegaKeep <- omega[keepIdx]
XiKeep <- xi[keepIdx]

GibbsSamplerResults2 <- data.frame(omega = OmegaKeep,
                                   theta = ThetaKeep,
                                   xi = XiKeep)

saveRDS(GibbsSamplerResults2, "GibbsSamplerHeston2.rds")

GibbsSamplerResultsHeston <- read_rds("GibbsSamplerHeston2.rds")

WAIC <- test %>%
  select(Date, dVt, AnnualizedRV) %>%
  mutate(Vt_1 = dVt + AnnualizedRV,
         i = row_number()) %>%
  expand_grid(GibbsSamplerResultsHeston%>% mutate(b = 1:keep)) %>%
  mutate(mu = AnnualizedRV + theta*(omega - AnnualizedRV),
         sigma2 = xi*sqrt(AnnualizedRV)) %>%
  mutate(Ly_i = dnorm(Vt_1, mu, sqrt(sigma2)))

lppd <- WAIC %>% group_by(i) %>% summarise(BB = 1/n() * sum(Ly_i)) %$% sum(log(BB), na.rm = T)

pWAIC <-  WAIC %>%
  group_by(i) %>%
  mutate(meanLogl = 1/keep * sum(log(Ly_i))) %>%
  summarise(inside = sum((log(Ly_i) - meanLogl)^2)) %$% sum((1/(keep-1)) * inside, na.rm = T)

WAIC_return <- -2 * lppd + 2 * pWAIC

