#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
#### Asymptotic normality of parameter estimation
#######################################################################
#######################################################################
source("/Users/Buritica/Dropbox/Thèse/git/index_regular_variation/IndexofRV.R")
source("/Users/Buritica/Dropbox/Thèse/git/Auxiliar_functions/random_paths.R")
#######################################################################
#######################################################################
library(evd)
library(extRemes)
library(fExtremes)
library(boot)
library(ggplot2)
library(stable)
sapply(1:3, function(k) stable.set.tolerance(k,1e-200))
#######################################################################
#######################################################################
#######################################################################
## Asymptotic normality for parameter estimation from declustered procedure in
## "Inference for clusters of extreme values" - Ferro, C.A.T. and Segers, J.
## MC simulation scale and shape parameters
paramet <- NULL
N <- 0
n <- 8000
while(N < 500){ 
  path         <- data.frame("sample" = abs( arima.sim(n = n, list(ar=0.8, ma=0), rand.gen=function(n) rt(n,df=1)) ))# ARMAX1(0.5, n))
  #path         <- data.frame("sample" = rfrechet(n))
  look2        <- decluster(path$sample, threshold = quantile( path$sample , 0.95))
  fitDC        <- fevd(y, data = data.frame("y" = c(look2)), threshold = quantile( path$sample , 0.95), type = "GP")
  #y1          <- fevd(sample, path, type = "GP", threshold = quantile( path$sample , 0.95) )
  paramet            <- rbind(paramet,fitDC$results$par)
  N <- N+1
  print(N)
}
## Boxplot and histogram 1/xi
par(mfrow = c(1,2))
paramet <- as.data.frame(paramet)
boxplot(paramet$shape, main = "boxplot shape"); abline(h=1, col = "red")
hist(paramet$shape, main = "histogram shape") ; abline(v=1, col = "red")

boxplot(paramet$scale, main = "boxplot scale")
hist(paramet$scale, main = "histogram scale")

## normality test 
ks.test( (paramet$shape - mean(paramet$shape))/sd(paramet$shape), 'pnorm') ## shape : ok n >>
ks.test( (paramet$scale - mean(paramet$scale))/sd(paramet$scale), 'pnorm') ## scale 
qqnorm((paramet$shape - mean(paramet$shape))/sd(paramet$shape)); abline(a=0,b=1)
qqnorm((paramet$scale - mean(paramet$scale))/sd(paramet$scale)); abline(a=0,b=1)

#ks.test(rnorm(1000),'pnorm')
#######################################################################
#######################################################################
#######################################################################
## Asymptotic normality for the ML estimator of the stable program.
##
##
##
N    <- 0 
thet <- 0.5
n    <- 8000
bu   <- 60
a    <- 1
box <- data.frame("a" = NULL, "b" =NULL, "g" = NULL, "d" = NULL, "Q50" = NULL)
while(N < 250){
  sample <- abs( arima.sim(n = n, list(ar=0.8, ma=0), rand.gen=function(n) rt(n,df=1)) )#ARMAX1((1-thet), n)abs(rstable(n,alpha=1))#  #squaredARCH(lambdaV[1,2],n); a <- alphaestimator(sample) #
  suma   <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
  fit    <- stable.fit.mle.restricted(suma, c(0,1,0,0), c(0,1,0,0))
  quan   <- sapply(c(50), function(k) qstable( (1-bu/(k*365.25)) , 
                                               alpha = fit[1], beta = 1, gamma = fit[3], delta = fit[4]))^(1/a)
  box    <- rbind(box, c(fit, quan))
  N <- N+1
  print(N)
}
real<- qt((1-1/(50*365.25)), df=1) 
## Boxplots
names(box) <- c("a", "b", "g", "d", "Q50")
par(mfrow=c(2,4))
boxplot(box$a, main = "alpha");hist(box$a, main = "alpha")
boxplot(box$g, main = "gamma");hist(box$g, main = "gamma")
boxplot(box$d, main = "delta");hist(box$d, main = "delta")
boxplot( log(box$Q50), main = "log Q50");abline(h = log(real), col="red")
hist(log(box$Q50), main = "log Q50");abline(v = log(real), col="red")


## test
names(box) <- c("a", "b", "g", "d", "Q50")
par(mfrow=c(2,2))

qqnorm((box$a - mean(box$a))/sd(box$a)); abline(0,1)
qqnorm((box$g - mean(box$g))/sd(box$g)); abline(0,1)
qqnorm((box$d - mean(box$d))/sd(box$d)); abline(0,1)
qqnorm((log(box$Q50) - mean(log(box$Q50)))/sd( log(box$Q50) )); abline(0,1)

ks.test( (box$a - mean(box$a))/sd(box$a), 'pnorm')
ks.test( (box$g - mean(box$g))/sd(box$g), 'pnorm')
ks.test( (box$d - mean(box$d))/sd(box$d), 'pnorm')
ks.test( (log(box$Q50) - mean(log(box$Q50)))/sd( log(box$Q50) ), 'pnorm')


#######################################################################
#######################################################################
#######################################################################
## Asymptotic normality for the stable blocks estimator 
##   with automatic block length selection.
##   Fails
N    <- 0 
thet <- 0.5
n    <- 8000
a    <- 1
box <- data.frame("a" = NULL, "b" =NULL, "g" = NULL, "d" = NULL, "bu" = NULL)
while(N < 1000){
  sample <- abs( arima.sim(n = n, list(ar=0.5, ma=0), rand.gen=function(n) rt(n,df=1)) )#ARMAX1((1-thet), n) #squaredARCH(lambdaV[1,2],n); a <- alphaestimator(sample) #
  bu     <- 32
  flag   <- 0
  paramet <- NULL
  while( bu < 256 && flag < 20 ){
    suma        <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
    fit         <- stable.fit.mle.restricted(suma, c(0,1,0,0), c(0,1,0,0))
    bb1         <- stable.fit.mle.restricted(suma, c(1,1,0,0), c(1,1,0,0))
    stat           <- -2*(stable.loglik(suma,bb1)-stable.loglik(suma,fit))
    if(pchisq(stat,1) < 0.95){
      paramet <- rbind(paramet, c( pchisq(stat,1),bu) )
      flag    <- flag + 1
    }  
    bu <- bu +1
  } 
  if(bu < 256){
    min_stat <- which.min(paramet[,1])
    bu       <- paramet[min_stat,2]
  }
  suma   <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
  fit    <- stable.fit.mle.restricted(suma, c(1,1,0,0), c(1,1,0,0))
  quan   <- sapply(c(50), function(k) qstable( (1- bu/(k*365.25)) , 
                                               alpha = 1, beta = 1, gamma = fit[3], delta = fit[4]))^(1/a)
  box    <- rbind(box, c(fit, bu, quan))
  N <- N+1
  print(N)
}

names(box) <- c("a", "b", "g", "d", "bu", "Q50")
par(mfrow=c(2,4))
boxplot(box$g, main = "gamma");hist(box$g, main = "gamma")
boxplot(box$d, main = "delta");hist(box$d, main = "delta")
boxplot(box$bu, main = "bu");hist(box$bu, main = "bu")
boxplot( log(box$Q50), main = "log Q50");abline(h = log(real), col="red")
hist(log(box$Q50), main = "log Q50");abline(v = log(real), col="red")

## test
par(mfrow=c(1,3))
qqnorm((box$g - mean(box$g))/sd(box$g)); abline(0,1)
qqnorm((box$d - mean(box$d))/sd(box$d)); abline(0,1)
qqnorm((log(box$Q50) - mean(log(box$Q50)))/sd( log(box$Q50) )); abline(0,1)

ks.test( (box$g - mean(box$g))/sd(box$g), 'pnorm')
ks.test( (box$d - mean(box$d))/sd(box$d), 'pnorm')
ks.test( (log(box$Q50) - mean(log(box$Q50)))/sd( log(box$Q50) ), 'pnorm')


