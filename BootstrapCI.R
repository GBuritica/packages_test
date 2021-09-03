#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
#### Bootstrap CI
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
sapply(1:3, function(k) stable.set.tolerance(k,1e-500))
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
while(N < 500){
  sample <- abs(rstable(n,alpha=1))# abs( arima.sim(n = n, list(ar=0.8, ma=0), rand.gen=function(n) rt(n,df=1)) )#ARMAX1((1-thet), n) #squaredARCH(lambdaV[1,2],n); a <- alphaestimator(sample) #
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


#### Assesing the ML estimator. 
#### It seems to overestimate the Return level
N    <- 0 
thet <- 0.5
n    <- 8000
a    <- 1
box <- data.frame("a" = NULL, "b" =NULL, "g" = NULL, "d" = NULL, "bu" = NULL)
while(N < 500){
    sample <- ARMAX1((1-thet), n) #squaredARCH(lambdaV[1,2],n); a <- alphaestimator(sample) #
    bu     <- 32
    flag   <- 0
    paramet <- NULL
    while( bu < 128 && flag < 20 ){
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
    if(bu < 128){
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
qqnorm(((box$Q50) - mean((box$Q50)))/sd( (box$Q50) )); abline(0,1)


ks.test( (box$g - mean(box$g))/sd(box$g), 'pnorm')
ks.test( (box$d - mean(box$d))/sd(box$d), 'pnorm')
ks.test( ((box$Q50) - mean((box$Q50)))/sd( (box$Q50) ), 'pnorm')




## Assesing bootstrap
library(ecdfHT)
N <- 0 ; thet <- 0.5; n <- 1000; bu <- 60; a <-1
sample <- ARMAX1((1-thet), n); a <- 1
sample <- squaredARCH(lambdaV[1,3],n);a <- 1.2#alphaestimator(sample)
s <- 1 ; i <- 3; a <- 1/0.21364509
sample <- cbind( pre[index[,s],(i+1)] , pre[index[,s],(i+2)] , pre[index[,s],(i+3)] ) 
sample <- na.omit(sample)
n      <- length(sample[,1] )
sample <- sapply(1:n, function(k) max(sample[k,]) ); a <- 1/0.21364509; 
par(mfrow=c(1,1))


bu   <- 1
flag <- TRUE
while( bu < 200 && flag ){
  suma        <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
  fit         <- stable.fit.mle.restricted(suma, c(0,1,0,0), c(0,1,0,0))
  bb1         <- stable.fit.mle.restricted(suma, c(1,1,0,0), c(1,1,0,0))
  t           <- -2*(stable.loglik(suma,bb1)-stable.loglik(suma,fit))
  print(paste0(pchisq(t,1), " ", fit[1]))
  if(pchisq(t,1) < 0.95) flag <- FALSE
  if(flag) bu <- bu +1
}  
  quan   <- sapply(c(100,300), function(k) qstable( (1- bu/(k*365.25)) , 
                                                    alpha = bb1[1], beta = bb1[2], gamma = bb1[3], delta = bb1[4]))
  box    <- data.frame("a" = NULL, "b" =NULL, "g" = NULL, "d" = NULL, "bu" = NULL, "Q100" = NULL, "Q300" = NULL)
  box    <- rbind(box, c(fit, bu, quan) )
  qqplot( sort(suma) , qstable( 1:length(suma)/(1+length(suma)) ,
                               alpha = bb1[1], beta = bb1[2], gamma = bb1[3], delta = bb1[4]))
  
  t.info        <-  ecdfHT( suma, show.axes=FALSE )
  ecdfHT.axes( t.info, x.labels=c(-50,-5,0,5,50), y.labels=c(.001,.01,.1,.5,.9,.99,.999),
               show.vert.gridlines=TRUE, show.horiz.gridline=TRUE, lty=2 )
  q1            <- qstable(t.info$ecdf, alpha=bb1[1], beta=bb1[2], gamma=bb1[3], delta=bb1[4])
  ecdfHT.draw( t.info, q1, t.info$ecdf, col='red', show.ci=TRUE)
  
  
  for(j in 1:500){
    sampleB <- rstable(floor(n/bu), alpha=1, beta=1, gamma=bb1[3], delta=bb1[4])
    fitB    <- stable.fit.mle.restricted(sampleB, c(0,1,0,0), c(0,1,0,0))
    #fitB   <- stable.fit.mle.restricted(sampleB, c(1,1,0,0), c(1,1,0,0))
    quanB   <- sapply(c(100,300), function(k) qstable( (1- bu/(k*365.25)) , 
                                                      alpha = 1, beta = 1, gamma = fitB[3], delta = fitB[4]))
    box     <- rbind(box, c(fitB,bu, quanB))
    print(j)
  }
  
names(box) <- c("a", "b", "g", "d", "bu", "Q100", "Q300")
##### Assessing bootstrap : 
## a parameter
fit <- bb1
p1 <- ggplot(box, aes( ( (a - fit[1]) - mean(a - fit[1]) )/sd((a - fit[1]) ) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : a* - a^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (a - fit[1]) - mean(a - fit[1]) )/sd((a - fit[1]) ) ) ) +
          ggtitle("Bootstrap : a* - a^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
## b parameter
p1 <- ggplot(box, aes(( (b - fit[2]) - mean((b - fit[2])))/sd((b - fit[2])) ) ) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : b* - b^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (b - fit[2]) - mean((b - fit[2])))/sd((b - fit[2])) ) ) +
  ggtitle("Bootstrap : b* - b^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
## g parameter
p1 <- ggplot(box, aes( ( (g - fit[3]) - mean(g - fit[3]))/sd(g - fit[3]) ) ) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : g* - g^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (g - fit[3]) - mean((g - fit[3])))/sd((g - fit[3])) ) ) +
  ggtitle("Bootstrap : g* - g^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
## d parameter
p1 <- ggplot(box, aes( ( (d - fit[4]) - mean((d - fit[4])))/sd((d - fit[4])) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : d* - d^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (d - fit[4]) - mean((d - fit[4])))/sd((d - fit[4])) ) ) +
  ggtitle("Bootstrap : d* - d^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
## Q100 parameter
p1 <- ggplot(box, aes( ( Q100- quan[1] - mean(Q100- quan[1]) )/sd((Q100- quan[1]) ) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue") +ggtitle("Bootstrap : R100* - R100^")
p2 <- ggplot(box) + stat_qq(aes( sample=( Q100- quan[1] - mean(Q100- quan[1]) )/sd((Q100- quan[1]) ) ) ) +
  ggtitle("Bootstrap : R100* - R100^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
## ok. 


##### However, with the log transformation I get normality ...  This is kind of understandable.
p1 <- ggplot(box, aes( ( ( log(Q100)- log(quan[1]) ))/sd(( log(Q100)- log(quan[1]) ) ) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+
  ggtitle("Bootstrap : R100* - R100^") 
p2 <- ggplot(box) + stat_qq(aes( sample=( ( log(Q100)- log(quan[1]) ) - mean( log(Q100)- log(quan[1]) ) )/sd(( log(Q100)- log(quan[1]) ) ) ) ) +
  ggtitle("Bootstrap : R100* - R100^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
y <- ( log(box$Q100) - log(quan[1])  - mean( log(box$Q100)- log(quan[1]) ))/sd(( log(box$Q100)- log(quan[1]) )) 

ks.test( y, "pnorm")
y <- ( (box$g - fit[3]) - mean(box$g - fit[3]))/sd(box$g - fit[3]) 






##### Theoretical confidence intervals
N <- 0 ; thet <- 0.2; n <- 3500;  a <-1; s <- 2; i <- 4;name <- "NE";bu   <- 50
{
sample <- ARMAX1((1-thet), n)
sample <- cbind( pre[index[,s],(i+1)] , pre[index[,s],(i+2)] , pre[index[,s],(i+3)] ) 
sample <- na.omit(sample)
n      <- length(sample[,1] )
sample <- sapply(1:n, function(k) max(sample[k,]) ); a <- alphaestimator(sample); par(mfrow=c(1,1))

flag <- TRUE
while( bu < 150 && flag ){
  suma        <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
  fit         <- stable.fit(suma)
  bb1         <- stable.fit.mle.restricted(suma, c(1,0,0,0), c(1,0,0,0))
  t           <- -2*(stable.loglik(suma,bb1)-stable.loglik(suma,fit))
  if(pchisq(t,1)< 0.95) flag <- FALSE
  if(flag) bu <- bu +1
}
#suma <- rgev(100,xi=0.7,mu=0,beta=1)# 
#suma <- rstable(100,alpha=1.5,beta=1,gamma = 1000,delta = 43000)
fit  <- stable.fit.mle(suma ); n2 <- length(suma); x <- 1:n2/(1+n2)
#fit <- stable.fit(suma); n2 <- length(suma); x <- 1:n2/(1+n2)
q <- NULL
stat <- function(data,i){
  b    <- data[i]
  #b    <- rstable(n2,alpha = fit[1], beta = fit[2], gamma = fit[3], delta = fit[4]) 
  #fit1 <- stable.fit.mle.restricted(b, c(1,1,0,0), c(1,1,0,0))
  fit1  <- stable.fit(b)
  return(qstable(  x , 
                    alpha = fit1[1], beta = fit1[2], gamma = fit1[3], delta = fit1[4]))
}
t     <- boot(suma, statistic = stat, R=500 )
bound <- sapply(1:length(x), function(k) boot.ci(t, type ="perc", index = k )$percent[4:5])
              
quan   <- sapply(c(100,300), function(k) qstable( (1- bu/(k*365.25)) , 
                                                  alpha = fit[1], beta = fit[2], gamma = fit[3], delta = fit[4]))
box    <- data.frame("a" = NULL, "b" =NULL, "g" = NULL, "d" = NULL, "bu" = NULL, "Q100" = NULL, "Q300" = NULL)
box    <- rbind(box, c(fit, bu, quan) )
#qqplot(qstable( 1:length(suma)/(1+length(suma)) ,
#                              alpha = fit[1], beta = fit[2], gamma = fit[3], delta = fit[4]), sort(suma) )

#t.info        <-  ecdfHT( suma, show.axes=FALSE )
#ecdfHT.axes( t.info, x.labels=c(-50,-5,0,5,50), y.labels=c(.001,.01,.1,.5,.9,.99,.999),
#             show.vert.gridlines=TRUE, show.horiz.gridline=TRUE, lty=2 )
#q1            <- qstable(t.info$ecdf, alpha=fit[1], beta=fit[2], gamma=fit[3], delta=fit[4])
#ecdfHT.draw( t.info, q1, t.info$ecdf, col='red', show.ci=TRUE)
## do your own qqplt
qgamma  <- qstable(  x , 
                   alpha = 1,beta=1, gamma = fit[3], delta = fit[4])
qgamma  <- qstable(  x , 
                     alpha = fit[1], beta = fit[2], gamma = fit[3], delta = fit[4])
fqgamma <- dstable(  qgamma , 
                     alpha = fit[1], beta = fit[2], gamma = fit[3], delta = fit[4])
upbound <- (qgamma + qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n2)))
lbound  <- (qgamma - qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n2)))

#bound <- sapply(1:length(x), function(k) c(qgamma[k]+quantile( (q[,k] - qgamma[k]),0.025), qgamma[k]+ quantile((q[,k] - qgamma[k]),0.975) ))
upbound2 <- bound[1,] 
lbound2  <- bound[2,]

plot(sort(suma^(1/a)), qgamma^(1/a), pch =16, cex = 0.5,
     xlab = "Obserseved", ylab = "Theoretical", main = name, 
     ylim =c( quantile(sort(suma^(1/a)),0),quantile(sort(suma^(1/a)),1) ),
     xlim = c( quantile(sort(suma^(1/a)),0),quantile(sort(suma^(1/a)),1) ))
abline(0,1, lty =2, col ="grey")
points( sort(suma^(1/a)), upbound2^(1/a) , col = "blue", pch =16 , cex = 0.5 )
points( sort(suma^(1/a)), lbound2^(1/a), col = "blue", pch =16 , cex = 0.5 )  
segments( sort(suma^(1/a))[1:(length(x)-1)], upbound2[1:(length(x)-1)]^(1/a),
          sort(suma^(1/a))[2:length(x)],upbound2[2:length(x)]^(1/a), col = "blue")
segments( sort(suma^(1/a))[1:(length(x)-1)], lbound2[1:(length(x)-1)]^(1/a),
          sort(suma^(1/a))[2:length(x)],lbound2[2:length(x)]^(1/a), col = "blue")
polygon(c(sort(suma^(1/a)),rev(sort(suma^(1/a))) ),c(upbound2^(1/a),rev(lbound2^(1/a)) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)

lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))

}
## Other plots
{
plot(qgamma, sort(suma), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5, 
     ylim =c( 0,quantile(sort(suma),0.95) ), xlim = c( 0,quantile(sort(suma),0.95) ))
abline(0,1, lty =2, col ="grey")
lines( qgamma, upbound )  
lines( qgamma, lbound )  
polygon(c(qgamma,rev(qgamma) ),c(upbound,rev(lbound) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))

plot(qgamma, sort(suma), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5, 
     ylim =c( 0,quantile(sort(suma),1) ), xlim = c( 0,quantile(sort(suma),1) ))
abline(0,1, lty =2, col ="grey")
lines( qgamma, upbound2 )  
lines( qgamma, lbound2 )  
polygon(c(qgamma,rev(qgamma) ),c(upbound2,rev(lbound2) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)


plot(qgamma, sort(suma), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5, 
     ylim =c( 0,quantile(sort(suma),0.95) ), xlim = c( 0,quantile(sort(suma),0.95) ))
abline(0,1, lty =2, col ="grey")
lines( qgamma, upbound )  
lines( qgamma, lbound )  
polygon(c(qgamma,rev(qgamma) ),c(upbound,rev(lbound) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))

par(mfrow=c(2,2))

plot(qgamma, sort(suma), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5)
abline(0,1, lty =2, col ="grey")
lines( qgamma, upbound )  
lines( qgamma, lbound )  
polygon(c(qgamma,rev(qgamma) ),c(upbound,rev(lbound) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))


plot(sort(suma), qgamma, pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5, 
     ylim =c( 0,quantile(sort(suma),1) ), xlim = c( 0,quantile(sort(suma),1) ))
abline(0,1, lty =2, col ="grey")
lines( sort(suma), upbound2 )  
lines(sort(suma), lbound2 )  
polygon(c(sort(suma),rev(sort(suma)) ),c(upbound2,rev(lbound2) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)

lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))


plot(sort(suma^(1/a)), qgamma^(1/a), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5, 
     ylim =c( 0,quantile(sort(suma^(1/a)),0.95) ), xlim = c( 0,quantile(sort(suma^(1/a)),0.95) ))
abline(0,1, lty =2, col ="grey")
lines( sort(suma^(1/a)), upbound2^(1/a) )  
lines(sort(suma^(1/a)), lbound2^(1/a) )  
polygon(c(sort(suma^(1/a)),rev(sort(suma^(1/a))) ),c(upbound2^(1/a),rev(lbound2^(1/a)) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)

lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))



plot(sort(suma^(1/a)), qgamma^(1/a), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = name, cex=0.5, 
     ylim =c( 0,quantile(sort(suma^(1/a)),1) ), xlim = c( 0,quantile(sort(suma^(1/a)),1) ))
abline(0,1, lty =2, col ="grey")
lines( sort(suma^(1/a)), upbound2^(1/a) )  
lines(sort(suma^(1/a)), lbound2^(1/a) )  
polygon(c(sort(suma^(1/a)),rev(sort(suma^(1/a))) ),c(upbound2^(1/a),rev(lbound2^(1/a)) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)

lm <- lm(qgamma~sort(suma))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))


}
## nice plot

sample <- rnorm(200,1)
n2 <- 200; x <- 1:n2/(1+n2)
qgamma  <- qnorm(  x ,1)
fqgamma <- dnorm(  qgamma ,1)
upbound <- (qgamma + qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n2)))
lbound  <- (qgamma - qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n2)))
plot(qgamma, sort(sample), pch =16,
     xlab = "Obserseved", ylab = "Theoretical", main = "SUD", cex=0.5, 
     #xlim =c(-10,10), ylim= c(-10,10)
)
abline(0,1, lty =2, col ="grey")
lines( qgamma, upbound )  
lines( qgamma, lbound )  
polygon(c(qgamma,rev(qgamma) ),c(upbound,rev(lbound) ), 
        col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
lm <- lm(qgamma~sort(sample))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
mtext(paste0("Stable fit bu = ", bu))
## nice plot

library(extRemes)
qqplot(rstable(100,alpha=1.5,beta=0.5,gamma=0.5,delta=0.5),
       qstable(1:100/101, alpha=1.5,beta=0.5,gamma=0.5,delta=0.5))




       
       