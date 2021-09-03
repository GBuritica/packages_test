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
library(gridExtra)
library(ecdfHT)
sapply(1:3, function(k) stable.set.tolerance(k,1e-100))
#######################################################################
#######################################################################
#######################################################################
## Bootstrap analysis for the stable blocks method
##
#######################
## parameters
N    <- 0 
thet <- 0.5
n    <- 8000
bu   <- 60
a    <- 1
## sample 
sample <- abs( arima.sim(n = n, list(ar=0.8, ma=0), rand.gen=function(n) rt(n,df=1)) ); a <- 1
flag <- 0
paramet <- NULL
## ratio test for block length
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
};print(bu)

## ML parameter estimation
suma   <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
bb1    <- stable.fit.mle.restricted(suma, c(1,1,0,0), c(1,1,0,0))
quan   <- sapply(c(100), function(k) qstable( (1- bu/(k*365.25)) , 
                                                  alpha = bb1[1], beta = bb1[2], gamma = bb1[3], delta = bb1[4]))
box    <- data.frame("a" = NULL, "b" =NULL, "g" = NULL, "d" = NULL, "bu" = NULL, "Q100" = NULL, "Q300" = NULL)
box    <- rbind(box, c(bb1, bu, quan) )

## assesing stable fit with ecdfHT
t.info        <-  ecdfHT( suma, show.axes=FALSE )
ecdfHT.axes( t.info, x.labels=c(-50,-5,0,5,50), y.labels=c(.001,.01,.1,.5,.9,.99,.999),
             show.vert.gridlines=TRUE, show.horiz.gridline=TRUE, lty=2 )
q1            <- qstable(t.info$ecdf, alpha=bb1[1], beta=bb1[2], gamma=bb1[3], delta=bb1[4])
ecdfHT.draw( t.info, q1, t.info$ecdf, col='red', show.ci=TRUE)

### Simulate parametric bootstrap sample from fitted parameters.
for(j in 1:500){
  sampleB <- rstable(floor(n/bu), alpha=1, beta=1, gamma=bb1[3], delta=bb1[4])
  fitB    <- stable.fit.mle.restricted(sampleB, c(0,1,0,0), c(0,1,0,0))
  quanB   <- sapply(c(100), function(k) qstable( (1- bu/(k*365.25)) , 
                                                     alpha = 1, beta = 1, gamma = fitB[3], delta = fitB[4]))
  box     <- rbind(box, c(fitB,bu, quanB))
  print(j)
}
names(box) <- c("a", "b", "g", "d", "bu", "Q100")
##### Assessing bootstrap : 
## a parameter
fit <- bb1
p1 <- ggplot(box, aes( ( (a - fit[1]) - mean(a - fit[1]) )/sd((a - fit[1]) ) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : a* - a^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (a - fit[1]) - mean(a - fit[1]) )/sd((a - fit[1]) ) ) ) +
  ggtitle("Bootstrap : a* - a^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
ks.test( ( (box$a - fit[1]) - mean(box$a - fit[1]) )/sd((box$a - fit[1]) ), 'pnorm')

## g parameter
p1 <- ggplot(box, aes( ( (g - fit[3]) - mean(g - fit[3]))/sd(g - fit[3]) ) ) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : g* - g^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (g - fit[3]) - mean((g - fit[3])))/sd((g - fit[3])) ) ) +
  ggtitle("Bootstrap : g* - g^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
ks.test( ( (box$g - fit[3]) - mean(box$g - fit[3]) )/sd((box$g- fit[3]) ), 'pnorm')

## d parameter
p1 <- ggplot(box, aes( ( (d - fit[4]) - mean((d - fit[4])))/sd((d - fit[4])) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue")+ggtitle("Bootstrap : d* - d^")
p2 <- ggplot(box) + stat_qq(aes( sample=( (d - fit[4]) - mean((d - fit[4])))/sd((d - fit[4])) ) ) +
  ggtitle("Bootstrap : d* - d^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
ks.test( ( (box$d - fit[4]) - mean(box$d - fit[4]) )/sd((box$d - fit[4]) ), 'pnorm')

## Q100 parameter
p1 <- ggplot(box, aes( ( Q100- quan[1] - mean(Q100- quan[1]) )/sd((Q100- quan[1]) ) )) + geom_histogram(aes(y =..density..), bins =40) +stat_function(fun = dnorm, col ="blue") +ggtitle("Bootstrap : R100* - R100^")
p2 <- ggplot(box) + stat_qq(aes( sample=( Q100- quan[1] - mean(Q100- quan[1]) )/sd((Q100- quan[1]) ) ) ) +
  ggtitle("Bootstrap : R100* - R100^") +geom_abline(intercept=0,slope=1)
grid.arrange(p1,p2,ncol=2)
ks.test( ( (box$Q100 - quan[1]) - mean(box$Q100 - quan[1]) )/sd((box$Q100 - quan[1]) ), 'pnorm')
## 


##### Theoretical confidence intervals
N <- 0 ; thet <- 0.2; n <- 8000;  a <- 5; name <- "ARMAX";bu   <- 50
{
  ### block length selection
  sample <- abs( arima.sim(n = n, list(ar=0.8, ma=0), rand.gen=function(n) rt(n,df=a)) )#ARMAX1((1-thet), n, al = a)
  paramet <- NULL
  flag    <- 0
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
  
  suma        <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
  fit         <- stable.fit.mle.restricted(suma, c(1,1,0,0), c(1,1,0,0))
  n2          <- length(suma); x <- 1:n2/(1+n2)
  
  
  ### parametric bootstrap
  ran.gen.stable <- function(data, param) rstable(length(data), alpha=1,beta=1,gamma=param$gam,delta=param$del)
  ## Bootstrap Statistic
  statisticQ     <- function(data, RV=RP){
    d            <- data 
    st           <- stable.fit.mle.restricted(d, c(0,1,0,0), c(0,1,0,0)) 
    res          <- qstable(p= x,  alpha=  st[1], beta= 1, gamma = st[3], delta=st[4])
    return(res) 
  }

  y1             <-  boot(suma, statistic = statisticQ, R = 50, sim="parametric",ran.gen=ran.gen.stable,
                          mle = list(gam=fit[3],del=fit[4]) )  
  
  
  q     <- NULL
  bound <- sapply(1:length(x), function(k) boot.ci(y1, type ="perc", index = k )$percent[4:5])
  box    <- rbind(box, c(fit, bu) )

  ## Theoretical bounds
  qgamma  <- qstable(  x ,  alpha = 1,beta=1, gamma = fit[3], delta = fit[4])
  fqgamma <- dstable( qgamma , alpha = fit[1], beta = fit[2], gamma = fit[3], delta = fit[4])
  
  ## CI for quantiles 
  upbound <- (qgamma + qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n2)))
  lbound  <- (qgamma - qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n2)))
  
  #bound <- sapply(1:length(x), function(k) c(qgamma[k]+quantile( (q[,k] - qgamma[k]),0.025), qgamma[k]+ quantile((q[,k] - qgamma[k]),0.975) ))
  
  ## bootstrap CI
  upbound2 <- bound[1,] 
  lbound2  <- bound[2,]
  
  ## linear regression
  lm <- lm(qgamma~sort(suma^(1/a)))
  
  
  ## plot bootstrap CI
  plot(sort(suma^(1/a)), qgamma^(1/a), pch =1, cex = 0.5,
       xlab = "Obserseved", ylab = "Theoretical", main = name, 
       ylim = c(  quantile(sort(suma^(1/a)),0),quantile(sort(suma^(1/a)),1) ),
       xlim = c(  quantile(sort(suma^(1/a)),0),quantile(sort(suma^(1/a)),1) ))
  abline(0,1, lty =1, col ="grey")
  
  lines( sort(suma^(1/a)), upbound2^(1/a) , col = "grey", pch =16 , lty=2)
  lines( sort(suma^(1/a)), lbound2^(1/a), col = "grey", pch =16 , lty=2 ) 
  
  #abline( (lm$coefficients) ,lty =2, col ="orange", cex = 0.5)
  legend("topleft", legend= c("1-1 line"), col =c("grey"), lty=2 , title = "95% bootstrap CI ")
  mtext(paste0("Stable fit bu = ", bu))
  
  ## Plot quantile theoretical CI
  plot(sort(suma^(1/a)), qgamma^(1/a), pch =1, cex = 0.5,
       xlab = "Obserseved", ylab = "Theoretical", main = name, 
       ylim = c(  quantile(sort(suma^(1/a)),0),quantile(sort(suma^(1/a)),1) ),
       xlim = c(  quantile(sort(suma^(1/a)),0),quantile(sort(suma^(1/a)),1) ))
  abline(0,1, lty =1, col ="grey")
  
  lines( sort(suma^(1/a)), upbound^(1/a) , col = "grey", pch =16 , lty=2 )
  lines( sort(suma^(1/a)), lbound^(1/a), col = "grey", pch =16 , lty=2) 
  
  #polygon(c(sort(suma^(1/a)),rev(sort(suma^(1/a))) ),c(upbound^(1/a),rev(lbound^(1/a)) ), 
  #        col =  adjustcolor("darkgrey", alpha.f = 0.30), border = NA)
  
  #abline( (lm$coefficients) ,lty =2, col ="orange")
  legend("topleft", legend= c("1-1 line"), col =c("grey"), lty=2 , title = "95% CI theoretical quantiles")
  mtext(paste0("Stable fit bu = ", bu))
  
}


##### Example theoretical CI for quantiles: normal distribution
n <- 50
sample <- rnorm(n,1)
x <- 1:n/(1+n)
qgamma  <- qnorm(  x ,1)
fqgamma <- dnorm(  qgamma ,1)
upbound <- (qgamma + qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n)))
lbound  <- (qgamma - qnorm(0.975)*sqrt(x*(1-x))/(fqgamma*sqrt(n)))
plot(qgamma, sort(sample), pch =1,
     xlab = "Obserseved", ylab = "Theoretical", main = "qqnorm", cex=0.5) 
     #xlim =c(-10,10), ylim= c(-10,10))
abline(0,1, lty =1, col ="grey")
lines( qgamma, upbound, lty=2,col="grey" )  
lines( qgamma, lbound , lty=2,col="grey" )  
lm <- lm(qgamma~sort(sample))
abline( (lm$coefficients) ,lty =2, col ="orange")
legend("topleft", legend= c("1-1 line", "regression"), col =c("grey","orange"), lty=2 , title = "95% confidence bands")
## nice plot




       
       