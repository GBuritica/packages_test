# extRemes
## Parameter estimation from 'GP' declustered procedure :
    
    fevd(y, data = data.frame("y" = c(look2)), threshold = quantile( path$sample , 0.95), type = "GP")
 
### Comments:
For the ARMAX, AR(1) model 
 - Shape parameter normally distributed.
 - Scale parameter not normally distributed : histrogram, qqplot, ks.test provided.
 
 For the iid models : Frechet, Burr, t-student ok.
  - Shape parameter normally distributed.
  - Scale parameter normally distributed.
  
  
  # Stable 
  ## Parameter estimation ML method.
  
    suma   <- sapply(1:floor(n/bu), function(k) sum( sample[((k-1)*bu +1):(k*bu)]^a ))
    fit    <- stable.fit.mle.restricted(suma, c(0,1,0,0), c(0,1,0,0))
    
### Comments:  
For the ARMAX, AR(1) models: stable, gamma, delta parameters normally distributed.

