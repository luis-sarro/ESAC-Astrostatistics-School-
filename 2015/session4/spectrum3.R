
library(mvtnorm)
lambda <- seq(1810,1950,length.out = 300)  
line_wavelengths <- c(1855,1863,1892,1900,1908)
sigma <- c(3,3,3,3,3)
heights <- c(0.4,0.4,0.3,0.2,0.4)
snr=3
set.seed(222)

cont_const <- 0.03
cont_slope <- 0

generate_data <- function(cont_const,cont_slope)
{
  x <- seq(0,1,length.out = length(lambda))
  flux_c <- cont_const+x*cont_slope

  flux_l <- rep(0,length(lambda))
  for (i in 1:length(line_wavelengths))
  {
    flux_l <- flux_l+heights[i]*dnorm(lambda,line_wavelengths[i],sigma[i])
  }
  
  flux <- flux_c+flux_l
  noiseless <- flux
  uncert <- flux/snr

  for (i in 1:length(lambda))
  {
    flux[i] <- rnorm(1,flux[i],uncert[i])
  }
  
  plot(lambda,flux,ty="l")
  lines(lambda,noiseless,col="blue", lwd=2)

  return(flux)
}

flux <- generate_data(0.03,0.0)

generate_noiseless <- function(cont_const,cont_slope,line_wavelengths,sigma,heights)
{
  x <- seq(0,1,length.out = length(lambda))
  flux_c <- cont_const+x*cont_slope
  
  flux_l <- rep(0,length(lambda))
  for (i in 1:length(line_wavelengths))
  {
    flux_l <- flux_l+heights[i]*dnorm(lambda,line_wavelengths[i],sigma[i])
  }
  flux <- flux_c+flux_l
  return(flux)
}

#######################################################
#                       3
#######################################################
# Same as 2 but inferring also the continuum level
# That is three parameters
#######################################################

log_likelihood.3 <- function(h4,w,cont,observed)
{
  
  h <- c(h4*2,h4*2,h4*3/2,h4,h4*2) # This is the key difference
  sigma <- rep(w,length(line_wavelengths))
  noiseless <- generate_noiseless(cont,0,line_wavelengths,sigma,h)
#  lines(lambda,noiseless,col=i)

  #  likelihood <- dmvnorm(observed,noiseless,covariance)
  #  return(likelihood)
  diff<- as.matrix(observed-noiseless,ncol=1)
  log_likelihood <- -(1/2)*t(diff)%*%solve(covariance)%*%diff
  
  return(log_likelihood)
}

# Values for snr=10
#trial.hs <- seq(0.185,0.21,length.out = 21)
#trial.w <- seq(2.85,3.1,length.out = 21)
#trial.cont <- seq(0.0285,0.03,length.out = 21)
# Values for snr=3
trial.hs <- seq(0.1,0.3,length.out = 21)
trial.w <- seq(2.0,4.0,length.out = 21)
trial.cont <- seq(0.01,0.05,length.out = 21)

loglik <- array(NA,dim=c(length(trial.hs),length(trial.w),length(trial.cont)))
uncert <- flux/snr
covariance <- diag(x=uncert^2)
for (i in 1:length(trial.hs))
{
  for (j in 1:length(trial.w))
  {
    for (k in 1:length(trial.cont))
    {
      #    plot(lambda,flux,ty="l")
      #     print(paste(i," ",trial.hs[i], " ", j," ",trial.w[j]))
      loglik[i,j,k] <- log_likelihood.3(trial.hs[i],trial.w[j],trial.cont[k],flux) 
      #    print(loglik[i,j])
    }
  }
}

image(trial.hs,trial.cont,exp(loglik[,7,]))
contour(trial.hs,trial.cont,exp(loglik[,7,]),add=T,nlevels = 10)


