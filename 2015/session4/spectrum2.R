
library(mvtnorm)
lambda <- seq(1810,1950,0.2)  
line_wavelengths <- c(1855,1863,1892,1900,1908)
sigma <- c(3,3,3,3,3)
heights <- c(0.4,0.4,0.3,0.2,0.4)
snr=10
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
#                       2
#######################################################
# Same as 1 with unknown variance (same for all lines)
# That is two parameters
#######################################################

log_likelihood.2 <- function(h4,w,observed)
{
  h <- c(heights[1],heights[2],heights[3],h4,heights[5])
  sigma <- rep(w,length(line_wavelengths))
  noiseless <- generate_noiseless(0.03,0,line_wavelengths,sigma,h)
#  lines(lambda,noiseless,col=i)

  #  likelihood <- dmvnorm(observed,noiseless,covariance)
  #  return(likelihood)
  diff<- as.matrix(observed-noiseless,ncol=1)
  log_likelihood <- -(1/2)*t(diff)%*%solve(covariance)%*%diff
  
  return(log_likelihood)
}

trial.hs <- seq(0.175,0.22,length.out = 21)
trial.w <- seq(2.85,3.1,length.out = 21)
loglik <- matrix(NA,nrow=length(trial.hs),ncol = length(trial.w))
uncert <- flux/snr
covariance <- diag(x=uncert^2)
for (i in 1:length(trial.hs))
{
  for (j in 1:length(trial.w))
  {
#    plot(lambda,flux,ty="l")
#     print(paste(i," ",trial.hs[i], " ", j," ",trial.w[j]))
    loglik[i,j] <- log_likelihood.2(trial.hs[i],trial.w[j],flux) 
#    print(loglik[i,j])
      }
}


image(trial.hs,trial.w,exp(loglik))
contour(trial.hs,trial.w,exp(loglik),add=T,nlevels = 10)
plot(trial.hs,exp(loglik),ty="l")


