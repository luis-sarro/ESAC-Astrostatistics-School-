
rm(list=ls())
library(stats4)
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

noiseless <- generate_noiseless(cont_const,cont_slope,line_wavelengths,sigma,heights)

#######################################################
#                       MLE
#######################################################
# Full problem

log_likelihood.all <- function(h1,h2,h3,h4,h5,w1,w2,w3,w4,w5,cont)
{
  sigma <- c(w1,w2,w3,w4,w5)
  heights <- c(h1,h2,h3,h4,h5)
  #print(sigma)
  #print(heights)
  noiseless <- generate_noiseless(cont,0,line_wavelengths,sigma,heights)
  #lines(lambda,noiseless,col=i)
  
  diff<- as.matrix(flux-noiseless,ncol=1)
#  lines(lambda,diff,col="red")  
#  print(mean(abs(diff)))
#  print(t(diff)%*%diff)
#  print(solve(covariance)%*%diff)  
    log_likelihood <- (1/2)*t(diff)%*%solve(covariance)%*%diff
  #print(log_likelihood)  
  return(log_likelihood)
}



uncert <- noiseless/snr
covariance <- diag(x=uncert^2)

hstart=c(0.4,0.4,0.3,0.2,0.4)
wstart=c(3,3,3,3,3)
contstart = 0.03
mle <- mle(log_likelihood.all,list(h1=hstart[1],h2=hstart[2],h3=hstart[3],h4=hstart[4],h5=hstart[5],
                                   w1=wstart[1],w2=wstart[2],w3=wstart[3],w4=wstart[4],w5=wstart[5],
                                   cont=contstart))

sigma_mle = c(coef(mle)[6:10])
heights_mle = c(coef(mle)[1:5])
const_mle = coef(mle)[11]

mle_noiseless = generate_noiseless(const_mle,0,line_wavelengths,sigma_mle,heights_mle)      
lines(lambda,mle_noiseless,col="orange",lwd=2)

print(logLik(mle))
print(coef(mle))

print(log_likelihood.all(h1=heights[1],h2=heights[2],h3=heights[3],h4=heights[4],h5=heights[5],
                         w1=sigma[1],w2=sigma[2],w3=sigma[3],w4=sigma[4],w5=sigma[5],
                         cont=cont_const))
print(log_likelihood.all(h1=heights_mle[1],h2=heights_mle[2],h3=heights_mle[3],h4=heights_mle[4],h5=heights_mle[5],
                         w1=sigma_mle[1],w2=sigma_mle[2],w3=sigma_mle[3],w4=sigma_mle[4],w5=sigma_mle[5],
                         cont=const_mle))
