require("rstan")

lambda <- seq(1810,1950,length.out = 100)  
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

flux <- generate_data(0.03,0.0)
noiseless <- generate_noiseless(cont_const,cont_slope,line_wavelengths,sigma,heights)

uncert <- noiseless/snr
covariance <- diag(x=uncert^2)
N <- length(lambda)

# Store the data and coefficients in a list to be passed to stan
data <- list(N=N,nlines=5,center=line_wavelengths,lambda=lambda,observed=flux,Sigma=covariance)

# Giving a good initial guess helps convergence
initf <- function() 
  {
  list(h1= 0.4, h2=0.4, h3=0.3, h4=0.2, h5=0.4,
    w1= 3, w2=3, w3=3, w4=3, w5=3,
    cont_const=0.03,cont_slope=0)
}     

# Start sampling
fit <- stan(file="spectrum.mod",data=data,init=initf,iter=500,chains=3)
#fit <- stan(fit,data=data,init=initf,iter=200,chains=5)

fit.man <- extract(fit,inc_warmup=T)



