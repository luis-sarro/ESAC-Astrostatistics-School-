

# The likelihood in Eq. (1) of Bailer-Jones (2015)
likelihood <- function(p,d,s_p)
  {
    if (s_p < 0) print("Error. Parallaxes must be positive.")
    const = 1.0/(sqrt(2.0*pi)*s_p)
    XTSX = (1.0/(sqrt(2.0)*s_p))*(p-(1.0/d))
    lik <- const*exp(-XTSX^2)
    return(lik)
}

# Let us do some examples
s_p <- 0.01 # Imagine we have 10 milliarcsecond uncertainties
d.true <- 100 # Measured in parsecs
p.true <- 1/d.true # in arcseconds

# Let us evaluate the likelihood for some potential values of the parallax
nr.values <- 200
values <- seq(p.true*0.5,p.true*1.5,length.out=nr.values)
tmp <- rep(NA,nr.values)
for (i in 1:length(values))
{
tmp[i] <- likelihood(values[i],d.true,s_p)  
}
plot(values,tmp, ty="l", xlab="Parallax", ylab="Likelihood")
# Is this normalized as a function of the measurement?
require(caTools)
trapz(values,tmp)
# Yes (to within numeric precision)

# Now... in terms of distance...
plot(1.0/values,tmp, ty="l", xlab="Distance in parsecs", ylab="Likelihood")

# Play with the values of the uncertainty

# What is the expected uncertainty for Gaia measurements?
# From GAIA ASTROMETRIC SCIENCE PERFORMANCE
# POST-LAUNCH PREDICTIONS
# J.H.J. de Bruijne, K.L.J. Rygl and T. Antoja
# http://arxiv.org/pdf/1502.00791.pdf
Gaia_s_p <- function(VminusI,G)
{
  z = apply(cbind(10^(0.4*(12.09-15)),10^(0.4*(G-15))), 1, max)
  s_p = sqrt(-1.631+680.766*z+32.732*z^2)*(0.986+(1-0.986)*VminusI)
  return(s_p)
}
# Plot the curve
G <- seq(6,21,0.1)
VminusI <- G/G 
plot(G,Gaia_s_p(VminusI,G),ty="l",xlab="G magnitude", ylab="parallax uncertainty in microarcseconds")
# In microarcseconds!

# Example I
VminusI <- 1.0
s_p <- Gaia_s_p(VminusI,20)/10.0^6 # To convert to arcseconds  
d.true <- 100 # Measured in parsecs
p.true <- 1/d.true
nr.values <- 200
values <- seq(p.true*0.5,p.true*1.5,length.out=nr.values)
tmp <- rep(NA,nr.values)
for (i in 1:length(values))
{
  tmp[i] <- likelihood(values[i],d.true,s_p)  
}
plot(values,tmp, ty="l", xlab="Parallax", ylab="Likelihood")
require(caTools)
trapz(values,tmp)

plot(1.0/values,tmp, ty="l", xlab="Distance in parsecs", ylab="Likelihood")

# Unfortunately, not everything is at a distance of 100 parsecs...
# Example II
VminusI <- 1.0
s_p <- Gaia_s_p(VminusI,20)/10.0^6 # To convert to arcseconds  
d.true <- 4000 # Measured in parsecs
p.true <- 1/d.true
nr.values <- 200
values <- seq(p.true*0.5,p.true*1.5,length.out=nr.values)
tmp <- rep(NA,nr.values)
for (i in 1:length(values))
{
  tmp[i] <- likelihood(values[i],d.true,s_p)  
}
plot(values,tmp, ty="l", xlab="Parallax", ylab="Likelihood")
require(caTools)
trapz(values,tmp)

values <- seq(-10*p.true,p.true*10,length.out=nr.values)
tmp <- rep(NA,nr.values)
for (i in 1:length(values))
{
  tmp[i] <- likelihood(values[i],d.true,s_p)  
}
plot(values,tmp, ty="l", xlab="Parallax", ylab="Likelihood")
require(caTools)
trapz(values,tmp)

plot(1.0/values,tmp, ty="l", xlab="Distance in parsecs", ylab="Likelihood")

filter <- (1.0/values > 0)
plot(1.0/values[filter],tmp[filter], ty="l", xlab="Distance in parsecs", ylab="Likelihood")
# And basically, any distance beyond 4kp is likely.

# OK, now let go bayesian.
unnormalized.posterior <- function(prior,likelihood)
  {
  u_posterior <- prior*likelihood
  return(u_posterior)
}

prior.iu<-function(r) 
  {
  if (r <= rlim) return(1/rlim) else return(0)
}


# And, in the philosophy of the paper, let us draw conclusions
# not from a single measurement, but from large samples of sources 
# In order to fix ideas, let us generate samples up to a maximum distance
rlim=1000 #in parsecs
# Recall from session1, how we can draw samples from a 1D PDF
icdf.const.dens <- function(x)
{
  if(!exists("rlim") | is.na(rlim)) print("Error: maximum radius undefined")
  r <- ((rlim^3)*x)^(1.0/3.0) # This is the prescription for a constant stellar volume density
  return(r)
  }

gen.data <- function(n,icdf.spatial.dens)
{
  x <- runif(n)
  r <- icdf.spatial.dens(x)
    return(r)
}

r_sample <- gen.data(1000,icdf.const.dens)

#Run simulation
f <- seq(0,1,0.01) # This is sigma_p/p (see Bailer-Jones, 2015)
# This takes ~5-10 minutes: true_r_sample <- gen.data(1000000,icdf.const.dens)
true_r_sample <- gen.data(1000,icdf.const.dens) # So we limit the number of samples: do not expect the same plots
true_p_sample <- 1/true_r_sample
measured_p_sample <- rep(NA, length(true_r_sample))

bias <- rep(NA,length(f))
sd <- rep(NA,length(f))
for (j in 1:length(f)) # This f is sigma/p
{
  s_p_sample <- f[j]*true_p_sample
  for (i in 1:length(true_p_sample)){
    measured_p_sample[i] <- rnorm(1,true_p_sample[i],s_p_sample[i])
  }
  estimated_r_sample <- 1.0/measured_p_sample # For flat unbounded priors, there is no need to compute the posterior
  filter <- measured_p_sample > 0
  scaled_residuals <- (estimated_r_sample[filter]-true_r_sample[filter])/true_r_sample[filter]
  bias[j] <- mean(scaled_residuals)
  sd[j] <- sd(scaled_residuals)
}

plot(f,bias,ty="l",ylim=c(0,500),ylab="Bias (black) and Std. Dev. (Blue)")
lines(f,sd,ty="l",ylim=c(0,500),col="blue")

# Second experiment
f <- seq(0,1,0.01) # This is sigma_p/p (see Bailer-Jones, 2015)
# This takes ~5-10 minutes: true_r_sample <- gen.data(1000000,icdf.const.dens)
true_r_sample <- runif(1000,0,1000) # So we limit the number of samples: do not expect the same plots
true_p_sample <- 1/true_r_sample
measured_p_sample <- rep(NA, length(true_r_sample))

bias <- rep(NA,length(f))
sd <- rep(NA,length(f))
for (j in 1:length(f)) # This f is sigma/p
{
  s_p_sample <- f[j]*true_p_sample
  for (i in 1:length(true_p_sample)){
    measured_p_sample[i] <- rnorm(1,true_p_sample[i],s_p_sample[i])
  }
  estimated_r_sample <- 1.0/measured_p_sample
  estimated_r_sample[measured_p_sample <0 | estimated_r_sample > rlim] <- rlim
#  filter <- measured_p_sample > 0
  filter <- rep(TRUE,length(measured_p_sample))
  scaled_residuals <- (estimated_r_sample[filter]-true_r_sample[filter])/true_r_sample[filter]
  bias[j] <- mean(scaled_residuals)
  sd[j] <- sd(scaled_residuals)
}

plot(f,bias,ty="l",ylim=c(0,5),ylab="Bias (black) and Std. Dev. (Blue)")
lines(f,sd,ty="l",ylim=c(0,5),col="blue")

# Proposed experiment

# 1 Generate a sample from a constant volume density
# 2 Code a prior that represents constant density (Eq. 9 in Bailer-Jones, 15)
# 3 Compute the (unnormalized) posterior and plot it for several choices of d, NSR and/or measured parallaxes
# 4 Code the values of the mode of the posterior, and compute them for the sample. (Eqs. 16 in Bailer-Jones, 15)
# 5 Plot the bias and variance 

# 1' Generate a sample from an exponentially decreasing density (Eq. 17 in Bailer-Jones, 15; L=10^3)
# 2' Code a prior that represents this exponentially decreasing density 
# 3' Compute the (unnormalized) posterior (Eq 18) and plot it for several choices of d, NSR and/or measured parallaxes
# 4' Compute the roots of the posterior as in page 13-14 of Bailer-Jones, 2015.
# 5' Plot the bias and variance 
# 6' Play with the value of L

