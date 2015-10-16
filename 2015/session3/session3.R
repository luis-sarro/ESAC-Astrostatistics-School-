setwd("~/Escritorio/ESAC-Astrostats/2015/session3")

# In this session, we will be introducing priors for the first time.
# So, first, a word of warning about priors:
# In the straight line example, we had two parameters: m and c
# What are the implications of a flat prior on m? 
# Let us sample from a flat (but proper) prior
m.min=0
m.max=1.0e3
# We define the prior...
prior.m <- function(m)
{
  pi <- rep(NA,length(m))
  filter.out <- m<0| m> m.max
  filter.in <- !filter.out
  pi[filter.out] <- 0
  pi[filter.in] <- 1/(m.max-m.min)
  return(pi)
}
# ... and plot it
m <- seq(-100,1100,length.out = 1200)
plot(m,prior.m(m),ty="l")

# Then, following what we learnt in session 1, we define the cdf of the prior...
cdf.prior <- function(m)
{
  cdf <- apply(cbind(apply(cbind(m,m.min),1,max)-m.min,m.max),1,min)/m.max
  return(cdf)
}
plot(m,cdf.prior(m),ty="l")

# ... and the inverse of the cdf
icdf.prior <- function(cdf)
{
  m <-cdf*m.max
  return(m)
}

# Now we are ready to generate samples from the prior...
cdf.m <- runif(10000,0,1)
m.sample <- icdf.prior(cdf.m)

# Let us define a function so we can plot the sample from this prior
straight.line <- function(x,m,c)
{
  y <- c+m*x
  return(y)
}

x <- matrix(seq(-1,1,length.out=100),nrow=10000,ncol=100,byrow=T)
y <- matrix(NA,nrow=10000,ncol=100)
for (i in 1:dim(x)[1])
{
  y[i,] <- straight.line(x[i,],m[i],0)
}
plot(x[1,],y[1,],ty="l",xlim=c(-1,1),ylim=c(-1,1))
for (i in 2:dim(x)[1])
{
  lines(x[i,],y[i,],col=rgb(0,0,0,0.5)) # Set transparency with fourth value so we can better discern where do lines accumulate
}

# Lesson to be learnt: even naïve uniform priors are informative priors, and it is important to recall that
# MLE are Bayesian modes for uniform priors.

# Now, let us put in practice what we have learnt.
rm(list=ls()) # To avoid problems 
load("endSession2.RData")

# OK. So let us try the new tools we have just discovered...
# First, let us compute the AIC values for the three models fitted so far
AIC <- rep(NA,3)
AIC[1] <- -2*logLik(mle.order1Poly)+2*2
AIC[2] <- -2*logLik(mle.order25Poly)+2*25
AIC[3] <- -2*logLik(mle.sine)+2*1
print(AIC)
# This suggests that we select the sinusoidal model

# Then, we try with the BIC
BIC <- rep(NA,3)
BIC[1] <- -2*logLik(mle.order1Poly)+log(25)*2
BIC[2] <- -2*logLik(mle.order25Poly)+log(25)*25
BIC[3] <- -2*logLik(mle.sine)+log(25)*1
print(BIC)
# Which again suggests that we select the sinusoidal model

# Finally, let us see what happens with the DIC
# Order 1 polynomial 
# Define a prior on theta instead of m
pi.theta <- function(theta) return(1/(pi/2))
pi.c <- function(c) return(1/(2-1)) # ... and a harmless flat prior on c.
# Define the likelihood for theta and c (instead of m and c)
loglik2.theta <- function(theta,c) 
{
  m <- tan(theta)
  predicted <- x*m+c
  loglik <- dmvnorm(y.noise,predicted,Sigma,log=T) #This time, we do it right, without the minus sign
  return(loglik)
}
# Let us compute the posterior on a grid of parameter values
np <- 100 # Number of points per axis in grid
theta.grid <- seq(0,pi/2,length.out = np)
dtheta <- ((pi/2)-0/np) # The grid step in theta
c.grid <- seq(0,2,length.out = np)
dc <- (2-0)/np # The grid step in c
u.posterior<-matrix(NA,np,np) # Define an empty matrix

# ... and fill it in
for (i in 1:np)
{
  for (j in 1:np)
  {
    u.posterior[i,j] <- loglik2.theta(theta.grid[i],c.grid[j])*pi.theta(theta.grid[i])*pi.c(c.grid[j]) 
  }
}

image(theta.grid,c.grid,exp(u.posterior))
Z <- sum(exp(u.posterior)*dtheta*dc)
posterior <- exp(u.posterior)/Z
print(sum(posterior*dtheta*dc))

#posterior averaged deviance
PAD <- 0
for (i in 1:np)
{
  for (j in 1:np)
  {
    PAD <- PAD + posterior[i,j]*loglik2.theta(theta.grid[i],c.grid[j])*dtheta*dc 
  }
}

# Deviance at posterior mean
# First compute the posterior mean
PM <- c(0,0)
for (i in 1:np)
{
  for (j in 1:np)
  {
    PM <- PM + c(posterior[i,j]*theta.grid[i]*dtheta*dc,posterior[i,j]*c.grid[j]*dtheta*dc) 
  }
}

# but also the mode... 
idcs.max <- which(posterior==max(posterior),arr.ind = T)
post.mode <- c(theta.grid[idcs.max[1]],c.grid[idcs.max[2]])
# ... and the MLE
max.lik <-  c(atan(coef(mle.order1Poly)[1]),coef(mle.order1Poly)[2])

# Visually inspect values
print(PM)
print(post.mode)
print(max.lik)

DPM <- max(loglik2.theta(PM[1],PM[2]),loglik2.theta(post.mode[1],post.mode[2]),loglik2.theta(max.lik[1],max.lik[2]))

Bayesian.complexity <- -2*(PAD-DPM)  
print(Bayesian.complexity)
DIC.linear <- -2*DPM +2* Bayesian.complexity
print(DIC.linear)

# Same for the sine wave
pi.omega <- function(omega) return(1/(2000))
m <- 20000
omega.grid <- seq(0,2000,length.out = m)
domega <- ((2000-0)/m)
u.posterior<- rep(NA,m)

for (i in 1:m)
{
  u.posterior[i] <- -loglik4(omega.grid[i])*pi.omega(omega.grid[i]) 
}

plot(omega.grid,exp(u.posterior),ty="l")
Z <- sum(exp(u.posterior)*domega)
posterior <- exp(u.posterior)/Z
print(sum(posterior*domega))

#posterior averaged deviance
PAD <- 0
for (i in 1:m)
{
  PAD <- PAD + posterior[i]*(-loglik4(omega.grid[i]))*domega 
}

# Deviance at posterior mean

PM <- c(0)
for (i in 1:m)
{
  PM <- PM + posterior[i]*omega.grid[i]*domega 
}

idx.max <- which(posterior==max(posterior))
post.mode <- omega.grid[idx.max]
max.lik <-  coef(mle.sine)

DPM <- max(loglik4(PM),loglik4(post.mode),loglik4(max.lik))

Bayesian.complexity <- -2*(PAD-DPM)  
print(Bayesian.complexity)

DIC.sine <- -2*DPM +2*Bayesian.complexity
print(DIC.sine)

# Proposed exercise: compute the Bayes factor for both the linear model 
# (2D integration) and for the sinusoidal model (1D integration)

