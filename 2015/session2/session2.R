#set.seed(10)

setwd("~/Escritorio/ESAC-Astrostats/2015/session2")
rm(list=ls()) # Remove any previous variable
require("mvtnorm") # Load library for multivariate normals

load("dataset.RData") # Load the dataset
n <- length(x) # How many points?
sigma <- 0.1 # Define measurement uncertainties
plot(x,y,pch=15,col="blue")
y.noise <- y

# Define the log of the Likelihood 
loglik <- function(x.obs,y.obs,m,c,Sigma)
{
  predicted <- x.obs*m+c
  loglik <- dmvnorm(y.obs,predicted,Sigma,log=T) # Note the argument log=T
  return(loglik)
  }

Sigma <- diag(rep(sigma,length(x))) # In this case, we assume iid, but the general case is a full covariance matrix

# What would be the log-likelihood for the true parameters?
m <- 1
c <- 1
loglik(x,y,m,c,Sigma)

# And for any other set of parameters?
# Let us define a grid of values
m.grid <- seq(0,2,length.out = 100)
c.grid <- seq(0,2,length.out = 100)
loglik.grid <- matrix(NA,length(m.grid),length(c.grid))
# ... compute the log-likelihood on the grid
for (i in 1:length(m.grid))
{
  for (j in 1:length(c.grid))
  {
    loglik.grid[i,j] <- loglik(x,y.noise,m.grid[i],c.grid[j],Sigma)  
  }
}
# And plot it
par(mar=c(5,5,1,8))
image(m.grid,c.grid,exp(loglik.grid)) # Plot the likelihood, not the log-lik. Note the exp function!
require(fields)
image.plot(add=T,legend.only=TRUE, zlim= range(exp(loglik.grid)),col=heat.colors(12),legend.width = 3)

# Where is the maximum?
# There is one R package that will search the parameter space for the maximum
require(stats4)
# But we need to define *minus!* the log-likelihood
# and only in terms of the parameters:
loglik2 <- function(m,c)
{
  predicted <- x*m+c
  loglik <- -dmvnorm(y.noise,predicted,Sigma,log=T) # In reality, this is MINUS the LogLikelihood
  return(loglik)
}
# Once we have this, we can use the mle function to search for the maximum
mle <- mle(loglik2,start=list(m=1,c=1))
# And print the results
print(mle)
coef(mle)
logLik(mle)

# PROPOSED EXERCISE: Do all of the above for the case with severely underestimated uncertainties.

# OK.
# Now, let us do the same with a much more complex model: a polynomial of the 24th degree 
loglik3 <- function(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
                    m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,
                    m21,m22,m23,m24,c)
{
  m <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24)
  predicted <- c
  for (i in 1:length(m))
  {
    predicted <- predicted + x^i*m[i]
  }
  loglik <- -dmvnorm(y.noise,predicted,Sigma,log=T) # In reality, this is MINUS the LogLikelihood
  return(loglik)
}

mle <- mle(loglik3,list(m1=0, m2=0, m3=0, m4=0,m5=0,m6=0,m7=0,m8=0,m9=0,m10=0,
                        m11=0, m12=0, m13=0, m14=0,m15=0,m16=0,m17=0,m18=0,m19=0,m20=0,
                        m21=0, m22=0, m23=0, m24=0,c=0.))
print(mle)
coef(mle)
logLik(mle)

# The log-likelihood is higher! Why?
# Let's plot the model that we just fitted:
plot(x,y.noise,pch=15,col="orange",xlim=c(-.1,1.1),ylim=c(.8,2.2))
# First, we have to define a function so that we can plot the model for any set of parameters
mle.model <- function(x,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,
                      m11,m12,m13,m14,m15,m16,m17,m18,m19,
                      m20,m21,m22,m23,m24,c) 
  {
  m <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24)
  predicted <- c
  for (i in 1:length(m))
  {
    predicted <- predicted + x^i*m[i]
  }
  return(predicted)
  }

# Let us sort the vectors in increasing order of the abscissa. 
  p <- order(x)
  lines(x[p],mle.model(x,coef(mle)[1],coef(mle)[2],coef(mle)[3],coef(mle)[4],coef(mle)[5],
                     coef(mle)[6],coef(mle)[7],coef(mle)[8],coef(mle)[9],coef(mle)[10],
                     coef(mle)[11],coef(mle)[12],coef(mle)[13],coef(mle)[14],coef(mle)[15],
                     coef(mle)[16],coef(mle)[17],coef(mle)[18],coef(mle)[19],coef(mle)[20],
                     coef(mle)[21],coef(mle)[22],coef(mle)[23],coef(mle)[24],coef(mle)[25])[p],col="blue")

# But wait... should there not be a solution that passes exactly through the 25 points?

A <- cbind(rep(1,length(x)),x,x^2,x^3,x^4,x^5, x^6,x^7,x^8,x^9,x^10,
           x^11,x^12,x^13,x^14,x^15, x^16,x^17,x^18,x^19,x^20,
           x^21,x^22,x^23,x^24)  
b = y.noise
theta <- solve(A,b)
# Unfortunately, the machine precision is not sufficient to capture the subtle differences between columns 
# once we go to the 24th power.

# So let's make it simpler to illustrate the point:
xb <- x[1:13]
A <- cbind(rep(1,length(xb)),xb,xb^2,xb^3,xb^4,xb^5, xb^6,xb^7,xb^8, xb^9,xb^10,xb^11,xb^12)  
b = y.noise[1:13]
theta <- solve(A,b)

# We plot the model with unconstrained limits in the axis...
x.tmp <- seq(0,1,length.out = 2000)
plot(x.tmp,mle.model(x.tmp,theta[2],theta[3],theta[4],theta[5],theta[6],theta[7],theta[8],theta[9],
                     theta[10],theta[11],theta[12],theta[13],0,0,0,0,0,0,0,0,0,0,0,0,theta[1]),col="red",ty="l")
points(xb,b,pch=15,col="blue")
# And then we scale the plot so we can see what's happening
plot(x.tmp,mle.model(x.tmp,theta[2],theta[3],theta[4],theta[5],theta[6],theta[7],theta[8],theta[9],
                     theta[10],theta[11],theta[12],theta[13],0,0,0,0,0,0,0,0,0,0,0,0,theta[1]),col="red",ty="l",ylim=c(0.8,2.2))
points(xb,b,pch=15,col="blue")
# Yes! We nailed it! Did we not?

# If we had solved the linear system of 25 eqs, and since we have "nailed" the data, the lok-lik would have been 
print(dmvnorm(rep(0,25),rep(0,25),Sigma,log=T))
# This is maximum likelihood that can possibly be achieved with any model because it passes exactly through all of the points.

# Lessons learnt:
# 1 When the parameter space is large and the likelihood is multimodal, there is no guarantee 
# that the search algorithms will reach a global maximum. 
# 2 Even if we could, it may not be desirable!

# So let's use Occam's razor. Ahhhh, so the problem was solved: the second model is overly 
#complicated because it has 25 parameters.

# Is this true?
# No, it is not. Take a look at this 1 parameter model
loglik4 <- function(omega)
{
  predicted <- 1.5+sin(2*pi*x*omega)
  loglik <- -dmvnorm(y.noise,predicted,Sigma,log=T) # In reality, this is MINUS the LogLikelihood
  return(loglik)
}

mle <- mle(loglik4,list(omega=1000))
print(mle)
coef(mle)
logLik(mle)
# Oooooops! Same log-likelihood, which is the maximum attainable, with only one parameter. 
# This is definitely unbeatable. 

# Let us plot the solution at the values of x that we have measured
# Define a function as the model at the maximum likelihood estimate of the parameter.
mle.model <- function(x,omega) return(1.5+sin(2*pi*x*omega))
# And now, plot it
plot(x,y.noise,pch=15,col="orange",xlim=c(-.1,1.1),ylim=c(.8,2.2))
p <- order(x)
lines(x[p],mle.model(x,coef(mle))[p],col="blue")

# Now, let us plot it for any arbitrary point
x2 <- seq(0,1,length.out = 10000)
y2 <- mle.model(x2,1000)
lines(x2,y2,pch=".")

# Finally, lets us look at the likelihood for any value of the parameter
test <- seq(0.,2000,length.out = 21700)
mles <- rep(NA,length(test))
for (i in 1:length(test))
{
  mles[i] <- -loglik4(test[i]) # The minus is meant to counteract the minus in the definition of loglik4
}
plot(test,mles,ty="l")
plot(test,mles,ty="l",ylim=c(-10,10))
# We already knew the existence of the omega=1000 maximum.
# But... what is the other one
test[which(mles>-5)]
test[1]
mles[1]
# So the second best solution is a horizontal line... 
# Idea to be kept in mind. There is only an extremely narrow range of values that yield significant 
# likelihoods, and virtually no other omega (outside that range and the so-called aliases) is anyway 
#near acceptable.

mle.order1Poly <- mle(loglik2,list(m=0.,c=0.))
mle.order25Poly <- mle(loglik3,list(m1=0, m2=0, m3=0, m4=0,m5=0,m6=0,m7=0,m8=0,m9=0,m10=0,
                                    m11=0, m12=0, m13=0, m14=0,m15=0,m16=0,m17=0,m18=0,m19=0,m20=0,
                                    m21=0, m22=0, m23=0, m24=0,c=0.))
mle.sine <- mle

#Also...
plot(x,y.noise)
n <- length(x)
x <- runif(n)
sigma <- 0.1
eps <- rnorm(n,0,sigma)
y.noise <- 1*x+1+eps
points(x,y.noise,pch=15,col="blue")
abline(1,1)

print(-loglik2(coef(mle.order1Poly)[1],coef(mle.order1Poly)[2])) # The minus sign counteracts the minus in the function definition
print(-loglik4(coef(mle.sine)))# The minus sign counteracts the minus in the function definition

# Second idea to keep in mind:
# The expected predictive density for new data is much worse for the sinusoidal model

load("dataset.RData")
y.noise <-y
save.image("endSession2.RData")

