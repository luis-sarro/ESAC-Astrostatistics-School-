setwd("~/Escritorio/ESAC-Astrostats/2015/session1")

######################################################################
# Block 1
# The ratio distribution
######################################################################
# What would be the probability density of a random variable
# defined as the quotient of two measured fluxes with 
# Gaussian uncertainties?
# Let us generate data and look at the empirical PDF.
######################################################################

mu.x <-10 # The mean flux in line 1
mu.y <-20 # The mean flux in line 2
sigma.x <- 1 # The uncertainty in the measurement of flux 1
sigma.y <- 4 # The uncertainty in the measurement of flux 2
X <- rnorm(1000000,mu.x,sigma.x) # Generate 1 million measurements
Y <- rnorm(1000000,mu.y,sigma.y) # of each
Z <- X/Y # And compute the resulting ratio

hist(Z,100,col="orange") # Draw a histogram
dens.Z <- density(Z,n=1000) # Calculae the kernel based estimate of the PDF
plot(dens.Z$x, dens.Z$y, ty="l") # Draw the kernel based estimate of the PDF
#rug(Z) # Overplot the data points

# The analytical expression for the ratio distribution is defined in terms of the following auxiliary functions
erf <- function(z) 2 * pnorm(z * sqrt(2)) - 1
a <- function(z,sigma.x,sigma.y) sqrt((z^2/sigma.x^2)+(1/sigma.y^2))
b <- function(z,mu.x,mu.y,sigma.x,sigma.y) (mu.x*z/sigma.x^2)+(mu.y/sigma.y^2) 
c2 <- function(mu.x,mu.y,sigma.x,sigma.y) (mu.x^2/sigma.x^2)+(mu.y^2/sigma.y^2) 
d <- function(z,mu.x,mu.y,sigma.x,sigma.y) exp(
  (b(z,mu.x,mu.y,sigma.x,sigma.y)^2-(c2(mu.x,mu.y,sigma.x,sigma.y)*a(z,sigma.x,sigma.y)^2))/(2*a(z,sigma.x,sigma.y)^2))
# as...
ratioOfNormals <- function(z,mu.x,mu.y,sigma.x,sigma.y)
{
  p.z <- 
    (b(z,mu.x,mu.y,sigma.x,sigma.y)*d(z,mu.x,mu.y,sigma.x,sigma.y)/a(z,sigma.x,sigma.y)^3)*
    (1/((sqrt(2*pi)*sigma.x*sigma.y)))*
    (erf(b(z,mu.x,mu.y,sigma.x,sigma.y)/a(z,sigma.x,sigma.y))-erf(-b(z,mu.x,mu.y,sigma.x,sigma.y)/a(z,sigma.x,sigma.y)))
    +(1/((a(z,sigma.x,sigma.y)^2)*pi*sigma.x*sigma.y))*exp(-c2(mu.x,mu.y,sigma.x,sigma.y)/2)
  return(p.z)
  }

# So now, let us overplot this analytical expression
z <- seq(min(Z),max(Z),length.out=1000) # First we generate the x values
pdf <- ratioOfNormals(z,mu.x,mu.y,sigma.x,sigma.y)# Then, we evaluate the analytical expression 

# We normalise it...
require(caTools) # We need this library to approximate integrals as the sum of the area of rectangles 
area <- trapz(z,pdf)
# And finally we overplot:
plot(dens.Z$x,dens.Z$y)
lines(z,pdf/area,lwd=3,col="blue")

########################################################################
# What if the normal distributed random variables are centred around 0?
# The Cauchy distribution
# https://en.wikipedia.org/wiki/Cauchy_distribution

mu.x <-0
mu.y <-0
sigma.x <- 1
sigma.y <- 4
X <- rnorm(10000,mu.x,sigma.x)
Y <- rnorm(10000,mu.y,sigma.y)
Z <- X/Y

dens.Z <- density(Z,n=1000)
print(trapz(dens.Z$x,dens.Z$y)) # 1000 values is obviously not enough. Not even 10000000 would be.
plot(dens.Z$x,dens.Z$y,ty="l")
rug(Z)
range(Z)

plot(density(Z,from=-5,to=5),xlim=c(-5,5))
x <- seq(-5,5,length.out=200)
gamma <- sigma.x/sigma.y
x0 <- 0
y <- (1/pi)*(gamma/((x-x0)^2+gamma^2))

cauchy <- dcauchy(x, location = 0, scale = sigma.x/sigma.y, log = FALSE)
lines(x,cauchy,col="blue")

######################################################################
# Block 2
# The generalized chi-squared distribution
######################################################################
# What would be the probability density of a random variable
# defined as the sum of the squares of several gaussian distributed 
# random variables?
# Let us generate data and look at the empirical PDF.
######################################################################
mu.x <-10 # Suppose x and y are velocities derived from proper motions 
mu.y <-10 # and distances (and hence have large --*non-gaussian*, but let us overlook this-- errors)
mu.z <- 10 # and the third velocity is obtained as a radial velocity
sigma.xy <- matrix(c(30,5,5,30),byrow=T,nrow=2) # Correlated errors
sigma.z <- 20
require("mvtnorm") # We need this to deal with multivariate normal PDFs
propmot <- rmvnorm(1000000,c(mu.x,mu.y),sigma.xy)
Z <- rnorm(1000000,mu.z,sigma.z)

V2 <- (propmot[,1]^2+propmot[,2]^2+Z^2)
dens.V2 <- density(V2)
plot(dens.V2$x,dens.V2$y,ty="l")

V <- sqrt(V2)
dens.V <- density(V)
plot(dens.V$x,dens.V$y,ylim=c(0,0.07),ty="l")
lines(seq(0,120,length.out = 200),dnorm(seq(0,150,length.out = 200),25,7),col="blue",lwd=2)

# This is the generalized chi-squared distribution
# We can sample from the (possibly non-central) chi-squared distributions in R with
ranvar <- rchisq(100000,3,0)
dens.ranvar <- density(ranvar)
plot(dens.ranvar$x,dens.ranvar$y,ty="l")
rug(ranvar)
# Or get the pfg values with 
x <- seq(0,25,length.out = 200)
y <- dchisq(x,3,0)
lines(x,y,lwd=2,col="blue")

######################################################################
# Block 3
# Mistery
######################################################################
# What would be the probability density of the difference between the 
# true mean of a normal random variable and an estimated mean from a
# small sample.
# Let us generate data and look at the empirical PDF.
######################################################################

n <- 5 # Let us do 10^4 experiments, each with 5 measurements
mu <- 5 # Define the mean as 5
sd <- 2 # And the standard deviation as 2 
t <- NULL
for (i in 1:10^4)
{
x <- rnorm(n,mu,sd) # Generate 5 measurements from out gaussian
norm.diff <- sqrt(n)*(mean(x)-mu)/sd(x) # Calculate the difference from the true mean (scaled by the std. deviation and sqrt(5))
t <- c(t,norm.diff) # append it to the vector t
}
plot(density(t))
rug(t)
plot(density(t),xlim=c(-5,5))

# # Student's t
# https://en.wikipedia.org/wiki/Student's_t-distribution
x <- seq(-10,10,length.out=100)
y <- dt(x,df=n-1)
lines(x,y,col="blue",lwd=3)

# Snedecor F
# https://en.wikipedia.org/wiki/F-distribution

f <- rf(10000,4,6)
plot(density(f),ty="l")
rug(f)
y <- df(x,4,6)
lines(x,y,col="blue",lwd=3)

# Generating random numbers from a giver distribution

# 1 For common R distributions...
# The normal distribution
x <- seq(-5,5,length.out = 100)
data <- rnorm(1000000,0,1) # play with n and check how many samples we need before it really looks gaussian
plot(density(data))
lines(x,dnorm(x,0,1),col="blue")
#The uniform distribution
x <- seq(0,1,length.out = 100)
data <- runif(1000000,0,1) # play with n and check how many samples we need before it really looks flat
h <-hist(data,50)
# The chi-squared distribution
data <- rchisq(1000000,1) # play with n and check how many samples we need before it really looks chi-squared
x <- seq(0,max(data),length.out = 200000)
# the bw computed by density as default  is too broad
plot(density(data,from=-.1,to=5,n=2^15),xlim=c(0,5),ylim=c(0,40))
rug(data[1:100])
lines(x,dchisq(x,1),col="blue")
plot(density(data,from=-.1,to=5,bw=0.0001,n=2^15),xlim=c(0,5),ylim=c(0,40))
lines(x,dchisq(x,1),col="blue")
plot(density(data,from=-.1,to=5,bw=0.0001,n=2^15),xlim=c(0,5),ylim=c(0,0.1))
lines(x,dchisq(x,1),col="blue")
plot(density(data,from=-.1,to=5,bw=0.0001,n=2^15),xlim=c(0,0.02),ylim=c(0,40))
lines(x,dchisq(x,1),col="blue")
# Student's t distribution
x <- seq(-5,5,length.out = 100)
data <- rt(1000000,3) 
plot(density(data,n=100000),xlim=c(-5,5))
lines(x,dt(x,3),col="blue")

# 2 Using Cumulative distribution functions
# Let us sample from the uniform...
# We already know how to: use runif
# But for more complex, maybe non-analytical PDFs...

# What is a CDF?

# Let us draw a Gaussian PDF
x <- seq(-5,5,length.out = 300)
y <- dnorm(x)
plot(x,y,ty="l",ylim=c(0,1))

# Now, let us calculate the area under the curve up to a given xval
xval <- -2
p <- c(xval,dnorm(xval))
points(p[1],p[2],pch=16,col="orange")
p <- rbind(p,c(xval,0))
lines(p,col="orange",lwd=3)
filter <- x < xval
lines(x[filter],dnorm(x)[filter],ty="h")
area <- trapz(x[filter],dnorm(x)[filter])
p <- c(xval,area)
points(p[1],p[2],pch=16,col="red")
# We can do this for several xvals

# Well, the CDF is the curve  that represents this area as a function of xval 
y2 <- seq(0.000001,0.999999,length.out = 300)
x2 <- qnorm(y2)
lines(x2,y2,lwd=3,col="blue")

# The algorithm to generate random numbers from a distributions proceeds
# by generating values yval from a uniform distribution between 0 and 1, and then 
# adding to the sample the xval that connects to that yval under the inverse cdf

# Let us try with an example: the uniform distribution between 0 and 10
x <- seq(-1,11,length.out = 10000)
pdf <- dunif(x,0,10)
plot(x,pdf,ty="l")
# Step 1: Generate values of the CDF
cdf <- rep(NA,length(pdf))
for (i in 1:length(x)) cdf[i] <- trapz(x[1:i],pdf[1:i])
plot(x,cdf, ty="l")
lines(x,pdf,col="blue")

# Now, let us interpolate the inverse CDF
# First, let us look to the inverse CDF
y <- x
x <- cdf
plot(x,y,ty="l") # Ooops. This has vertical asymptotes
# Now we interpolate the inverse CDF
f <- splinefun(x,y)

# Finally, generate a random sample of size n from the uniform distribution between 0 and 1
n <- 100000
yvals <- runif(n,0,1)
# ... and evaluate the inverse CDF at these points
xvals <- f(yvals)

hist(xvals)

# Can you try with normal distribution?
