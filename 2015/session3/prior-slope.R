m=1
c=0

# sample from flat (but proper) prior
m.min=0
m.max=1.0e3
pi = 1/(m.max-m.min)

cdf.prior <- function(m)
{
  cdf <-m/m.max
  return(cdf)
}

icdf.prior <- function(cdf)
{
  m <-cdf*m.max
  return(m)
}

cdf.eq <- seq(0,1,length.out=1000) 
m.eq <- icdf.prior(cdf.eq)
plot(m.eq,cdf.eq,ty="l")

cdf.m <- runif(10000,0,1)
m <- icdf.prior(cdf.m)

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
  lines(x[i,],y[i,],col=rgb(0,0,0,0.5))
}



