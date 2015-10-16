
set.seed(10)

omega <- 1000
m <- 1
c <- 1
x <- seq(0,1,length.out=100000)
y <- 1.5+sin(2*pi*x*omega)
#plot(x,y,ty="l")
filter <- abs(y - (1*x+1)) < 4.0e-1
x.1 <- x[filter]
y.1 <- y[filter]
#points(x.1,y.1,pch=15,col="orange")
hist(y.1-m*x.1-c)
x.tmp <- seq(-.5,.5,length.out = 200)
lines(x.tmp,90*dnorm(x.tmp,0,0.1))
n <- 25
counter <- 0
filter2 <- rep(FALSE,length(x.1))
norm <- dnorm(0,0,0.1)
for (i in 1:length(x.1))
{
  k <- runif(1)
  diff<-1-dnorm(y.1[i],m*x.1[i]+c,0.1)/norm
  if (k > diff) filter2[i] <- TRUE
}
sum(filter2)
x.2 <- x.1[filter2]
y.2 <- y.1[filter2]
hist(y.2-m*x.2-c)
x.tmp <- seq(-.5,.5,length.out = 200)
lines(x.tmp,400*dnorm(x.tmp,0,0.1))

n <- 25
indices <-sample(1:length(x.2),n)
x.3 <- x.2[indices]
y.3 <- y.2[indices]
plot(x.3,y.3,pch=15,col="orange")
x <- x.3
y <- y.3
save(x,y,file="dataset.RData")
