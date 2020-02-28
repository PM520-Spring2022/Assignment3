# example taken from http://www.mas.ncl.ac.uk/~ndjw1/teaching/sim/metrop/indep.r

# metropolis-hastings independence sampler for a
# gamma based on normal candidates with the same mean and variance
set.seed(37)
gamm<-function (n, a, b) 
{
  mu <- a/b   # the mean of the gamma distribution
  sig <- sqrt(a/(b * b))   # the stadard deviation of the gamma distn
  vec <- vector("numeric", n)   # this is where we are going to put the random variables we generate
  x <- a/b 
  vec[1] <- x # We arbitrarily start the MCMC process at the mean
  for (i in 2:n) {
    can <- rnorm(1, mu, sig)
    hprob <- min(1, (dgamma(can, a, b)/dgamma(x,a,b))/(dnorm(can, mu, sig)/dnorm(x, mu, sig)))   # where is the q term here?
    u <- runif(1)
    if (u < hprob) 
      x <- can
    vec[i] <- x
  }
  return (vec)
}


vec<-gamm(10000,2.3,2.7)
#vec<-gamm(10000,0.1,0.01)
par(mfrow=c(1,1))
plot(ts(vec))
plot(ts(vec)[1:1000],type='l')
acf(vec)
hist(vec[1000:10000],30)
curve(dgamma(x,2.3,2.7),from=0,to=3)
#curve(dgamma(x,0.1,0.01))
par(mfrow=c(1,1))
# end

# Try out this R code for different shape and scale parameters. In particular, try it for a=0.1, b=0.01. 
# Notice how there is a bit of a problem with the sampler getting "stuck" at very small values
# values. (To see what it should look like use this command: curve(dgamma(x,0.1,0.01)) ).
# Modify the code to keep track of acceptance probabilities and plot the acf
# When does the sampling scheme do worst? 
# How might you modify this sampling scheme to be more efficient?
