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


vec1<-gamm(10000,2.3,2.7)
vec2<-gamm(10000,2.3,2.7)
plot(ts(vec1))
plot(ts(vec1)[1:1000],type='l')
acf(vec1)
hist(vec1[1000:10000],30)
plot(ts(vec2))
plot(ts(vec2)[1:1000],type='l')
acf(vec2)
hist(vec2[1000:10000],30)
curve(dgamma(x,2.3,2.7),from=0,to=3)

# convert to mcmc objects, with a burn-in 
MCMC1<-mcmc(vec1,start=1000)
MCMC2<-mcmc(vec2,start=1000)

# combine different mcmc chain objects to an mcmc list.
Combined<-mcmc.list(list(MCMC1,MCMC2))

# gelman functions are 
gelman.plot(Combined) # for plots
print(gelman.diag(Combined)) # for diagnostic values




vec3<-gamm(10000,0.1,0.01)
vec4<-gamm(10000,0.1,0.01)
plot(ts(vec3))
plot(ts(vec3)[1:1000],type='l')
acf(vec3)
hist(vec3[1000:10000],30)
plot(ts(vec4))
plot(ts(vec4)[1:1000],type='l')
acf(vec4)
hist(vec4[1000:10000],30)
curve(dgamma(x,0.1,0.01))

MCMC3<-mcmc(vec3,start=1000)
MCMC4<-mcmc(vec4,start=1000)

# combine different mcmc chain objects to an mcmc list.
Combined2<-mcmc.list(list(MCMC3,MCMC4))

# gelman functions are 
gelman.plot(Combined2) # for plots
print(gelman.diag(Combined2)) # for diagnostic values



# Try out this R code for different shape and scale parameters. In particular, try it for a=0.1, b=0.01. 
# Notice how there is a bit of a problem with the sampler getting "stuck" at very small values
# values. (To see what it should look like use this command: curve(dgamma(x,0.1,0.01)) ).
# Modify the code to keep track of acceptance probabilities and plot the acf
# When does the sampling scheme do worst? 
# How might you modify this sampling scheme to be more efficient?
