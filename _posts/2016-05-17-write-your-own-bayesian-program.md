---
layout: post
title:  "Write your own Bayesian MCMC phylogenetic program with R"
---
This is an R tutorial I wrote for a workshop on Bayesian phylogenetics at Bristol University in 2016. The tutorial first introduces a Bayesian problem where we calculate, analytically, the posterior distribution of the molecular branch length between two species. Then the tutorial moves onto calculating the posterior by MCMC sampling. Here you can learn MCMC concepts such as proposal densities, proposal window size, mixing, convergence, efficiency, autocorrelation, and effective sample size. By tweaking the parameters you can see how the efficiency of the MCMC and the parameter’s ESS are affected. At the end there is an exercise where you can modify the MCMC program to convert it to a two-parameter model where you estimate the divergence time and the molecular evolutionary rate between the two species. Send me an email if you would like to see the answers to the exercises. The workshop solves example 7.1 (p. 216) from Ziheng Yang’s book “Molecular Evolution: A Statistical Approach” by OUP (see also example 6.4, p.188, and problems 7.1 and 7.5, p.260 in the book).

```r
# Workshop "Bayesian methods to estimate species divergence times"
# July 2015
# Mario dos Reis, Ziheng Yang and Phil Donoghue

# A Bayesian phylogenetics inference program

rm(list=ls())
# ##########################################
# Simple example, calculated analytically
# ##########################################

# The data are a pairwise sequence alignment of the 12s RNA gene sequences
# from human and orangutan.
# * The alignment is n=948 nucleotides long.
# * We observed x=90 differences between the two sequences.
# (Example from p.216, Molecular Evoltion: A Statistical Approach)

# Our aim is to estimate the molecular distance d, in substitutions per site,
# between the two sequences. Because the two species are closely related (i.e.,
# the alignment has less than 10% divergence) we use the Jukes and Cantor (1969)
# nucleotide substitution model (JC69).

# The Bayesian estimate of d, is given by the posterior distribution:
# p(d | x) = C * p(d) * p(x | d).

# p(d) is the prior on the distance.
# p(x | d) is the (Jukes-Cantor) likelihood, i.e. the probability of observing
# the data x given d, n. C is a normalising constant, C = 1 / ∫ p(d) p(x | d) dd

# likelihood under JC69 model, i.e. L = p(x | d):
L <- function(d, x, n) (3/4 - 3*exp(-4*d/3)/4)^x * (1/4 + 3*exp(-4*d/3)/4)^(n - x)

# plot the likelihood:
curve(L(x, 90, 948), from=0, to=0.2, n=500, xlab="distance, d", ylab="likelihood, p(x | d)")

# The maximum likelihood estimate is (the highest point in the curve)
d.mle <- -3/4 * log(1 - 4/3 * 90/948) # [1] 0.1015060
abline(v=d.mle, lty=2)
title(expression(hat(d) == 0.1015))

# For the prior on d, we use an exponential distribution with mean = 0.2
# p(d) = exp(-d/mu) / mu; mu = 0.2.
d.prior <- function(d, mu) exp(-d/mu) / mu

# un-normalised posterior:
d.upost <- function(d, x, n, mu) d.prior(d, mu=mu) * L(d, x=x, n=n)

# To calculate the normalising constant C, we need to integrate the un-normalised posterior:
integrate(d.upost, lower=0, upper=Inf, x=90, n=948, mu=0.2, abs.tol=0)
# 5.167762e-131 with absolute error < 1.9e-135
C  <- 1 / 5.167762e-131

# The posterior distribution function (when x=90, n=948 and mu=0.2)is thus
d.post <- function(d) C * d.upost(d=d, x=90, n=948, mu=0.2)

# Plot posterior distribution:
curve(d.post(x), from=0, to=0.2, n=500, xlab="d", ylab="density")
# Add prior:
curve(d.prior(x, mu=0.2), add=TRUE, lty=2)

# EXERCISE 1: Try adding the likelihood curve to the plots above.
# HINT: You may want to scale the likelihood by a constant like C * 3

# The posterior mean of d is found by integrating ∫ d * p(d | x) dd
f <- function(d) d * d.post(d)
integrate(f, lower=0, upper=Inf, abs.tol=0)
# 0.1021246 with absolute error < 5.5e-06
# Add it to the plot:
abline(v=0.1021246, lty=3) # posterior mean

# Add MLE:
abline(v=d.mle, lty=3, col="red")

# The posterior distribution can be used to answer questions such as
# what is the probability that d > 0.1?
integrate(d.post, lower=0.1, upper=Inf)
# 0.5632372 with absolute error < 1.6e-08, i.e. 56.3%

# The 95% equal-tail credibility interval is (0.08191, 0.12463)
integrate(d.post, lower=0, upper=0.08191) # 2.5%
integrate(d.post, lower=0, upper=0.12463) # 97.5%

# EXERCISE 2 (VERY HARD): Can you work out how to calculate the 99% CI?
# HINT: You need to find the upper bounds where the integral is 0.005
# and 0.995 respectively. The integrate and optim functions in R may be useful.


# #############################################
# The same example analysed using
# Markov chain Monte Carlo (MCMC)
# #############################################

# We now obtain the posterior distribution by MCMC sampling.
# In most practical problems, constant C cannot be calculated (either
# analytically or numerically), and so the MCMC algorithm becomes necessary.

# Draft MCMC algorithm:
# 1. Set initial state for d.
# 2. Propose a new state d* (from an appropriate proposal density).
# 3. Accept or reject the proposal with probability
#      min(1, p(d*)p(x|d*) / p(d)p(x|d))
#      If the proposal is accepted set d=d*, otherwise d=d.
# 4. Save d.
# 5. Go to step 2.
d.mcmcf <- function(init.d, N, w) {
   # init.d is the initial state
   # N is the number of 'generations' the algorithm is run for.
   # w is the 'width' of the proposal density.
   	d <- numeric(N+1) # Here we will keep the visited values of d.
   	d[1] <- init.d
   	acc.prop <- 0 # acceptance proportion
   	for (i in 1:N) {
        # here we use a uniform density with reflection to propose a new d*
   		d.prop <- abs(d[i] + runif(1, -w/2, w/2))
   		#d.prop  <- abs(d[i] + rnorm(1, 0, w)) # a normal proposal
   		p.ratio <- d.upost(d.prop, x=90, n=948, mu=0.2) /
                  d.upost(d[i], x=90, n=948, mu=0.2)
   		alpha <- min(c(1, p.ratio))
                # if ru < alpha accept proposal:
   		if (runif(1, 0, 1) < alpha) { d[i+1] <- d.prop; acc.prop  <- acc.prop + 1 }
                # else reject it:
   		else d[i+1] <- d[i]
   		# Can you think of a way to print the progress of the MCMC to the screen?
   	}
    # print out the proportion of times the proposal was accepted
   	print(c("Acceptance proportion:", acc.prop/N))
   	return (d) # return vector of d visited during MCMC
}

# Test the algorithm
d.1 <- d.mcmcf(0.5, 100, 0.08)

# plot output:
plot(d.1, type='l', xlab="generation", ylab="d") # This is the 'trace' plot of d.
# the chain wiggles about spending time at different values of d in
# proportion to their posterior density.

# Add 95% CI as reference:
abline(h=c(0.08191, 0.12463), lty=2, col="red")
#text(locator(1), "Burn-in phase")
#text(locator(1), "Stationary phase")

# It may be useful to caculate the time it takes to run x generations:
system.time(d.mcmcf(0.5, 1e4, 0.08))
# ~ 0.6s in my laptop, that means 1e6 would take ~60s or 1min.

# Get a very long sample:
d.l <- d.mcmcf(0.5, 1e4, 0.08)
plot(d.l, type='l')
# The trace plot looks like a dense, hairy caterpillar!
# This is good!

# Remove initial state and burn-in
d.l <- d.l[-(1:101)]

hist(d.l, prob=TRUE, n=50, xlab="d", main="") # The sample histogram from the MCMC
curve(d.post(x), add=TRUE, col="red", lwd=2)  # The true posterior distribution.

# Alternatively:
plot(density(d.l, adj=1.5)); rug(d.l)
curve(d.post(x), add=TRUE, col="red", lwd=2)

# EXERCISE 3: What does adj in the density function do? Try different values
# such as adj=0.1 or adj=1.0.
# How long (how many generations, N) do you need to run the MCMC so that
# the density plot is almost identical to the true distribution?

# Now we can get the posterior mean, and 95% CI from the sample:
mean(d.l) # should be close to 0.1021246
quantile(d.l, prob=c(0.025, 0.975)) # should be close to (0.08191, 0.12463)

# Effect of step size, w:
d.1 <- d.mcmcf(0.5, 100, w=0.08) # moderate step size
plot(d.1, type='b', col="black", pch=19, cex=.5) # chain moves about

d.2 <- d.mcmcf(0.5, 100, w=0.01) # tiny step size
lines(d.2, type='b', col="green", pch=19, cex=.5) # baby steps, chain moves slowly

d.3 <- d.mcmcf(0.5, 100, w=1) # large step size
lines(d.3, type='b', col="blue", pch=19, cex=.5) # chain gets stuck easily

# EXERCISE 4: The optimal acceptance proportion is about 30%. Try various step sizes, w,
# and see the effect on the acceptance proportion. What is the optimal w?
# Finding the best step size is called fine-tuning.


# ##########################################
# Autocorrelation, efficiency, effective
# sample size and thinning:
# ##########################################

# States in a MCMC sample are autocorrelated, because each new state is
# either the previous state or a modification of it.

x <- d.l[-length(d.l)] # remove the last state
y <- d.l[-1] # remove the first state
plot(x, y, xlab="current state", ylab="next state", main="lag, k=1")
cor(x, y) # this is the autocorrelation for the d sample for lag 1, ~50-60%

# We can use the acf function to calculate the autocorrelation for various lags:
d.acf <- acf(d.l)  # The correlation is lost after a lag of k=10
d.acf$acf # the actual acf values are here.

# The efficiency of an MCMC chain is closely related to the autocorrelation.
# Intuitively, if the autocorrelation is high, the chain will be inefficient,
# i.e. we will need to run the chain for a long time to obtain a good
# approximation to the posterior distribution.
# The efficiency of a chain is defined as:
# eff = 1 / (1 + 2(r1 + r2 + r3 + ...))
# where ri is the correlation for lag k=i.

eff <- function(acf) 1 / (1 + 2 * sum(acf$acf[-1]))
eff(d.acf)  # ~ 30% What does that mean?

# Efficiency is closely related to the effective sample size (ESS),
# which is defined as: ess = N * eff
N <- length(d.l) # 9,900
length(d.l) * eff(d.acf) # ~ 2,700
# That is, our MCMC sample of size 9,900 has as much information about d as an
# independent sample of size ~2,700. In other words, 30% efficiency means
# the chain is 30% as efficient as an independent sample.

# EXERCISE 5: Run a long chain (say N=1e4) with step size w=0.01.
# Calculate the efficiency and the ESS.
# Do it at home: Try various values of w. What value of w gives the highest ESS?
#                Plot estimated efficiency vs. w.

# Usually, to save hard drive space, the MCMC is 'thinned', i.e. rather than
# saving every single value sampled, we only save every m-th value. The resulting
# sample will be smaller (thus saves hard drive space) but still has nearly as much info
# as the total sample. In our case, we may reduce the sample to ~1/4 the original size:
d.reduced <- d.l[seq(from=1, to=N, by=4)]
mean(d.reduced) # still close to 0.1021246
quantile(d.reduced, c(0.025, 0.975)) # still close to (0.08191, 0.12463)


# ################################################
# MCMC diagnostics.
# How do we know our results are reliable?
# ################################################

# (1) Run the chain from different starting values, the chains
#     should converge to the same region in parameter space:
d.high  <- d.mcmcf(0.4, 1e2, 0.1)
plot(d.high, col="red", ty='l', ylim=c(0, 0.4))

d.medium <- d.mcmcf(0.1, 1e2, 0.1)
lines(d.medium, col="black", ty='l')

d.low <- d.mcmcf(0, 1e2, 0.1)
lines(d.low, col="blue", ty='l')

abline(h=c(0.08191, 0.12463), lty=2)

# (2) Compare the histograms ,and summary statistics obtained
#     from several long, independent runs:
d.l1 <- d.mcmcf(0.4, 1e3, 0.1)[-(1:101)]
d.l2 <- d.mcmcf(0.0, 1e3, 0.1)[-(1:101)]

plot(density(d.l1, adj=1), col="red")
lines(density(d.l2, adj=1), col="blue")

mean(d.l1); mean(d.l2)
# relative difference
abs(mean(d.l1) - mean(d.l2)) / mean(c(d.l1, d.l2)) * 100

# I cannot stressed enough how important it is that you run all your
# MCMC analyses at least twice!
# You cannot assess the robustness of a MCMC if it has been run just
# once!

# Try again using N=1e5 (will take a little to finish)
d.l1 <- d.mcmcf(0.4, 1e5, 0.1)[-(1:101)]
d.l2 <- d.mcmcf(0.4, 1e5, 0.1)[-(1:101)]

# (3) Calculate ESS. Ideally, we should aim for ESS > 1,000, when
#     summary statistics such as the posterior mean and 95% CI can be
#     accurately estimated. ESS > 10,000 offers superb performance, but this
#     is rarely achieved in phylogenetic problems. ESS < 100 is poor:
length(d.l1) * eff(acf(d.l1))  # > 22,000 (nice!)

# (4) Plot running means:
rm1 <- cumsum(d.l1) / (1:length(d.l1))
rm2 <- cumsum(d.l2) / (1:length(d.l2))
plot(rm1, ty='l', ylab="posterior mean", xlab="generation"); lines(rm2, col="red")

# There are more diagnostics tests and statistics that can be applied, for
# example the potential scale reduction statistic (p. 242, Molecular Evolution:
# A statistical).

# R packages such as CODA have been written to perform MCMC
# diagnostics, and they may prove quite useful.

# ################################################
# Bayesian MCMC - Molecular dating
# ################################################

# Now consider estimating the time of divergence between human and orangutan
# using the data above, x=90, n=948.
# The molecular distance is the product of geological time, t, and the
# molecular rate, r: d = r * t.

# Our JC69 likelihood based on r and t is:
L.rt <- function(r, t, x, n) L(d=r*t, x=x, n=n)

# The likelihood is now a 2D surface:
r <- seq(from=2e-3, to=6e-3, len=50)
t <- seq(from=10, to=40, len=50)

# creates a table of all combinations of r and t values:
rt <- expand.grid(r=r, t=t)

z <- L.rt(r=rt$r, t=rt$t, x=90, n=948)

# The likelihood surface
image(r, t, matrix(z, ncol=50), xlab="r (s/My)", ylab="t (Ma)")
contour(r, t, matrix(z, ncol=50), nlevels=5, add=TRUE)

# The likelihood looks like a flat, curve mountanous ridge
# The top of the ridge is located at d.mle = 0.1015060 = r * t
curve(d.mle/x, add=TRUE, col="blue", lwd=2, lty=2)

# The molecular data are informative about d only, but not about r or t
# separately. We say that r and t are confounded parameters. Any combination
# of r and t that satisfy 0.1015060 = r * t are ML solutions.

# We can make a neat 3D plot of the likelihood:
persp(r, t, matrix(z, ncol=50), zlab="likelihood", theta=40, phi=40)

# The Bayesian method provides a useful methodology to estimate t and r through
# the use priors.
# The fossil record suggests that the common ancestor of human-orang lived
# 33.7-11.2 Ma (see Benton et al. 2009 in the Timetree of Life book)
# Thus the prior on t is a normal distribution with mean 22.45 (the midpoint
# of the fossil calibraiton) and sd 5.6, so that the 95% range is roughly
# 33.7 to 11.2.
t.prior <- function(t, mu=22.45, sd=5.6) dnorm(t, 22.45, 5.6)
curve(t.prior(x), from=5, to=40, n=5e2, xlab="t (Ma)", ylab="density")

# For the rate, r, we use an exponential distribution with mean
# mu = 0.10 / 22.45 = 0.00445 (based on the distance MLE and the midpoint
# of the calibration)
r.prior <- function(r, mu=0.10/22.45) dexp(r, 1/mu)
curve(r.prior(x), from=0, to=.02, xlab="r (s/My)", ylab="density")

# We can plot the joint prior:
z1 <- r.prior(rt$r) * t.prior(rt$t)
image(r, t, matrix(z1, ncol=50), xlab="r (s/My)", ylab="t (Ma)")
contour(r, t, matrix(z1, ncol=50), nlevels=5, add=TRUE)

# The unnormalized posterior is thus
rt.upost <- function(r, t, x=90, n=948) {
  r.prior(r) * t.prior(t) * L.rt(r, t, x, n)
}

# We can plot the un-normalised posterior:
z2 <- rt.upost(rt$r, rt$t)
image(r, t, matrix(z2, ncol=50), xlab="r (s/My)", ylab="t (Ma)")
contour(r, t, matrix(z2, ncol=50), nlevels=5, add=TRUE)
curve(d.mle/x, add=TRUE, col="blue", lty=2, lwd=2)

# EXERCISE 6: Can you plot the un-normalise posterior as a 3D surface using persp?

# We are now ready to build a two-parameter MCMC!

# The draft algorithm is as follows:
# 1. Set initial states for t and r.
# 2. Propose a new state r* (from an appropriate proposal density).
# 3. Accept or reject the proposal with probability
#      min(1, p(r*)p(t)p(x|d*) / p(r)p(t)p(x|d))
#      If the proposal is accepted set r=r*, otherwise r=r
# 4. Repeat 2 and 3 for t (i.e. propose t*, accept or reject).
# 5. Save r and t.
# 6. Go to step 2.

# EXERCISE 7: Create a function called rt.mcmcf based on the d.mcmcf function.
# Run your newly written phylogenetics MCMC program! You will need two step sizes,
# w1 and w2, for r and t. Vary the step sizes and check the performance of the
# MCMC chain. Run the analysis several times, and calculate posterior means and
# 95% CIs for r and t. Calculate efficiency and ess. A 3D density histogram of the
# MCMC posterior can be obtained using the kde2d in the MASS package.

# EXERCISE 8: Try using a normal density (rather than uniform) with reflection as
# the proposal distribution for r and t. Does it make a difference?

# EXERCISE (VERY HARD) 9: The data are now a longer alignment with x=900 and n=9480.
# When you have so much data, the likelihood is soo tiny that the computer has
# trouble calculating it. In this case it is best to work with the log-likelihood:
# lL.d <- function(d, x, n)  x * log(3/4 - 3*exp(-4*d/3)/4) +
#                            (n - x) * log(1/4 + 3*exp(-4*d/3)/4)
# You then need to calculate the log-posterior:
# log p(d | x) = log p(d) + log p(x | d)

# An MCMC algorithm can be constructed that works only with the log-prior and
# the log-likelihood. The proposal ratio is then
# exp (log p(d*) + log p(x | d*) - log p(d) - log p(x | d))

# Can you write an MCMC program to calculate the posterior of d for such large
# dataset? All Bayesian phylogenetic software work with the log-prior and log-likelihood.

# EXERCISE (VERY HARD) 10: The normalising constant in this case is a double integral
# C = 1 /  ∫∫ p(r) p(t) p(x | d) dr dt
# Can you calculate C using R's cubature package (the normal integrate function
# does not work, you need special multidimensional integration packages)?
# The cubature package is not a default package, so you need to use install.packages
# to get it.
```
