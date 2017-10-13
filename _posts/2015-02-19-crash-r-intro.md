---
layout: post
title:   "A crash introduction to R"
author: "Mario dos Reis"
---

I wrote this R tutorial several years ago for people working in my lab, and then completely forgot about it. Apparently, the tutorial had been doing the rounds, and was found to be quite useful, so I have decided to post it here. It is a very brief tutorial meant for people who have used R a little bit and that consider themselves R beginners. Lines that start with a # are comments and are ignored by R. Lines in black are R commands. You should cut and paste them into the R console and see their effect.

```r
# --- Crash introduction to R ---
#
# Mario dos Reis, June 2008

# This tutorial was prepared with R v 2.5.1 on Linux
# Updated Feb 2015, Tested with R v 3.0.2 on Ubuntu Linux 14.04.2

# First, R can be used as a simple calculator
2 + 5

2 * 10 / pi

# `pi' is a built in variable, and you can probably guess its value
pi # [1] 3.141593
2 * pi^2 # [1] 19.73921
exp(1) # Exponential function
log(10) # Natural logarithm
log(0) # [1] -Inf
2e4 # Exponential notation: 2x10^4

# The built in variable `Last.value' stores the last value computed
5^2 # [1] 25
sqrt(.Last.value) # [1] 5

# Vectors are the simplest data structure in R
w <- 1:10 # this generates a vector with all the numbers from 1 to 10
x <- rnorm(1000) # this generates 1,000 pseudo-random numbers from a normal distribution
y <- rnorm(1000, mean=3, sd=1)

# `x' and `y' are numerical vectors of length 1,000
is.vector(x) # [1] TRUE
length(x) # [1] 1000

# A random cloud of points
plot(x, y)

# A typical boxplot
boxplot(x, y)

# or ... (can you spot the difference?)
boxplot(x, y, names=c("x", "y"), main="A Boxplot")

# vectors can be easily concatenated
z <- c(x, y)
hist(z, n=50) # a histogram

# Performing basic statistics is straightforward
mean(x)
sd(x) # standard deviation
summary(x) # various summary statistics
cor(x, y) # correlation between `x' and `y'

# R can easily handle matrices
# This generates a 10x10 matrix of pseudo-random Poisson numbers
x <- matrix(rpois(100, lambda=5), ncol=10)
# summary statistics still apply
mean(x)
range(x) # maximum and minimum values in `x'

x[1,] # the first row of the matrix

# EXERCISE: can you output the 5th column of the matrix?

# How about calculating the mean for each row in `x'?
# R has control flow structures that alow to do this easily:
means <- numeric(10) # this creates a vector of zeroes of length 10
for (i in 1:nrow(x)) means[i] <- mean(x[i,])

# EXERCISE: Try using functions like `sd' or `range'

# R has several built in datasets for didactic purposes
# One of such datasets is `women'
women

# Women is a data frame
is.data.frame(women) # [1] TRUE

# Data frames are special types of `lists', which are more general
# data structures than vectors or matrices
names(women)
women$height # this is a vector of heights
women$weight # this is a vector of weights

# EXERCISE: try the mean() and summary() functions on the women dataset

# A simple regression model
plot(women$height, women$weight)

# This looks like a linear relationship, is it?
women.lm <- lm(weight ~ height, data=women)
# Function `lm' fits a linear model of weight on height
# The term `weight ~ height' is called a formula, this is the
# shorthand for the statistical model
# `weight = alpha + beta * height + error ~ N(0, sigma)'

# adds the estimated regression line to the plot
abline(women.lm, col="red")

# `women.lm' is an R object that contains the output of the `lm'
# function. `women.lm' is also a list.
summary(women.lm)
anova(women.lm) # The classical ANOVA table
names(women.lm)

# Individual items within `women.lm' can be extracted with the $ sign
women.lm$residuals # the residuals obtained after fitting the model

plot(women.lm$residuals ~ women$height)
abline(h=0, lty=2) # adds a dashed, horizontal line to the plot
# It looks like a linear model is a very bad representation of the data!

# For comparison
x <- 1:100
y <- 2*x + rnorm(100, sd=10)
plot(y ~ x)
xy.lm <- lm(y ~ x)
plot(xy.lm$residuals ~ x); abline(h=0, lty=2)

# Getting help is easy
help(women) # This gives a description of what is in the women dataset
?women # Has the same effect as above
?lm # A good explanation of what the `lm' function does
help.start() # Opens up your web browser with the R help files

# EXERCISE: Go through the excellent examples provided in the `lm'
# documentation.

# EXERCISE: Try `plot(women.lm)'. Can you obtain all the plots simultaneously?
# Check help(par) and look for the `mfrow' parameter. R has many graphical
# parameters that can be set to gain fine control of all the plots. HINT: try
# par(mfrow=c(2,2)) and then plot(women.lm). Try also `plot(xy.lm)'.

example(lm) # This will run the example in the `lm' documentation automatically

# If you don't know the name of the function you're looking for, try
# `help.search'
help.search("gamma distribution")

# R has fast, buit-in pseudo random number generators for any distribution you
# can think of: binomial, multinomial, beta, gamma, poisson, etc, etc, etc.
x <- matrix(rcauchy(10000), ncol=100) # The Cauchy distribution
dim(x) # x has 10,000 elements divided into 100 rows and 100 columns
means <- apply(x, 1, mean)
plot(means, type='l') # The mean of the Cauchy distribution is not defined!

# Let's increased the sample size
x <- matrix(rcauchy(1e6), ncol=1e3)
means <- apply(x, 1, mean)
plot(means, type='l') # The sample mean wiggles erratically!

# EXERCISE: what does the `apply' function does? Find out using the `R'
# documentation. `apply' is much more efficient than using `for' loops
# for traversing along the rows or columns of a matrix. It can dramatically
# increase the speed of code originally written with for loops! Other useful
# functions are `sapply', `lapply', `aggregate' and `replicate'.

# EXERCISE: repeat the example above and convince yourself that for the Normal
# distribution the sample mean quickly converges to zero with sample size, while
# the Cauchy distribution will never converge!

# R has also `array' data structures. Arrays are multidimensional matrices
# with an arbitray number of dimensions!
cube <- array(1:27, dim=c(3,3,3)) # this creates a cube of 3x3x3

cube[3,3,3] # [1] 27

# The contents of `cube' can also be accessed linearly
cube[16] # This is the same as cube[1,3,2]

require(MASS) # This loads the MASS package with its corresponding functions
# we need this for the galaxies and geyser datasets

# Example from p.129 in Venables and Ripley (2002) Modern Applied Statistics with S.
# 4th ed. Springer
# This is the bible for `R' and `S' (S is R's grandaddy!)

gal <- galaxies/1e3 # The speed of galaxies

hist(gal, n=20) # A histogram of galaxy speeds

# We can do something much nicer
plot(x = c(0, 40), y = c(0, 0.18), type="n", bty="l",
xlab="Velocity of galaxy (1000Km/s)", ylab="density")
rug(gal) # adds the data points as a rug at the bottom of the plot
lines(density(gal, width="SJ-dpi", n=256), lty=1)
lines(density(gal, width="SJ", n=256), lty=3)
abline(h=0, lty=2)

# Once you have a nice plot you want to keep
dev.print() # will send the plot to the network printer
dev.print(pdf, file="myplot.pdf") # will save the plot as a pdf file
dev.copy2eps(file="myplot.eps") # will save it as an eps file

# EXERCISE: try `?dev2bitmap' for saving the graphics device as a bitmap file.
# Try `?dev.list' for a series of functions about graphic devices and what they
# do.

# The Old Faithful Geyser (Example from the book above, p.131)
# We'll plot the waiting time between geyser eruptions vs. the previous
# wating time between eruptions
geyser2 <- data.frame(as.data.frame(geyser)[-1,], # removes the first line in geyser
pduration = geyser$duration[-299])
attach(geyser2) # when attaching, the internal data frame names become visible
plot(pduration, waiting)
# without attaching, the command above would have been
# plot(waiting ~ pduration, data=geyser2)

# In a similar fashion to the example above, we can fit a two dimensional kernel
# to the data in order to vizualise the bivariate distribution
geyser.k2d <- kde2d(pduration, waiting, n=50)
image(geyser.k2d)
contour(geyser.k2d, add=TRUE)
detach(geyser2)

# EXERCISE: try contour on its own. HINT: set add=FALSE (its default value!)

# We can also plot a neat perspective of the bivariate distribution
persp(geyser.k2d, phi=30, theta=20, d=5,
xlab="Previous duration", ylab="Waiting time", zlab="")

# EXERCISE (HARD): check the documentation for `persp' first. Get rid of the
# bounding box and the facet edges. Colour the facets using the same coloring
# scheme as in the `image' example above (HINT: check out `heat.colors'), add
# light and shades.

# Sometimes we need to calculate confidence intervals for estimates from arbitrary
# distributions. For example, the galaxies example above shows a weird distribution
# and using normality assumptions to calculate a CI is not appropriate. Bootstrap
# is a neat way to calculate non-parametric CI's and it is very easy to do in R.

set.seed(101) # It is useful to specify the random set as to be able to reproduce
# the sampling.
m <- 1000; gal.boot <- numeric(m)
for (i in 1:m) gal.boot[i] <- median(sample(gal, replace=T))

plot(density(gal.boot))
quantile(gal.boot, p=c(0.025, 0.975)) # 95% CI for the median

# R has very efficient built-in functions for bootstrapping
library(boot)
set.seed(101)
gal.boot2 <- boot(gal, function(x, i) median(x[i]), R = 1000)
gal.boot2
# The CI calculated above are known as the percentile CI's. They are not
# considered the best. `boot.ci' will calcualte a collection of CI's
# according to different methods. Please consult the book above for details
# about them.
boot.ci(gal.boot2)
plot(gal.boot2) # for some diagnostic plots

# Non-parametric regression with lowess:
help(GAGurine) # concetration of a chemical GAG in 314 childern
attach(GAGurine)
plot(GAG ~ Age, pch='+')
lines(lowess(GAG ~ Age, f=1/4), col="red", lwd=2)

# The lowess function fits a non-parametric regression line. It works sort like
# a running mean, but it has a more sophisticated algorith. The `f' parameter
# controls the influence of neighbouring values when fitting the regression at
# a particular point. EXERCISE: repeat the example above with `f=2/3'.

# EXERCISE: checkout `loess' a newer and more sophisticated version of `lowess'.
# It works in a different way.

# For more info about R check www.r-project.org
# There is a very large list of contributed packages to do any sort of
# analyses imaginable. Check http://www.stats.bris.ac.uk/R/ for an extensive
# list of packages. The R mailing list (available form the main R page above)
# provides an excellent forum to get answers to anything about R.

# A few more advanced examples: (added Feb 2015)

# The `curve' function can be used to plot mathematical functions:
curve(dnorm(x, mean=0, sd=1), from=-4, to=4) # plots the normal distribution

# R adds a bit of blank space between what is being plotted, and the plotting box
# you can remove the space using `xaxs' and `yaxs'. y-axis labels can be set to be
# horizontal with `las', we can change the box with `bty' and we can change
# the axis labels and add a title:
curve(dnorm(x, mean=0, sd=1), from=-4, to=4, xaxs="i", yaxs="i", las=1,
bty="l", ylab="Density", main="The normal distribution")
# you can add additional curves to the same plot:
curve(dnorm(x, mean=0, sd=2), from=-4, to=4, add=TRUE, lty=2)

# more functions:
curve(sin(x), from=0, to=2*pi); abline(h=0, lty=2)
# (if you use a semicolon, you can put several R commands in the same line!)

curve(x^2, from=-1, to=1, xaxs="i", yaxs="i")

# R can calculate integrals
# define the x^2 function:
x2 <- function(x) x^2
integrate(x2, lower=0, upper=1) # 0.3333333 with absolute error < 3.7e-15

# integrate a statistical distribution
integrate(dnorm, lower=0, upper=Inf) # 0.5 with absolute error < 4.7e-05

# R can also calculate derivatives:
D(expression(x^2), name="x")
D(expression(y * x^2), name="x") # derivative with respect to x
D(expression(y * x^2), name="y") # with respect to y

# For any questions, comments, corrections or suggestions, please write to:
# mariodosreis@gmail.com
```

Hope you find it useful!
