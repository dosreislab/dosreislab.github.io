---
layout: post
title:  "bppr: Calibrating a BPP phylogeny to geological time"
author: "Mario dos Reis"
---

BPP estimates relative node ages (called tau's) in substitutions per site assuming the molecular clock. If we have information on the age of a node (for example, from the fossil record), or information about the molecular rate (for example, measurements of per-generation rate from parent-offspring sequencing), then it is possible to use a random sampling procedure to calibrate the phylogeny (see Angelis and dos Reis, 2015 and Yoder et al. 2016). Note that BPP also provides samples of the nucleotide diversities (aka theta's =4Nu). These can be converted into population sizes (N) by dividing by the per-generation rate (u), if the latter is available from a prior.

The package comes with two datasets, `hominids` and `microcebus`. Both contain MCMC samples (`mcmc.txt`) from BPP A00 analyses of hominid and _Microcebus_ spp. phylogenies. A BPP A00 analysis uses a fixed species tree to estimate the demographic parameters (the tau's and theta's, See Yang, 2015, for details). We will use the hominid samples to calibrate this phylogeny to geological time. You can check the posterior means for the tau's and thetas for the hominid dataset in R:

```R
# For the hominid dataset:
# Calculate posterior means
apply(bppr::hominids$mcmc, 2, mean)

# Calculate 95% CIs:
t(apply(bppr::hominids$mcmc, 2, quantile, probs=c(.025,.975)))

# If you have the coda package, you can get 95% HPD CIs:
coda::HPDinterval(coda::as.mcmc(bppr::hominids$mcmc))

# If you have ape installed you can plot the phylogeny:
ape::plot.phylo(bppr::hominids$tree)

```
**Exercise:** Repeat the calculations for the `microcebus` dataset.

For example, the estimated branch length from human to the human-chimp ancestor (`tau_7humanchimp`) is 4.06e-3 (95% CI: 3.96e-03  4.16e-03; 95% HPD: 3.96e-03  4.16e-03), and the nucleotide diversity for the human-chimp ancestral population (`theta_7humanchimp`) to be 6.09e-03 (95% CI: 5.66e-03  6.53e-03; 95% HPD: 5.65e-03  6.52e-03).

### 1. Calibrating using a prior on the rate

Sequencing studies place estimates of the per-generation mutation rate in human at around 1.2e-8 to 1.4e-8 (see Scally and Durbin, 2012, for a review). Estimates of the generation time in apes range around 20 to 30 years (Langergraber et al. 2012). We can use these values to construct prior densities on the per-generation mutation rate, u, and the generation time, g. From these we can then estimate the per-year rate (r = u / g), and use this estimate to convert the relative node ages (tau's) into geological divergence times. For the per-generation rate we will assume a gamma prior with mean 1.3e-8 and standard deviation 0.1e-8, which gives a 95% prior CI of roughly 1.1e-8 to 1.5e-8, to accommodate uncertainty in estimates of the rate. For the generation time, we will use a gamma density with mean 25, and standard deviation 2.5, which gives a 95% prior CI of roughly 20 to 30 years. Then function `bppr::msc2time.r` will carry out the random sampling procedure of Angelis and dos Reis (2015) to generate the times. In R:

```R
ape.time <- bppr::msc2time.r(bppr::hominids$mcmc, u.m = 1.3e-8, u.sd = .1e-8, g.m = 25, g.sd = 2.5)

# posterior means:
apply(ape.time, 2, mean)

# posterior HPDs:
coda::HPDinterval(coda::as.mcmc(ape.time))
```

For example, the estimate for the human-chimp divergence time (`t_7humanchip`) is 7.85e+06 (95% HPD: 5.95e+06 9.80e+06) which corresponds to 7.85 Ma (5.95 - 9.80 Ma). These estimates appear reasonable. Benton et al. (2015) give an estimate based on careful consideration of the fossil record of 6.5 to 10 Ma. Because we have a prior on the per-generation rate, the function will convert the theta's (=4Nu) into effective population sizes. For example, for the human-chimp ancestral population (`Ne_7humanchimp`) we get an estimate of 118 thousand individuals (95% HPD: 98.2, 137 K individuals). Note that the `hominids` dataset is based on an analysis of almost 15,000 loci (Angelis and dos Reis 2015, Burgess and Yang, 2008), and thus the estimates of tau's and theta's are reasonably precise.

## 2. Calibrating using a prior on a node age

We can use Benton et al. (2015) calibration of 6.5 to 10 Ma for the human-chimp divergence to recalibrate the ape phylogeny. This is done with the `bppr::msc2time.t` function, which uses a prior on a node age to recalibrate the phylogeny. Note that `bppr` can use a calibration on any node on the phylogeny, but only one calibration can be used at a time.

First, we calibrate using an uniform distribution between 6.5 and 10 Ma. In R:

```R
ape.time2 <- bppr::msc2time.t(bppr::hominids$mcmc, node.name="7humanchimp", calf=runif, min=6.5, max=10)

# posterior means:
apply(ape.time2, 2, mean)

# posterior HPDs:
coda::HPDinterval(coda::as.mcmc(ape.time2))
```

Our estimates here of the human-chimp divergence time are simply 6.5-10 as this is the prior age. However, we obtain estimates for the other node ages. For example, for the root of the phylogeny (node 5, crown apes), the estimated age is 28.0 Ma (HPD: 22.3, 33.6 Ma). In the previous rate-calibrated analysis the estimate was 26.6 Ma (HPD: 20.1, 33.1 Ma).

We can use other random distributions to construct the calibration (virtually any function that generates random deviates in R can be used). For example, Benton et al. (2015) suggest the maximum of 10 Ma be soft, that is, they acknowledge there is a substantial probability that the age of the human-chimp ancestor could be older. We can model this using a shifted log-normal distribution. The _shift_ is a number that moves the distribution to the left or right. Here we will use a shift of 6.5 Ma, corresponding to the minimum bound, and then we will adjust the parameters of the distribution so that the upper 90% limit of the distribution is roughly 10 Ma. In R:

```R
# A description of the shifted log-normal is given in the helpfile:
?bppr::ShiftedLognormal

# By trial and error, I chose values for the shifted log-normal
# that gave a 10% limit close to 10 Ma:
curve(bppr::dslnorm(x, shift=6.5, meanlog=0, sdlog=1), from=0, to=15, n=1e3)

# The 2.5% and 90% limits:
bppr::qslnorm(c(.025,.9), shift=6.5, meanlog=0, sdlog=1)
# [1]  6.640863 10.102224

ape.time3 <- bppr::msc2time.t(bppr::hominids$mcmc, node.name="7humanchimp", calf=bppr::rslnorm, shift=6.5, meanlog=0, sdlog=1)

# posterior means:
apply(ape.time3, 2, mean)

# posterior HPDs:
coda::HPDinterval(coda::as.mcmc(ape.time3))
```

The age of the root in this case is 27.6 Ma (95% HPD: 21.8, 39.7 Ma). Note the upper bound of the root age is now older. This is a consequence of the longer tail in the shifted log-normal calibration.

**Excercise:** Look into the `rate` estimates (this is the molecular rate per geological time) for the various time-calibrated analysis. Are they similar? Do the CIs or HPDs for the two analysis overlap?

**Excercise:** Yoder et al. (2016) describe priors for the generation time and per-generation rate for the mouse lemurs (_Microcebus_ spp.). Use these priors to calibrate dataset `microcebus` to geological time.

## Bibliography

* Angelis, K. and dos Reis, M. 2015. _The impact of ancestral population size and incomplete lineage sorting on Bayesian estimation of species divergence times_. Current Zoology, 61: 874–885.

* Benton et al. 2015. _Constraints on the timescale of animal evolutionary history_ Palaeontologica Electronica, 18.1.1FC.

* Burgess, R. and Z. Yang. 2008 _Estimation of hominoid ancestral population sizes under Bayesian coalescent models incorporating mutation rate variation and sequencing errors_. Molecular Biology and Evolution, 25: 1979-1994.

* Flouri T, Xiyun J, Rannala B, Yang Z. 2018. _Species tree inference with BPP using genomic sequences and the multispecies coalescent_. Mol Biol Evol.

* Langergraber et al. 2012. _Generation times in wild chimpanzees and gorillas suggest earlier divergence times in great ape and human evolution_. Proceedings of the National Academy of Sciences, 109: 15716-15721.

* Scally and Durbin. 2012. _Revising the human mutation rate: implications for understanding human evolution_. Nature Reviews Genetics volume 13, pages 745–753.

* Yang (2014) _Molecular Evolution: A Statistical Approach_. Oxford University Press.

* Yang, Z. 2015. _The BPP program for species tree estimation and species delimitation_. Current Zoology, 61: 854-865.

* Yoder, A. et al. 2016. _Geogenetic patterns in mouse lemurs (genus Microcebus) reveal the ghosts of Madagascar’s forests past_. Proceedings of the National Academy of Sciences, 113: 8049–8056.
