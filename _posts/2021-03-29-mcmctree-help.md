---
layout: post
title:   "Help, my MCMCtree analysis doesn't work!"
author: "Mario dos Reis"
---

I regularly receive emails from people who are struggling with their
molecular-clock dating analyses in MCMCtree and who do not understand what may
be going wrong. Here, I comment on some common issues you may encounter when
using MCMCtree and how to address them.

# Do you have a fossil calibration on the root?

MCMCtree requires a calibration density on the age of the root to be explicitly
provided, either in the tree file or in the control file. Not having a properly
set up root calibration is a common problem with MCMCtree. For example, consider
the following MCMCtree tree file with one calibration:

```
((A, B)'B(2, 3)', C);
```

Species A and B have a calibration with soft bounds between 2 and 3. If the time
unit is, say, 100 My, then this calibration means the divergence of A and B is
between 200 and 300 Ma. Note this tree does not have a root calibration, thus,
MCMCtree will need one to be specified in the control file. Make sure the
calibration in the control file is sensible. For example, if `'<1'`, which means
the root must be younger than 1, or in this case, 100 Ma, is in the control
file, then the analysis would be unsafe because this root age is younger and in
conflict with the calibration on A and B. An analysis with such calibrations
would produce non-sensical results.

**Recommendations:**
* Make sure you explicitly provide a root calibration.
* Run MCMC analysis with `usedata = 0` in the control file to
obtain the prior on times. Examine the resulting time prior and check the
prior age of the root makes sense and is not in conflict with the prior age of
other nodes.

# Are your calibrations consistent?

MCMCtree does not check for consistency of calibrations. Particularly, if using
L (minimum bound) or SN (Skew-Normal) or ST (Skew-t), it is possible for the
specified calibration densities to be strange. It is also possible that you may
place a calibration on a part of the tree for which that calibration was
unintended.

**Recommendations:**
* Use a graphical tree program (such as FigTree) to plot the tree file.
Calibrations are labels on nodes, thus, visualising node labels with your tree
program will let you inspect the calibrations and make sure they are placed
correctly. Check for silly errors such ``'B(2,1)'`` (the minimum bound is older
than the maximum) or for younger nodes having older calibrations than older
nodes.
* Use the `sn` package in R to plot skew-normal and skew-t calibrations. Make
sure the densities are placed in the intended time range and that they are not
in conflict with calibrations on other nodes.
* Use our `mcmc3r` R package[^1] to plot the B and L calibrations.
* Run MCMC analysis with `usedata = 0` in the control file to obtain the prior on
times. Examine the resulting time prior and check the prior ages on nodes makes
sense. Plot the prior densities on node ages in R and overlay them with the
analytical calibration densities plotted with the packages above.

# My analysis does not converge. What do I do?

When analysing very large datasets convergence can be quite slow. In an analysis
of 372 primate species, it took about 15 days and 8 parallel MCMC runs to achieve
convergence. With more species and longer alignments convergence may need even
longer runs.

**Recommendations:**
* Build a small dataset first. That is, start with 10 to 20 taxa and use a
relatively small alignment to calculate the gradient and Hessian (e.g. the in.BV
file). Then run your analysis and make sure you get sensible results.
* You can then increase data size progressively. Note calculation of the
approximate likelihood increases with the square of the number of species, and
calculation of the time and rate priors increase linearly with number of species.
This will give you a rough idea of how long your full data set should take to run.

# How do I choose the best clock model?

MCMCtree now implements sampling from power posteriors, which allows you to
perform Bayesian selection of clock model, tree topology and substitution model.
Unfortunately, the power posteriors must be calculated using exact likelihood and
the approximate method can not be used. See my [blog post](https://dosreislab.github.io/2017/10/24/marginal-likelihood-mcmc3r.html)
on Bayesian model selection with MCMCtree.

# Are there any tutorials to learn how to use MCMCtree?

I have written one book chapter and a set of tutorials:

* dos Reis and Yang (2019) **Bayesian molecular clock dating using genome-scale
datasets**. In: Anisimova (ed.) Evolutionary Genomics. Methods in Molecular
Biology, vol 1910. Humana, New York, NY. [DOI:
10.1007/978-1-4939-9074-0_10](https://doi.org/10.1007/978-1-4939-9074-0_10)
* dos Reis, Álvarez-Carretero and Yang (2017) **MCMCtree tutorials**. Available
online [here](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf)

# How do I cite MCMCtree?

A general citation for MCMCtree is

* Rannala and Yang. (2007) Inferring speciation times under an episodic molecular
clock. Systematic Biology, 56:453-466.

If you use the approximate likelihood method please cite

* dos Reis and Yang (2011) Approximate likelihood calculation for Bayesian
estimation of divergence times. Molecular Biology and Evolution,
28:2161-2172.

If you use the `mcmc3r` package for Bayesian selection of clock model please
cite

* dos Reis et al. (2018) Using phylogenomic data to explore the effects of
relaxed clocks and calibration strategies on divergence time estimation: Primates
as a test case. Systematic Biology, 67: 594–615.

-----
[^1]: [https://github.com/dosreislab/mcmc3r](https://github.com/dosreislab/mcmc3r)
