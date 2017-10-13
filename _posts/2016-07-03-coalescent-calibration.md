---
layout: post
title:  "How to calibrate coalescent trees to geological time"
author: "Mario dos Reis"
---

There has been much interest recently in using the multi-species coalescent to estimate species phylogenies from molecular data. The advantage of the method is that incomplete lineage sorting (the discrepancy between gene trees and the species tree), and ancestral polymorphism are accounted for during phylogenetic inference. Bayesian implementations of the method are computationally expensive, and are best suited for inference among closely related species (or populations).

The figure below shows an example of a phylogeny of mouse lemurs (_Microcebus_ spp.) estimated from RADseq (restriction site associated DNA sequencing) data using the program BPP (Yang, 2015), which implements the multi-species coalescent. Each node in the tree represents a speciation event, with the node ages given as numbers of substitutions per site. The blue bars represent the 95% credibility interval of the node age. For example, the molecular distance from the last common ancestor of _M. rufus_ and _M. berthae_ to the present is 1.29 × 10<sup>–4</sup> substitutions per site. If we knew the molecular substitution rate per year for mouse lemurs, we could calibrate the tree to geological time, that is, we could convert the node ages from units of substitution per site to units of real time. I'll explain how to do so in this post.

![](/assets/figs/microcebus-tree.png)  

# A simple method

The mouse lemur genome is roughly 2.8 Gbp long, with perhaps the majority of the genome evolving neutrally. The RADseq data used to build the phylogeny above represents a random sample of the genome, and thus we may assume that most RAD-fragments are evolving neutrally. Estimates of the de novo mutation rate in the human genome are around 1.2 × 10<sup>–8</sup> substitutions per site per generation (see Scally and Durbin, 2012, for a review). Given that mouse lemurs are also primates, we could use the human mutation rate as a proxy for the mouse lemur rate. Assuming a generation time of 4 years for mouse lemurs would roughly give us 1.2 × 10<sup>–8</sup> / 4 = 0.3 × 10<sup>–8</sup> substitutions per site per year for the lemur rate. Let _τ_ be the age of a node in the phylogeny in units of substitutions per site. Then it follows that _τ_ = _rt_, where _r_ is the mutation rate per year, and _t_ the age in years. For the _M. rufus_ and _M. berthae_ divergence we have that _τ_ = 1.29 × 10<sup>–4</sup>, and thus _t_ = _τ_/_r_ = 1.29 × 10<sup>–4</sup> / 0.3 × 10<sup>–8</sup> = 43,000 years. We could repeat the procedure to the other nodes in the tree to obtain all the ages of divergence. This is perhaps the simplest procedure to calibrate a coalescent tree to geological time.

The procedure outlined above has at least two serious limitations however. The first is that we used point estimates for _τ_ and the mutation rate, _r_, and we have thus ignored the uncertainties associated with these estimates. The second is that the human mutation rate is probably not a good proxy for the mouse lemur rate. We now look at a more sophisticated way to calibrate the tree that overcomes these issues.

# A full Bayesian approach

The following discussion assumes you have a basic grasp of Bayesian statistics, are familiar with BPP (or other similar Bayesian coalescent inference software), and with the R environment for statistical computing.

The program BPP uses an MCMC Bayesian approach to estimate the model parameters (i.e. the species phylogeny, the ancestral nucleotide diversities, _θ_ = 4_Nu_, and the _τ_'s) from a molecular sequence alignment. Once an MCMC BPP analysis has finished, one obtains a sample of the posterior distribution of the parameters of interest. The sample size is specified by the user when running the software. For example, for the tree above I analysed 80,662 RAD-fragments and took 20,000 samples from the MCMC (see Yoder et al., 2016, for details). In fact, because the analysis is so computationally intensive, I actually ran 8 parallel MCMC analysis, each with 20,000 samples. Here I discuss the results from only one of those analysis, but the results presented in the main paper by Yoder et al. were calculated using all 8 × 20,000 = 160,000 samples.

To obtain a full Bayesian estimate of the divergence times of the mouse lemurs we must first construct priors on the per-generation mutation rate, _u_; and the generation time in years, _g_. One can then obtain samples from these priors, which can then be combined with the samples of _τ_'s from BPP to obtain the desired divergence times (Angelis and dos Reis, 2015). Let _τ<sub>i</sub>_ be the _i_-th sample of a node age from the MCMC, and _u<sub>i</sub>_, and _g<sub>i</sub>_ the _i_-th sample obtained from the priors. Then, the _i_-th sample of the mutation rate per year is given by _r<sub>i</sub>_ = _u<sub>i</sub>_ / _g<sub>i</sub>_. Finally, the _i_-th sample of the divergence time in years is _t<sub>i</sub>_ = _τ<sub>i</sub>_ / _r<sub>i</sub>_. Thus, one simply repeats the procedure for as many times as the MCMC sample size to obtain the full posterior distribution of _t_. The whole procedure is then repeated for the other nodes in the phylogeny. Below I exemplify the procedure using the mouse lemur phylogeny.

**A prior on the mouse lemur generation time:** Considerations of the reproductive biology and lifespan in wild and captive populations of mouse lemurs give an estimate of 3 to 4.5 years for the generation time in this species (Yoder et al. 2016). Thus we can construct a gamma prior for _g_ with mean 3.75 years (the mid-point of 3–4.5) and standard deviation 0.375 years (so that 3.75 ± 2 × 0.375 equates to 3–4.5 years). Recall that the gamma distribution with parameters _α_ and _β_ has mean _α_/_β_ and variance _α_/_β_<sup>2</sup>. Thus, the prior for _g_ is

_f_ (_g_) = Gamma (_g_ \| _α_ = 100, _β_ = 26.67).  (1)

**A prior on the mouse lemur per-generation mutation rate:** Estimates of the mutation rate in the lab mouse are about 0.54 × 10<sup>–8</sup> substitutions per site per generation (Uchimura et al. 2015), while, as mentioned above, they are about 1.2 × 10<sup>–8</sup> substitutions per site per generation in human. It is unclear whether the per-generation mutation rate in the mouse lemurs should be similar to the lab mouse (which has a similar body size and lifespan), or to the human (a fellow Primate), thus, we construct a conservative gamma prior on _u_ between 0.5–1.2 (× 10<sup>–8</sup>) substitutions per site per generation. Using the same procedure as above, we get a mean of 0.87 (× 10<sup>–8</sup>) substitutions per site per generation, and a standard deviation of 0.165 (× 10<sup>–8</sup>) substitutions per site per generation. The prior on _u_ is then

_f_ (_u_) = Gamma (_u_ \| _α_ = 27.80, _β_ = 31.96).  (2)

**Converting the _τ_'s into geological times of divergence:** The table below summarises the procedure. The first column (_i_) indicates the sample number. The second column (_τ<sub>rb</sub>_) is the MCMC posterior sample obtained with BPP of the age of the  _M. rufus_ and _M. berthae_ divergence. The third column (_u_), is the sample of the per-generation mutation rate from the gamma prior of Eq. (1). I generated the sample of _u_ using R, that is, by running the R command `rgamma(2e4, 27.80, 31.96)`. The fourth column (_g_), is the sample of the generation time from the gamma prior of Eq. (2), also obtained with R. The fifth column (_r_) is the prior sample of the per-year mutation rate, obtained by combining the priors on _u_ and _g_. Finally, the last column (_t<sub>rb</sub>_) is the resulting sample of the posterior distribution of the geological time of divergence of _M. rufus_ and _M. berthae_, obtained by combining the prior on _r_ with BPP's posterior of _τ<sub>rb</sub>_. One can summarise the posterior sample of _t<sub>rb</sub>_ in the usual way, to obtain the posterior mean, variance, and 95% credibility interval (CI). The posterior mean is 57,545 years before present with 95% CI: 142,193–10,207 years. The table below can be uploaded into a program such as Tracer, to analyse the traces of _r_ and _t<sub>rb</sub>_, calculate the ESS, etc.

|_i_|_τ<sub>rb</sub>_ (× 10<sup>–4</sup>)|_u_ (× 10<sup>–8</sup>)|_g_|_r_ (= _u_/_g_ × 10<sup>–8</sup>)|_t<sub>rb</sub>_ (= _τ_/_r_)|
|:---:|:---:|:---:|:---:|:---:|:---:|
|1|1.632|0.742|4.109|0.181|90,372|
|2|1.531|0.903|4.033|0.224|68,378|
|3|1.566|1.124|2.956|0.380|41,165|
|. . .|. . .|. . .|. . .|. . .|. . .|
|20,000|1.599|0.758|3.673|0.206|77,473|

<br>
The advantage of the procedure just outlined is that one obtains the full, correct posterior distribution of the divergence times integrating all the uncertainties involved. Note that this procedure allows the use of arbitrary prior distributions for _u_ and _g_ (as long as one can obtain samples for them) without having to make any complex analytical calculations. The procedure can be easily extended to calculate the ancestral effective population sizes, _N_, from BPP's posterior sample of the nucleotide diversities (_θ_ = 4_Nu_). The reader should be able to work out how to do this easily. The method can also be trivially extended to deal with the output of other programs and methods (such as, for example, the pairwise sequentially Markovian coalescent or PSMC).

You can read more about our mouse lemur research, and see how estimating their times of divergence let us get a glimpse of Madagascar's forests past, at [The Washington Post](https://www.washingtonpost.com/news/speaking-of-science/wp/2016/07/18/this-absurdly-adorable-mouse-lemur-lets-scientists-travel-back-in-time/) and [National Geographic](http://voices.nationalgeographic.com/2016/07/19/ridiculously-cute-mouse-lemurs-hold-key-to-madagascars-past/).

The files and code needed to reproduce this tutorial are available from [github.com/mariodosreis/mousies](http://github.com/mariodosreis/mousies).

# References

Angelis K, and dos Reis M. (2015) **[The impact of ancestral population size and incomplete lineage sorting on Bayesian estimation of species divergence times.](http://cz.oxfordjournals.org/content/61/5/874)** Current Zoology, 61: 874–885.

Scally A, and Durbin R. (2012) **[Revising the human mutation rate: implications for understanding human evolution.](http://www.nature.com/nrg/journal/v13/n10/full/nrg3295.html)** Nature Reviews Genetics, 13: 745–753.

Uchimura A, et al. (2015) **[Germline mutation rates and the long-term phenotypic effects of mutation accumulation in wild-type laboratory mice and mutator mice.](http://genome.cshlp.org/content/early/2015/06/30/gr.186148.114.abstract)** Genome Research, 25: 1125–1134.

Yang Z. (2015)  **[The BPP program for species tree estimation and species delimitation.](http://cz.oxfordjournals.org/content/61/5/854)** Current Zoology, 61: 854–865.

Yoder AD, et al. (2016) **[Geogenetic patterns in mouse lemurs (genus _Microcebus_) reveal the ghosts of Madagascar’s forests past.](http://www.pnas.org/content/early/2016/07/13/1601081113.full)** Proceedings of the National Academy of Sciences, 113: 8049–8056.
