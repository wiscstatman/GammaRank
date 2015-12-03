# GammaRank

This R package calculates the probability that a set of independent gamma (or normal) variates attains a certain order.
For the Gamma case, it takes as input a vector of shape and scale parameters.  Of course, all orders are equally probable (1/n!)
if the parameters are shared among the variables, but the rank probability is much more difficult to compute in the general case.

The package deploys a special dynamic programming scheme described in a manuscript in preparation.

#Installation

`library(devtools)`

`install_github("wiscstatman/GammaRank")`
