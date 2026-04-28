# Instrumental Variables Quantile Regression (ivqr)

A fork of the `ivqr` package originally by Omkar A. Katta to add the following features:
- Support for HiGHS and ECOS, open source solvers, to avoid a srtict dependency on Gurobi
- Custom weights in quantile regression (mostly useful in cases such as a weighted bootstrap)

Note that the weak instrument tests require a solver that is capable of handling mixed integer conic quadratic resrtictions (micqr), which are currently only ECOS and Gurobi. 
Of the two, only Gurobi has a feasible solve time.

## Installation

```r
# install.packages("remotes") # install `devtools` if not already installed
remotes::install_github("zeyuz35/ivqr")
```
