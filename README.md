# HDCovMatTest

![GitHub R package version](https://img.shields.io/github/r-package/v/shawnsanto/HDCovMatTest)

Homogeneity tests of covariance matrices with high-dimensional longitudinal data

## Overview

This repository hosts an R package that implements the methodology developed in
"Homogeneity tests of covariance matrices with high-dimensional longitudinal 
data" by Ping-Shou Zhong,  Runze Li, and Shawn Santo. 

That paper deals with the detection and identification of changepoints among 
covariances of high-dimensional longitudinal data, where the number of features 
is greater than both the sample size and the number of repeated measurements. 
The proposed methods are applicable under general temporal-spatial dependence. 
A new test statistic is introduced for changepoint detection, and its asymptotic 
distribution is established. If a changepoint is detected, an estimate of the 
location is provided. The rate of convergence of the estimator is shown to 
depend on the data dimension, sample size, and signal-to-noise ratio. Binary 
segmentation is used to estimate the locations of possibly multiple 
changepoints, and the corresponding estimator is shown to be consistent under 
mild conditions. Simulation studies provide the empirical size and power of the 
proposed test and the accuracy of the changepoint estimator. An application to 
a time-course microarray dataset identifies gene sets with significant gene 
interaction changes over time.

## Install package `HDCovMatTest`

Package HDCovMatTest can be installed using the `devtools` package in R. 

```r
# install (if not available) and load the devtools package
if(!require(devtools)) {
    install.packages("devtools")
    library(devtools)
}

# install HDCovMatTest
install_github("shawnsanto/HDCovMatTest")

# load HDCovMatTest
library(HDCovMatTest)
```

## Example

Consider a small simulated data example with a change point at time 2.

```r
# set primary parameters
p <- 30 # data dimension
n <- 10 # number of individuals
TT <- 5 # total number of time points

# set secondary parameters for simulation and data generation
delta <- 0.35
k0 <- 2

m <- p + 20; L <- 3; w <- 0.2
```

Generate the data. Here `delta` is chosen to be large so we ensure a change 
point exists at time 2.

```r
set.seed(928)

gamma_1 <- gamma_2 <- matrix(0, p, m * L)
y <- array(0, c(p, n, TT))

for (i in 1:p) {
  for (j in 1:p) {
    dij <- abs(i - j)

    if (dij < (p * w)) {
      gamma_1[i, j] <- (dij + 1) ^ (-2)
      gamma_2[i, j] <- (dij + 1 + delta) ^ (-2)
    }
  }
}

Z <- matrix(rnorm(m * (TT + L - 1) * n), m * (TT + L - 1), n)

for (t in 1:k0) {
  y[, , t] <- gamma_1 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
}

for (t in (k0 + 1):TT) {
  y[, , t] <- gamma_2 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
}
```

Perform the testing procedure.

```r
test_covmat(y, n, p, TT, alpha = 0.01)
```

Results are provided in an R list.

```r
$reject # rejection flag
[1] 1

$estcp # estimated change point
[1] 2

$tstat # test statistic
[1] 3.501253

$pvalue # p-vale
[1] 0.0005785431

$corrmat # corresponding correlation matrix
          [,1]      [,2]      [,3]      [,4]
[1,] 1.0000000 0.8423203 0.8191572 0.7686094
[2,] 0.8423203 1.0000000 0.9790518 0.9195994
[3,] 0.8191572 0.9790518 1.0000000 0.9471566
[4,] 0.7686094 0.9195994 0.9471566 1.0000000
```

## Reference

Ping-Shou Zhong, Runze Li, Shawn Santo, Homogeneity tests of covariance matrices 
with high-dimensional longitudinal data, Biometrika, Volume 106, Issue 3, 
September 2019, Pages 619–634, https://doi.org/10.1093/biomet/asz011.
