# HDCovMatTest
### Homogeneity tests of covariance matrices with high-dimensional longitudinal data

#### Description

This repository hosts an R package that implements the methodology developed in
"Homogeneity tests of covariance matrices with high-dimensional longitudinal data". 
That paper considers detecting and identifying change points among covariances of high-dimensional
longitudinal data, where the number of features is greater than both the sample size
and the number of repeated measurements. The proposed methods are applicable under general
temporospatial dependence.

#### Install HDCovMatTest package

Package HDCovMatTest can be installed using the devtools package in R. Below is R code to obtain HDCovMatTest.

```R
# install (if not available) and load the devtools package
if(!require(devtools)){
    install.packages("devtools")
    library(devtools)
}

# install HDCovMatTest
install_github("shawnsanto/HDCovMatTest")

# load HDCovMatTest
library(HDCovMatTest)
```


#### Reference

Zhong, Li, and Santo (2018). Homogeneity tests of covariance matrices with high-dimensional longitudinal data. *Biometrika*.
