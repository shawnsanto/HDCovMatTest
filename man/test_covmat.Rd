% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_covmat.R
\name{test_covmat}
\alias{test_covmat}
\title{Test of high-dimensional covariance matrices with longitudinal data}
\usage{
test_covmat(y, n, p, TT, alpha = 0.01)
}
\arguments{
\item{y}{A high-dimensional longitudinal data set in the format of a three
dimensional array where the first coordinate is for features, the second
coordinate is for sample subjects, and the third coordinate is for time
repetitions. Thus, the dimension of y is \eqn{p x n x TT} where
\eqn{p} is the data dimension, \eqn{n} is the sample size, and \eqn{TT}
is the number of repetition times.}

\item{n}{The number of individuals.}

\item{p}{The data dimension.}

\item{TT}{The number of repetition times. It is recommended that
\eqn{TT \le 30} so the asymptotics in Zhong, Li, and Santo (2019)
can be applied.}

\item{alpha}{The type I error of the homogeniety test. The nominal level for
the quantile is computed as \eqn{1-alpha}. Suggested values for alpha
include 0.01 (default) and 0.05.}
}
\value{
The function returns a test result, estimated change point,
  test statistic, p-value, and correlation matrix. The output is provided in
  a list.
  \describe{
    \item{$reject}{Null hypothesis rejection indicator. A value of 1
      indicates the null hypothesis is rejected. The null hypothesis is that
      all the covariance matrices are equal across time.}
    \item{$estcp}{The first estimated change point provided the null
      hypothesis is rejected. This value will be 0 if the null hypothesis is
      not rejected.}
    \item{$teststat}{The test statistic.}
    \item{$pvalue}{The p-value.}
    \item{$corrmat}{The test statistic is a maximum of \eqn{TT-1}
      standardized statistics which quantifies the Frobenius norm of
      covariance matrices before and after time \eqn{t} for
      \eqn{t=1,...,TT-1}. The correlation matrix is the correlation matrix
      among the \eqn{TT-1} standardized statistics. For further details,
      see Zhong, Li, and Santo (2019).}
    }
}
\description{
This function implements a test procedure proposed by
Zhong, Li, and Santo (2019) for testing the homogeneity of covariance
matrices in high-dimensional longitudinal data. Temporal and spatial
dependence are allowed. The null hypothesis of the test is
that the covariance matrices at all repetition times are equal.
The methodology has been proved, tested, and shown to work well for
high-dimensional longtidudinal data where both \eqn{p} and \eqn{n} diverge,
but no explicit relationship between \eqn{p} and \eqn{n} is needed.
Furthermore, the current methodology works for longitudinal data with a small
number of repetition times, i.e., \eqn{TT < \infty}.
}
\examples{
# A testing example with a change point at time 2

# Set parameters
p <- 30; n <- 10; TT <- 5
delta <- 0.35
m <- p + 20; L <- 3; k0 <- 2; w <- 0.2

# Generate data
gamma_1 <- gamma_2 <- matrix(0, p, m * L)
y <- array(0, c(p, n, TT))
set.seed(928)

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
  y[, , t] <- gamma_1 \%*\% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
}

for (t in (k0 + 1):TT) {
  y[, , t] <- gamma_2 \%*\% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
}

test_covmat(y, n, p, TT, alpha = 0.01)
}
\references{
\emph{Ping-Shou Zhong, Runze Li, Shawn Santo, Homogeneity tests of
covariance matrices with high-dimensional longitudinal data, Biometrika,
Volume 106, Issue 3, September 2019, Pages 619–634,
https://doi.org/10.1093/biomet/asz011.}
}
\author{
\strong{Maintainer}: Shawn Santo \email{shawn.santo@duke.edu}

  Authors:
    \itemize{
      \item Ping-Shou Zhong
      \item Runze Li
      \item Shawn Santo
    }
}
