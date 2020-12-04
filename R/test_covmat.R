
#' Test of high-dimensional covariance matrices with longitudinal data
#'
#' This function implements a test procedure proposed by
#' Zhong, Li, and Santo (2019) for testing the homogeneity of covariance
#' matrices in high-dimensional longitudinal data. Temporal and spatial
#' dependence are allowed. The null hypothesis of the test is
#' that the covariance matrices at all repetition times are equal.
#' The methodology has been proved, tested, and shown to work well for
#' high-dimensional longtidudinal data where both \eqn{p} and \eqn{n} diverge,
#' but no explicit relationship between \eqn{p} and \eqn{n} is needed.
#' Furthermore, the current methodology works for longitudinal data with a small
#' number of repetition times, i.e., \eqn{TT < \infty}.
#'
#' @param y A high-dimensional longitudinal data set in the format of a three
#'   dimensional array where the first coordinate is for features, the second
#'   coordinate is for sample subjects, and the third coordinate is for time
#'   repetitions. Thus, the dimension of y is \eqn{p x n x TT} where
#'   \eqn{p} is the data dimension, \eqn{n} is the sample size, and \eqn{TT}
#'   is the number of repetition times.
#' @param n The number of individuals.
#' @param p The data dimension.
#' @param TT The number of repetition times. It is recommended that
#'   \eqn{TT \le 30} so the asymptotics in Zhong, Li, and Santo (2019)
#'   can be applied.
#' @param alpha The type I error of the homogeniety test. The nominal level for
#'   the quantile is computed as \eqn{1-alpha}. Suggested values for alpha
#'   include 0.01 (default) and 0.05.
#'
#' @return The function returns a test result, estimated change point,
#'   test statistic, p-value, and correlation matrix. The output is provided in
#'   a list.
#'   \describe{
#'     \item{$reject}{Null hypothesis rejection indicator. A value of 1
#'       indicates the null hypothesis is rejected. The null hypothesis is that
#'       all the covariance matrices are equal across time.}
#'     \item{$estcp}{The first estimated change point provided the null
#'       hypothesis is rejected. This value will be 0 if the null hypothesis is
#'       not rejected.}
#'     \item{$teststat}{The test statistic.}
#'     \item{$pvalue}{The p-value.}
#'     \item{$corrmat}{The test statistic is a maximum of \eqn{TT-1}
#'       standardized statistics which quantifies the Frobenius norm of
#'       covariance matrices before and after time \eqn{t} for
#'       \eqn{t=1,...,TT-1}. The correlation matrix is the correlation matrix
#'       among the \eqn{TT-1} standardized statistics. For further details,
#'       see Zhong, Li, and Santo (2019).}
#'     }
#'
#' @author \strong{Maintainer}: Shawn Santo \email{shawn.santo@@duke.edu}
#'
#'   Authors:
#'     \itemize{
#'       \item Ping-Shou Zhong
#'       \item Runze Li
#'       \item Shawn Santo
#'     }
#'
#' @references \emph{Ping-Shou Zhong, Runze Li, Shawn Santo, Homogeneity tests of
#' covariance matrices with high-dimensional longitudinal data, Biometrika,
#' Volume 106, Issue 3, September 2019, Pages 619–634,
#' https://doi.org/10.1093/biomet/asz011.}
#'
#' @export
#'
#' @examples
#' # A testing example with a change point at time 2
#'
#' # Set parameters
#' p <- 30; n <- 10; TT <- 5
#' delta <- 0.35
#' m <- p + 20; L <- 3; k0 <- 2; w <- 0.2
#'
#' # Generate data
#' gamma_1 <- gamma_2 <- matrix(0, p, m * L)
#' y <- array(0, c(p, n, TT))
#' set.seed(928)
#'
#' for (i in 1:p) {
#'   for (j in 1:p) {
#'     dij <- abs(i - j)
#'
#'     if (dij < (p * w)) {
#'       gamma_1[i, j] <- (dij + 1) ^ (-2)
#'       gamma_2[i, j] <- (dij + 1 + delta) ^ (-2)
#'     }
#'   }
#' }
#'
#' Z <- matrix(rnorm(m * (TT + L - 1) * n), m * (TT + L - 1), n)
#'
#' for (t in 1:k0) {
#'   y[, , t] <- gamma_1 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
#' }
#'
#' for (t in (k0 + 1):TT) {
#'   y[, , t] <- gamma_2 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
#' }
#'
#' test_covmat(y, n, p, TT, alpha = 0.01)

test_covmat <- function(y, n, p, TT, alpha = 0.01) {
  yvec <- matrix(y, 1, (p * n * TT))
	rejind <- 0
	khat <- 0

	if (TT > 1) {
  	 khat <- numeric(1)
  	 std_dk <- numeric(TT - 1)
  	 max_std_dk <- numeric(1)
  	 corr_mat <- numeric((TT - 1) * TT / 2)

  	 storage.mode(yvec) <- 'double';

  	 z <- .C('compute_test_cp_r',
  		        as.double(yvec),
          		as.integer(n),
          		as.integer(p),
          		as.integer(TT),
          		khat       = as.integer(khat),
          		std_dk     = as.double(std_dk),
          		max_std_dk = as.double(max_std_dk),
            	corr_mat   = as.double(corr_mat))

  	if (TT > 2) {
  	  offdiag <- z$corr_mat[-c((2 * (TT - 1) - (0:(TT - 2)) + 1) *
  	                            (0:(TT - 2)) / 2 + 1)]
  	  varvec <- matrix(1, TT - 1, TT - 1)
  	  varvec[lower.tri(varvec)] <- offdiag
  	  varvec[upper.tri(varvec)] <- t(varvec)[upper.tri(t(varvec))]
  	  corr_mat_full <- varvec

  	  quant <- mvtnorm::qmvnorm(1 - alpha, corr = varvec)$quantile

  	  pvalue <- 1 - mvtnorm::pmvnorm(lower = -Inf,
  	                                 upper = rep(z$max_std_dk, TT - 1),
  	                                 mean  = rep(0, TT - 1),
  	                                 corr  = varvec)[1]
  	}
	  if (TT == 2) {
	    quant <- stats::qnorm(1 - alpha)
	    pvalue <- 1-stats::pnorm(z$max_std_dk)
	  }
	  if (z$max_std_dk > quant) {
	    rejind <- 1
	    khat <- z$khat
	  }
	}

	result <- list(reject  = rejind,
	               estcp   = khat,
	               tstat   = z$max_std_dk,
	               pvalue  = pvalue,
	               corrmat = corr_mat_full)
	return(result)
}


