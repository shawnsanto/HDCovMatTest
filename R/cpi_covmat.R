
#' Estimate multiple change points among high-dimensional covariance matrices
#'
#' This function estimates multiple change points among covariance
#' matrices for high-dimensional longitudinal data. The change points are
#' identified using the testing procedure and binary segmentation approach
#' proposed in Zhong, Li, and Santo (2018).
#'
#' The methodology has been proved, tested, and shown to work well for
#' high-dimensional longtidudinal data where both \eqn{p} and \eqn{n} diverge,
#' but no explicit relationship between \eqn{p} and \eqn{n} is needed.
#' Furthermore, the current methodology works for longtudinal data with a small
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
#'   \eqn{TT \le 30} so the asymptotics in Zhong, Li, and Santo (2017)
#'   can be applied.
#' @param alpha The type I error of the homogeniety test. The nominal level for
#'   the quantile is computed as \eqn{1-alpha}. Suggested values for alpha
#'   include 0.01 (default) and 0.05.
#'
#' @return The function returns the estimated chane point(s), corresponding
#'     test statistic value(s), corresponding p-value(s), and a table that
#'     provides an identification for which time points have a homogeneous
#'     covariance structure. The output is a list.
#'   \describe{
#'     \item{$change_points}{The estimated change points. Order is based
#'       on the algorithm's binary segmentation approach.}
#'     \item{$teststats}{The test statistic(s) corresponding to the estimated
#'       change point(s).}
#'     \item{$pvalues}{The p-value(s) corresponding to the estimated
#'       change point(s).}
#'     \item{$covariance_id}{A table that indicates which covariance matrices
#'       are homogeneous given the estimated change point(s). For example, when
#'       TT = 5, a single change point identified at time 3 implies the
#'       covariance matrices for times 1, 2, and 3 are equal, but they are
#'       different from the covariance matrices that are equal at time points
#'       4 and 5.}
#'     \item{$note}{A comment that explains $covariance_id.}
#'       }
#'
#' @references \emph{Zhong, Li, and Santo (2018). Homogeneity tests of covariance
#'   matrices with high-dimensional longitudinal data. Biometrika.}
#'
#' @author \strong{Maintainer}: Shawn Santo \email{santosha@@stt.msu.edu}
#'
#'   Authors:
#'     \itemize{
#'       \item Ping-Shou Zhong
#'       \item Runze Li
#'       \item Shawn Santo
#'     }
#'
#' @export
#'
#' @examples
#' # A change point identification example with a change points at times 2 and 4
#'
#' # Set parameters
#' p <- 30; n <- 10; TT <- 5
#' delta <- 0.85
#' m <- p+20; L <- 3; k0 <- 2; k1 <- 4; w <- 0.2
#'
#' # Generate data
#' Gamma1 <- Gamma2 <- Gamma3 <- matrix(0, p, m * L)
#' y <- array(0, c(p, n, TT))
#' set.seed(928)
#'
#' for (i in 1:p){
#'   for (j in 1:p){
#'     dij <- abs(i - j)
#'
#'     if (dij < (p * w)){
#'       Gamma1[i, j] <- (dij + 1) ^ (-2)
#'       Gamma2[i, j] <- (dij + 1 + delta) ^ (-2)
#'       Gamma3[i, j] <- (dij + 1 + 2 * delta) ^ (-2)
#'     }
#'   }
#' }
#'
#' Z <- matrix(rnorm(m * (TT + L - 1) * n), m * (TT + L - 1), n)
#'
#' for (t in 1:k0){
#'   y[, , t] <- Gamma1 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
#' }
#' for (t in (k0 + 1):k1){
#'   y[, , t] <- Gamma2 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
#' }
#' for (t in (k1 + 1):TT){
#'   y[, , t] <- Gamma3 %*% Z[((t - 1) * m + 1):((t + L - 1) * m), ]
#' }
#'
#' cpi_covmat(y, n, p, TT, alpha = 0.01)


cpi_covmat <- function(y, n, p, TT, alpha = 0.01){

  multiCPs <- function(Yvec, n, p, TT_length, start, end){
    rejind <- 0
    khat <- 0

    if (TT_length > 1){
      khat <- numeric(1)
      stdDk <- numeric(TT_length - 1)
      maxstdDk <- numeric(1)
      CorrMat <- numeric((TT_length - 1) * TT_length / 2)

      nYvec <- Yvec[(start * n * p + 1):(end * n * p)]
      storage.mode(nYvec) <- 'double'

      z <- .C('testandCP4r',
             as.double(nYvec),
             as.integer(n),
             as.integer(p),
             as.integer(TT_length),
             khat = as.integer(khat),
             stdDk = as.double(stdDk),
             maxstdDk = as.double(maxstdDk),
             CorrMat = as.double(CorrMat))

      if (TT_length > 2){
        offdiag <- z$CorrMat[-c((2 * (TT_length - 1) - (0:(TT_length - 2)) + 1) *
                                  (0:(TT_length - 2)) / 2 + 1)]
        Varkvec <- matrix(1, TT_length - 1, TT_length - 1)
        Varkvec[lower.tri(Varkvec)] <- offdiag
        Varkvec[upper.tri(Varkvec)] <- t(Varkvec)[upper.tri(t(Varkvec))]
        cmat <- Varkvec

        qTn <- mvtnorm::qmvnorm(1 - alpha, corr = Varkvec)$quantile
        pvalue <- 1-mvtnorm::pmvnorm(lower = -Inf,
                                     upper = rep(z$maxstdDk, TT_length - 1),
                                     mean = rep(0, TT_length - 1),
                                     corr = Varkvec)[1]
      }

      if (TT_length == 2){
        qTn <- stats::qnorm(1 - alpha)
        pvalue <- 1-stats::pnorm(z$maxstdDk)
      }

      if (z$maxstdDk>qTn){
        rejind <- 1
        khat <- z$khat
      }
    }
    return(c(rejind, khat, z$maxstdDk, pvalue))
  }

  binarySearch <- function(Yvec, n, p, start, end){

    segmet_list <- c(start, end)
    found_list <- NULL
    test_stats <- NULL
    p_values <- NULL
    while (length(segmet_list) > 0){
      s1 <- segmet_list[1]
      e1 <- segmet_list[2]
      TT_len <- e1-s1

      if (TT_len > 1){
        cpvec <- multiCPs(Yvec, n, p, TT_len, s1, e1)

        if (cpvec[1] == 1) {
          k1 <- s1 + cpvec[2]
          segmet_list <- c(segmet_list, s1, k1, k1, e1)
          found_list <- c(found_list, k1)
          test_stats <- c(test_stats, cpvec[3])
          p_values <- c(p_values, cpvec[4])
        }
       }
      segmet_list <- segmet_list[-c(1, 2)]
    }
    results <- rbind(found_list, test_stats, p_values)
    rownames(results) <- NULL
    return(results)
  }

  Yvec <- matrix(y, 1, (p * n * TT))

  cpts <- binarySearch(Yvec, n, p, 0, TT)

  if (!is.null(cpts)){
    time_points <- c(1:TT)
    covariance_id <- NULL

    myletters <- expand.grid(LETTERS, LETTERS)
    myletters <- do.call(paste0, myletters)
    myletters_vec <- c(LETTERS, as.vector(t(matrix(myletters, 26, 26))))

    cp_identifier <- c(0, sort(cpts[1, ]), TT)

    for (i in 1:(length(cp_identifier) - 1)){
      covariance_id <- c(covariance_id,
                         rep(myletters_vec[i],
                             (cp_identifier[i + 1] - cp_identifier[i])))
    }

    note <- "The same character indicates the same population covariance matrix."

    covariance_id_table <- as.data.frame(rbind(time_points, covariance_id))
    rownames(covariance_id_table) <- c("Time Point", "Covariance ID")
    names(covariance_id_table) <- NULL

    return(list(change_points = cpts[1, ],
                teststats = cpts[2, ],
                pvalues = cpts[3, ],
                covariance_id = covariance_id_table,
                note = note))
  }else{
    return(cat("No change points found!"))
  }
}
