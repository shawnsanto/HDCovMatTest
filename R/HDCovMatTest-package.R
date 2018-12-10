#' tecoma: Test of high-dimensional covariance matrices with longitudinal data
#'
#' An implementation for testing the homogeneity of covariance matrices in
#' high-dimensional longitudinal data with temporospatial dependence. The null
#' hypothesis is that all covariance matrices are equal at each time point.
#' If the null hypothesis is rejected, estimates of the change points are provided.
#' A binary segmentation approach is applied to estimate multiple change points.
#'
#' \tabular{ll}{
#'   Package: \tab HDCovMatTest\cr
#'   Type: \tab package\cr
#'   Version: \tab 1.0.0\cr
#'   Date: \tab 2018-09-01\cr
#'   License: \tab GPL-2\cr
#'   }
#'
#' @section Functions:
#'
#' \itemize{
#'   \item test_covmat
#'   \item cpi_covmat
#' }
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
#' @references \emph{Zhong, Li, and Santo (2018). Homogeneity tests of covariance
#'   matrices with high-dimensional longitudinal data. Biometrika.}
#'
#' @name HDCovMatTest-package
#' @docType package
NULL
