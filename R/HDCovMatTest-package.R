#' HDCovMatTest: Homogeneity tests of covariance matrices with high-dimensional longitudinal data
#'
#' An implementation for testing the homogeneity of covariance matrices in
#' high-dimensional longitudinal data with temporospatial dependence. The null
#' hypothesis is that all covariance matrices are equal at each time point.
#' If the null hypothesis is rejected, estimates of the change points are
#' provided. A binary segmentation approach is applied to estimate multiple
#' change points.
#'
#' \tabular{ll}{
#'   Package: \tab HDCovMatTest\cr
#'   Type: \tab package\cr
#'   Version: \tab 1.0.5\cr
#'   Date: \tab 2020-12-04\cr
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
#' @name HDCovMatTest-package
#' @docType package
NULL
