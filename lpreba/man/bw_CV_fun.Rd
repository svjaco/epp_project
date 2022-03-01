\name{bw_CV_fun}
\alias{bw_CV_fun}
\title{CV optimal bandwidth}
\usage{
bw_CV_fun(X, Y, kernel = epanechnikov, bw_grid, degree = 1L,
          kernel_left = epanechnikov_left,
          boundary_left = NA, boundary_right = NA,
          boundary_adjustment = FALSE,
          plot = TRUE)
}
\arguments{
\item{X}{Data for the regressor (vector).}
\item{Y}{Data for the regressand (vector).}
\item{kernel}{Kernel (function). Default is \code{epanechnikov}.}
\item{bw_grid}{Bandwidth grid to find global minimum of CV error (vector).}
\item{degree}{Degree of the locally fitted polynomial (integer). Default is \code{1L}.}
\item{kernel_left}{Left boundary kernels (function). Default is \code{epanechnikov_left}.}
\item{boundary_left}{Lower boundary of the support of X (scalar). Default is \code{NA}.}
\item{boundary_right}{Upper boundary of the support of X (scalar). Default is \code{NA}.}
\item{boundary_adjustment}{Explicit boundary adjustment (boolean). Default is \code{FALSE}.}
\item{plot}{Plot CV error over bandwidth grid (boolean). Default is \code{TRUE}.}
}
\value{
CV optimal bandwidth (scalar) and plot of CV errors if \code{plot = TRUE}.
}
\description{
Selection of the bandwidth that yields the smallest leave-one-out
cross-validation (LOOCV) error over a supplied grid.
}
\examples{
m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Data for the regressor
m_X <- m_fun(X) # True values of regression function
epsilon <- rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Data for the regressand
bw_grid <- seq(0.05, 0.2, length.out = 200) # Bandwidth grid

# Local polynomial estimator
bw_CV_fun(X = X, Y = Y, kernel = epanechnikov, bw_grid = bw_grid, degree = 1L,
          kernel_left = epanechnikov_left,
          boundary_left = NA, boundary_right = NA,
          boundary_adjustment = FALSE,
          plot = TRUE)

# Boundary-adjusted Nadaraya-Watson estimator
bw_CV_fun(X = X, Y = Y, kernel = epanechnikov, bw_grid = bw_grid, degree = 0L,
          kernel_left = epanechnikov_left,
          boundary_left = 0, boundary_right = 1,
          boundary_adjustment = TRUE,
          plot = TRUE)

}
\seealso{
\code{\link{CV_error_fun}}
}