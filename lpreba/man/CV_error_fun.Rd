\name{CV_error_fun}
\alias{CV_error_fun}
\title{Leave-one-out cross-validation (LOOCV) error}
\usage{
CV_error_fun(X, Y, kernel = epanechnikov, bw, degree = 1L,
             kernel_left = epanechnikov_left,
             boundary_left = NA, boundary_right = NA,
             boundary_adjustment = FALSE)
}
\arguments{
\item{X}{Data for the regressor (vector).}
\item{Y}{Data for the regressand (vector).}
\item{kernel}{Kernel (function). Default is \code{epanechnikov}.}
\item{bw}{Bandwidth (scalar).}
\item{degree}{Degree of the locally fitted polynomial (integer). Default is \code{1L}.}
\item{kernel_left}{Left boundary kernels (function). Default is \code{epanechnikov_left}.}
\item{boundary_left}{Lower boundary of the support of X (scalar). Default is \code{NA}.}
\item{boundary_right}{Upper boundary of the support of X (scalar). Default is \code{NA}.}
\item{boundary_adjustment}{Explicit boundary adjustment (boolean). Default is \code{FALSE}.}
}
\value{
LOOCV error for bandwidth \code{bw} (scalar).
}
\description{
Leave-one-out cross-validation (LOOCV) error for a given bandwidth, either
for the local polynomial estimator or the boundary-adjusted
Nadaraya-Watson estimator.
}
\details{
The LOOCV error function is the mean squared prediction errors.
}
\examples{
m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Data for the regressor
m_X <- m_fun(X) # True values of regression function
epsilon <- rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Data for the regressand
bw <- 0.2 # Bandwidth

# Local polynomial estimator
CV_error_fun(X = X, Y = Y, kernel = epanechnikov, bw = bw, degree = 1L,
             kernel_left = epanechnikov_left,
             boundary_left = NA, boundary_right = NA,
             boundary_adjustment = FALSE)

# Boundary-adjusted Nadaraya-Watson estimator
CV_error_fun(X = X, Y = Y, kernel = epanechnikov, bw = bw, degree = 0L,
             kernel_left = epanechnikov_left,
             boundary_left = 0, boundary_right = 1,
             boundary_adjustment = TRUE)
}