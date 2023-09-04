\name{LP}
\alias{LP}
\title{Local polynomial estimator}
\usage{
LP(x, X, Y, kernel = epanechnikov, bw, degree = 1L)
}
\arguments{
\item{x}{Evaluation points (vector).}
\item{X}{Data for the regressor (vector).}
\item{Y}{Data for the regressand (vector).}
\item{kernel}{Kernel (function). Default is \code{epanechnikov}.}
\item{bw}{Bandwidth (scalar).}
\item{degree}{Degree of the locally fitted polynomial (integer). Default is \code{1L}.}
}
\value{
List containing:
\item{estimates}{Estimates for the regression function at the evaluation points (vector).}
\item{effective_kernels}{Effective kernels at the evaluation points (matrix).}
\item{slopes}{Estimates for the first derivative (slope) of the regression function at the evaluation points (vector).}
\item{curvatures}{Estimates for the second derivative (curvature) of the regression function at the evaluation points (vector).}
}
\description{
Local polynomial regression estimator of arbitrary degree.
}
\details{
The degree of the local polynomial estimator can be chosen. Two special
cases are the Nadaraya-Watson estimator (\code{degree = 0L}) and the
local linear estimator (\code{degree = 1L}).

The most important compactly supported (on the interval [-1, 1]) kernels
are available: \code{uniform}, \code{triangular}, \code{epanechnikov},
\code{biweight}, \code{triweight}, \code{tricube}, \code{cosine}.

The effective kernel at an evaluation point is the set of effectively
assigned weights in the smoothing process (see Hastie and Loader, 1993).

To obtain estimates for the first derivative of the regression function,
the degree has to be at least \code{1L}.
To obtain estimates for the second derivative of the regression function,
the degree has to be at least \code{2L}.
}
\examples{
m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Data for the regressor
m_X <- m_fun(X) # True values of regression function
epsilon <- rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Data for the regressand
bw <- 0.2 # Bandwidth
x <- seq(0, 1, length.out = n/2) # Evaluation points

output_LP <- LP(x = x, X = X, Y = Y, kernel = epanechnikov, bw = bw, degree = 1L)
estimates_LP <- output_LP$estimates
effective_kernels_LP <- output_LP$effective_kernels
slopes_LP <- output_LP$slopes
curvatures_LP <- output_LP$curvatures # Yields only NAs since degree >= 2L is required
}
\references{
Hastie, T. and C. Loader (1993).
``Local regression: Automatic kernel carpentry''.
\emph{Statistical Science} 8 (2), pp. 120--129.
}