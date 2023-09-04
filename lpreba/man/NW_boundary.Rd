\name{NW_boundary}
\alias{NW_boundary}
\title{Boundary-adjusted Nadaraya-Watson estimator}
\usage{
NW_boundary(x, X, Y, bw,
            kernel_interior = epanechnikov, kernel_left = epanechnikov_left,
            boundary_left, boundary_right)
}
\arguments{
\item{x}{Evaluation points (vector).}
\item{X}{Data for the regressor (vector).}
\item{Y}{Data for the regressand (vector).}
\item{bw}{Bandwidth (scalar).}
\item{kernel_interior}{Kernel used in the interior (function). Default is \code{epanechnikov}.}
\item{kernel_left}{Left boundary kernels (function). Default is \code{epanechnikov_left}.}
\item{boundary_left}{Lower boundary of the support of X (scalar).}
\item{boundary_right}{Upper boundary of the support of X (scalar).}
}
\value{
List containing:
\item{estimates}{Estimates for the regression function at the evaluation points (vector).}
\item{effective_kernels}{Effective kernels at the evaluation points (matrix).}
}
\description{
Explicit boundary adjustment via boundary kernels for the
Nadaraya-Watson estimator (local constant regression).
}
\details{
When applying the Nadaraya-Watson estimator with a compactly supported
(on the interval [-1, 1]) kernel to data of a regressor that is
compactly supported, boundary effects arise over the regions
[\code{boundary_left}, \code{boundary_left} + \code{bw}) and
(\code{boundary_right} - \code{bw}, \code{boundary_right}].
That is, the bias is of larger order compared to the interior.
The original order can be preserved by using special so-called boundary
kernels when applying Nadaraya-Watson at the boundaries
(e.g. Gasser and Müller, 1979).

``Smooth optimum boundary kernels'' (Müller, 1991, Table 1) are available:
\code{uniform_left}, \code{epanechnikov_left}, \code{biweight_left},
\code{triweight_left}.
Moreover, an alternative version of these kernels with greater asymptotic
efficiency at the cost of discontinuity at their endpoints
(Müller and Wang, 1994, Table 1) is available:
\code{uniform_alt_left}, \code{epanechnikov_alt_left}, \code{biweight_alt_left},
\code{triweight_alt_left}.

The effective kernel at an evaluation point is the set of effectively
assigned weights in the smoothing process (see Hastie and Loader, 1993).
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

output_NW_boundary <- NW_boundary(x = x, X = X, Y = Y, bw = bw,
                                  kernel_interior = epanechnikov, kernel_left = epanechnikov_left,
                                  boundary_left = 0, boundary_right = 1)
estimates_NW_boundary <- output_NW_boundary$estimates
effective_kernels_NW_boundary <- output_NW_boundary$effective_kernels
}
\references{
Gasser, T. and H.-G. Müller (1979).
``Kernel estimation of regression functions''.
In: \emph{Smoothing Techniques for Curve Estimation}.
Ed. by T. Gasser and M. Rosenblatt.
Berlin: Springer, pp. 23--68.

Müller, H.-G. (1991).
``Smooth optimum kernel estimators near endpoints''.
\emph{Biometrika} 78 (3), pp. 521--530.

Müller, H.-G. and J.-L. Wang (1994).
``Hazard rate estimation under random censoring with varying kernels and bandwidths''.
\emph{Biometrics} 50 (1), pp. 61--76.

Hastie, T. and C. Loader (1993).
``Local regression: Automatic kernel carpentry''.
\emph{Statistical Science} 8 (2), pp. 120--129.
}