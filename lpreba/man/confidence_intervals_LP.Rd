\name{confidence_intervals_LP}
\alias{confidence_intervals_LP}
\title{Asymptotic confidence intervals for the local polynomial estimator}
\usage{
confidence_intervals_LP(x, X, Y, kernel = epanechnikov, bw, degree = 1L, alpha = 0.05)
}
\arguments{
\item{x}{Evaluation points (vector).}
\item{X}{Data for the regressor (vector).}
\item{Y}{Data for the regressand (vector).}
\item{kernel}{Kernel (function). Default is \code{epanechnikov}.}
\item{bw}{Bandwidth (scalar).}
\item{degree}{Degree of the locally fitted polynomial (integer). Default is \code{1L}.}
\item{alpha}{Significance level (scalar). Default is \code{0.05}.}
}
\value{
List containing:
\item{confidence_intervals_lower}{Lower points for the (1 - alpha) asymptotic confidence intervals (vector).}
\item{confidence_intervals_upper}{Upper points for the (1 - alpha) asymptotic confidence intervals (vector).}
}
\description{
Pointwise asymptotic confidence intervals of arbitrary significance level
for the local polynomial estimator.
}
\details{
The confidence intervals account for the variance of the estimator but not
its bias, i.e. they are not bias-corrected.

The standard errors for the construction of the asymptotic confidence intervals
are computed according to Hansen (2022, Section 19.16).
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

output_confidence_intervals_LP <- confidence_intervals_LP(x = x, X = X, Y = Y, bw = bw,
							  alpha = 0.05)
confidence_intervals_lower <- output_confidence_intervals_LP$confidence_intervals_lower
confidence_intervals_upper <- output_confidence_intervals_LP$confidence_intervals_upper
}
\references{
Hansen, B. (2022).
\emph{Econometrics}.
New Jersey: Princeton University Press.
\url{https://www.ssc.wisc.edu/~bhansen/econometrics/Econometrics.pdf}.
Revised: 8/18/2021.
}