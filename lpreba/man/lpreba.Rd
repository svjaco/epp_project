\docType{package}
\name{lpreba}
\alias{lpreba}
\title{lpreba: Local polynomial regression with option for explicit boundary adjustment}
\description{
The package implements the local polynomial regression estimator
of arbitrary degree. Estimates for the regression function (conditional mean),
its first and second derivative can be obtained. Moreover, computation of
effective kernels (i.e. effectively assigned weights in the kernel smoothing
process) is provided. Different compactly supported kernels (Uniform,
Epanechnikov, etc.) are available. For local constant regression
(Nadaraya-Watson), explicit boundary adjustment via boundary kernels can be
conducted. Different boundary kernels are available. Additional functionality
includes bandwidth selection via cross-validation and the computation of
asymptotic confidence intervals.
}
\section{Functions}{
\tabular{ll}{
\code{LP} \tab Local polynomial estimator. \cr
\code{NW_boundary} \tab Boundary-adjusted Nadaraya-Watson estimator. \cr
\code{CV_error_fun} \tab Leave-one-out cross-validation (LOOCV) error. \cr
\code{bw_CV_fun} \tab CV optimal bandwidth. \cr
\code{confidence_intervals_LP} \tab Asymptotic confidence intervals for the local polynomial estimator.
}
}
\author{
Sven Jacobs.
}