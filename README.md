# lpreba: <ins>L</ins>ocal <ins>p</ins>olynomial <ins>r</ins>egression with option for <ins>e</ins>xplicit <ins>b</ins>oundary <ins>a</ins>djustment

## Description

The R package implements the local polynomial regression estimator of arbitrary degree. <br>
Estimates for the regression function (conditional mean), its first and second derivative can be obtained.
Moreover, computation of effective kernels (i.e. effectively assigned weights in the kernel smoothing process) is provided.
Different compactly supported kernels (Uniform, Epanechnikov, etc.) are available. <br>
For local constant regression (Nadaraya-Watson), explicit boundary adjustment via boundary kernels can be conducted.
Different boundary kernels are available. <br>
Additional functionality includes bandwidth selection via cross-validation and the computation of asymptotic confidence intervals.

## Functions

| Name                              | Description                                                        |
|-----------------------------------|--------------------------------------------------------------------|
| ``` LP() ```                      | Local polynomial estimator                                         |
| ``` NW_boundary() ```             | Boundary-adjusted Nadaraya-Watson estimator                        |
| ``` CV_error_fun() ```            | Leave-one-out cross-validation (LOOCV) error                       |
| ``` bw_CV_fun() ```               | CV optimal bandwidth                                               |
| ``` confidence_intervals_LP() ``` | Asymptotic confidence intervals for the local polynomial estimator |

## How to get started

```r
install.packages("devtools")

devtools::install_github("svjaco/epp_project/lpreba") # Install package
library(lpreba) # Load and attach package

?lpreba # Help file for package
# Use e.g. ?LP to access the documentation for the LP function
```

## Examples

```r
# Example data

set.seed(123) # Seed for reproducibility

m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Data for regressor
m_X <- m_fun(X) # True values of regression function
epsilon <- rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Data for regressand
bw <- 0.2 # Bandwidth
x <- seq(0, 1, length.out = n/2) # Evaluation points
bw_grid <- seq(0.05, 0.2, length.out = 200) # Bandwidth grid

# Example for LP function
output_LP <- LP(x = x, X = X, Y = Y, kernel = epanechnikov, bw = bw, degree = 1L)
estimates_LP <- output_LP$estimates
effective_kernels_LP <- output_LP$effective_kernels
slopes_LP <- output_LP$slopes
curvatures_LP <- output_LP$curvatures # Yields only NAs since degree >= 2L is required

# Example for NW_boundary function
output_NW_boundary <- NW_boundary(x = x, X = X, Y = Y, bw = bw,
                                  kernel_interior = epanechnikov, kernel_left = epanechnikov_left,
                                  boundary_left = 0, boundary_right = 1)
estimates_NW_boundary <- output_NW_boundary$estimates
effective_kernels_NW_boundary <- output_NW_boundary$effective_kernels

# Example for CV_error_fun function

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
             
# Example for bw_CV_fun function

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
          
# Example for confidence_intervals_LP function
output_confidence_intervals_LP <- confidence_intervals_LP(x = x, X = X, Y = Y, bw = bw,
							  alpha = 0.05)
confidence_intervals_lower <- output_confidence_intervals_LP$confidence_intervals_lower
confidence_intervals_upper <- output_confidence_intervals_LP$confidence_intervals_upper
```

## Application

To showcase the practical use of the ``lpreba`` package,
we consider the nonparametric estimation of treatment effects in the regression discontinuity design (RDD).

To build the final manuscript ``application.pdf``, we use the build system [pytask](https://github.com/pytask-dev/pytask).
First, an environment according to ``environment.yml`` needs to be created.
Using conda this can be achieved by running in the terminal:

```zsh
# cd into root of project
$ conda env create -f environment.yml
$ conda activate epp_project_sven_jacobs
$ pip install -e .
```

To install the ``lpreba`` package, open an R terminal by typing ``R`` in the terminal and run:

```r
devtools::install_github("svjaco/epp_project/lpreba")
```

The last step is:

```zsh
# cd into root of project
$ pytask
```

The file ``application.pdf`` is then located in the ``bld`` folder.

---

[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/svjaco/epp_project/blob/main/LICENSE)
[![R-CMD-check](https://github.com/svjaco/epp_project/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/svjaco/epp_project/actions/workflows/R-CMD-check.yaml)
