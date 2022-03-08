The goal of this project is to build an R package for local polynomial regression
with the option for explicit boundary adjustment via boundary kernels when local constant regression (Nadaraya-Watson) is applied.
The latter is not available in R.

Intended features of the package include:
* Computation of effective kernels, i.e. effectively assigned weights in the smoothing process.
* Estimates for the first and second derivative of the regression function.
* Different (compactly supported) kernel functions.
* Different boundary kernels.
* Bandwidth selection via cross-validation.
* Computation of asymptotic confidence intervals.

Moreover, we
* provide a motivation for the package; (to be done)
* provide the theoretical background for the package; (to be done)
* showcase the practical use of the package through a real-data application.
