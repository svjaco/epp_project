# Data for automated testing used by test files

set.seed(123) # Seed for reproducibility

n <- 100 # Sample size
X <- seq(0, 1, length.out = n) # Data for regressor
m_X <- sin(2*pi*X) # True values of regression function
epsilon <- stats::rnorm(n, sd = 0.25) # Error term
Y <- m_X + epsilon # Data for regressand
bw <- 0.2 # Bandwidth
df <- data.frame(Y, X) # Data frame, required for locpol
bw_grid <- seq(0.05, 0.2, length.out = 1000) # Bandwidth grid
interval <- c(min(bw_grid), max(bw_grid)) # Interval limits of bandwidth grid, required for locpol