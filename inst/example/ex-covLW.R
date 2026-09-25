set.seed(123)
X <- matrix(rnorm(100), nrow = 20, ncol = 5)

## Linear shrinkage; center X internally and use an effective sample size of N - 1
Sigma_linear <- covLW(X)

# Linear shrinkage; population mean is assumed to be known and equal to zero
Sigma_linear_zeromean <- covLW(X, k = 0, method = "linear")

## Linear-inverse shrinkage; internal centering
Sigma_lis <- covLW(X, method = "lis")

# Linear-inverse shrinkage; X has already been centered as one class
X_centered <- scale(X, center = TRUE, scale = FALSE)
Sigma_lis_centered <- covLW(X_centered, k = 1, method = "lis")
