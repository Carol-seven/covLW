set.seed(123)
X <- matrix(rnorm(100), nrow = 20, ncol = 5)

# Center X internally and use an effective sample size of N - 1
Sigma_hat <- covest(X)

# X has already been centered as one class
X_centered <- scale(X, center = TRUE, scale = FALSE)
Sigma_hat_centered <- covest(X_centered, k = 1)

# Population mean is assumed to be known and equal to zero
Sigma_hat_zeromean <- covest(X, k = 0)
