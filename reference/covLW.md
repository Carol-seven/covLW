# Ledoit-Wolf Covariance Matrix Estimation

Estimates a covariance matrix using a selected Ledoit-Wolf linear or
nonlinear shrinkage method.

## Usage

``` r
covLW(X, k = -1, method = "linear")
```

## Arguments

- X:

  A numeric matrix containing the data. Rows represent observations and
  columns represent variables. `X` must contain at least two rows and
  one column, and all entries must be finite.

- k:

  A numeric scalar (default = -1) controlling centering and the
  effective sample size.

  - `k < 0`: The columns of `X` are centered internally, `k` is set to
    1, and the effective sample size is `nrow(X) - 1`.

  - `k = 0`: No centering is performed and the effective sample size is
    `nrow(X)`; this case is appropriate when the population mean is
    known to be zero.

  - `k >= 1`: `X` is assumed to contain `k` classes that have already
    been centered separately, and the effective sample size is
    `nrow(X) - k`.

  Nonnegative values of `k` must be integers smaller than `nrow(X)`.

- method:

  A character string (default = "linear") specifying the method for
  covariance matrix estimation:

  1.  `"linear"`: Linear shrinkage toward a scaled identity matrix
      (Ledoit and Wolf 2004) .

  2.  `"lis"`: Linear-inverse shrinkage, nonlinear shrinkage derived
      under Stein's loss (Ledoit and Wolf 2022) .

## Value

A symmetric numeric matrix containing the selected Ledoit-Wolf
covariance matrix estimate.

## References

Ledoit O, Wolf M (2004). “A Well-Conditioned Estimator for
Large-Dimensional Covariance Matrices.” *Journal of Multivariate
Analysis*, **88**(2), 365–411.
[doi:10.1016/S0047-259X(03)00096-4](https://doi.org/10.1016/S0047-259X%2803%2900096-4)
.  
  
Ledoit O, Wolf M (2022). “Quadratic Shrinkage for Large Covariance
Matrices.” *Bernoulli*, **28**(3), 1519–1547.
[doi:10.3150/20-BEJ1315](https://doi.org/10.3150/20-BEJ1315) .

## Examples

``` r
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
```
