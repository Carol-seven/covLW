# Ledoit-Wolf Linear Shrinkage Covariance Estimator

Estimates a covariance matrix by shrinking the sample covariance matrix
toward a scaled identity matrix (Ledoit and Wolf 2004) .

## Usage

``` r
covLW(X, k = -1)
```

## Arguments

- X:

  A numeric matrix containing the data. Rows represent observations and
  columns represent variables. `X` must contain at least two rows and
  one column, and all entries must be finite.

- k:

  A numeric scalar controlling centering and the effective sample size.

  - `k < 0`: The columns of `X` are centered internally, `k` is set to
    1, and the effective sample size is `nrow(X) - 1`.

  - `k = 0`: No centering is performed and the effective sample size is
    `nrow(X)`; this case is appropriate when the population mean is
    known to be zero.

  - `k >= 1`: `X` is assumed to contain `k` classes that have already
    been centered separately, and the effective sample size is
    `nrow(X) - k`.

  Nonnegative values of `k` must be integers smaller than `nrow(X)`.

## Value

A symmetric numeric \\p \times p\\ matrix containing the Ledoit-Wolf
linear shrinkage estimate of the covariance matrix.

## Details

Let \\N\\ be the number of rows of `X`, \\p\\ its number of columns, and
\\n = N - k\\ the effective sample size after the treatment of `k`. The
sample covariance matrix is \\S = X^\top X / n\\. The shrinkage target
is the scaled identity matrix \\\widehat{m} I_p\\, where \\\widehat{m} =
\operatorname{tr}(S) / p\\ is the average sample variance.

Using the normalized squared Frobenius norm, the estimated squared
distance between the sample covariance matrix and the target is
\$\$\widehat{d}^2 = \frac{1}{p} \lVert S - \widehat{m}I_p
\rVert_F^2.\$\$

The estimator of the sampling error is \$\$\overline{b}^2 =
\frac{\sum\_{i=1}^{N} \lVert x_i x_i^\top \rVert_F^2 - n \lVert S
\rVert_F^2}{p n^2},\$\$ where \\x_i^\top\\ is row \\i\\ of `X`. It is
truncated to \\\widehat{b}^2 = \min\\\max(\overline{b}^2,
0),\widehat{d}^2\\\\, and \\\widehat{a}^2 = \widehat{d}^2 -
\widehat{b}^2\\.

The resulting covariance estimator is \$\$\widehat{S}^{\ast} =
\frac{\widehat{b}^2}{\widehat{d}^2}\widehat{m}I_p +
\frac{\widehat{a}^2}{\widehat{d}^2}S.\$\$

If \\\widehat{d}^2 = 0\\, the function returns the target matrix
directly.

In the formulas above, \\X\\ denotes the matrix after any centering
performed by the function.

The function does not perform class-specific centering when `k >= 1`.
Such centering must be completed before calling the function.

## References

Ledoit O, Wolf M (2004). “A Well-Conditioned Estimator for
Large-Dimensional Covariance Matrices.” *Journal of Multivariate
Analysis*, **88**(2), 365–411.
[doi:10.1016/S0047-259X(03)00096-4](https://doi.org/10.1016/S0047-259X%2803%2900096-4)
.

## Examples

``` r
set.seed(123)
X <- matrix(rnorm(100), nrow = 20, ncol = 5)

# Center X internally and use an effective sample size of N - 1
Sigma_hat <- covLW(X)

# X has already been centered as one class
X_centered <- scale(X, center = TRUE, scale = FALSE)
Sigma_hat_centered <- covLW(X_centered, k = 1)

# Population mean is assumed to be known and equal to zero
Sigma_hat_zeromean <- covLW(X, k = 0)
```
