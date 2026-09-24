#' Ledoit-Wolf Linear Shrinkage Covariance Estimator
#'
#' Estimates a covariance matrix by shrinking the sample covariance matrix
#' toward a scaled identity matrix \insertCite{ledoit2004well}{covLW}.
#'
#' @param X A numeric matrix containing the data.
#' Rows represent observations and columns represent variables.
#' \code{X} must contain at least two rows and one column, and all entries
#' must be finite.
#'
#' @param k A numeric scalar controlling centering and the effective sample size.
#' \itemize{
#' \item \code{k < 0}: The columns of \code{X} are centered internally,
#' \code{k} is set to 1, and the effective sample size is \code{nrow(X) - 1}.
#' \item \code{k = 0}: No centering is performed and the effective sample size
#' is \code{nrow(X)}; this case is appropriate when the population mean is known
#' to be zero.
#' \item \code{k >= 1}: \code{X} is assumed to contain \code{k} classes that
#' have already been centered separately, and the effective sample size is
#' \code{nrow(X) - k}.
#' }
#' Nonnegative values of \code{k} must be integers smaller than \code{nrow(X)}.
#'
#' @return A symmetric numeric \eqn{p \times p} matrix containing
#' the Ledoit-Wolf linear shrinkage estimate of the covariance matrix.
#'
#' @details
#' Let \eqn{N} be the number of rows of \code{X}, \eqn{p} its number of columns,
#' and \eqn{n = N - k} the effective sample size after the treatment of \code{k}.
#' The sample covariance matrix is \eqn{S = X^\top X / n}.
#' The shrinkage target is the scaled identity matrix \eqn{\widehat{m} I_p},
#' where \eqn{\widehat{m} = \operatorname{tr}(S) / p} is the average sample variance.
#'
#' Using the normalized squared Frobenius norm, the estimated squared distance
#' between the sample covariance matrix and the target is
#' \deqn{\widehat{d}^2 = \frac{1}{p} \lVert S - \widehat{m}I_p \rVert_F^2.}
#'
#' The estimator of the sampling error is
#' \deqn{\overline{b}^2 = \frac{\sum_{i=1}^{N} \lVert x_i x_i^\top \rVert_F^2
#' - n \lVert S \rVert_F^2}{p n^2},}
#' where \eqn{x_i^\top} is row \eqn{i} of \code{X}. It is truncated to
#' \eqn{\widehat{b}^2 = \min\{\max(\overline{b}^2, 0),\widehat{d}^2\}},
#' and \eqn{\widehat{a}^2 = \widehat{d}^2 - \widehat{b}^2}.
#'
#' The resulting covariance estimator is
#' \deqn{\widehat{S}^{\ast} = \frac{\widehat{b}^2}{\widehat{d}^2}\widehat{m}I_p
#' + \frac{\widehat{a}^2}{\widehat{d}^2}S.}
#'
#' If \eqn{\widehat{d}^2 = 0}, the function returns the target matrix directly.
#'
#' In the formulas above, \eqn{X} denotes the matrix after any centering
#' performed by the function.
#'
#' The function does not perform class-specific centering when \code{k >= 1}.
#' Such centering must be completed before calling the function.
#'
#' @references
#' \insertAllCited{}
#'
#' @example
#' inst/example/ex-covest.R
#'
#' @importFrom Rdpack reprompt
#'
#' @export

covest <- function(X, k = -1) {

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }

  if (nrow(X) < 2L || ncol(X) < 1L) {
    stop("`X` must have at least two rows and one column.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite, non-missing values.", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) {
    stop("`k` must be a finite numeric scalar.", call. = FALSE)
  }

  N <- nrow(X)
  p <- ncol(X)

  if (k < 0) {
    X <- scale(X, center = TRUE, scale = FALSE)
    k <- 1L
  } else if (k != floor(k) || k >= N) {
    stop("`k` must be a nonnegative integer smaller than `nrow(X)`.", call. = FALSE)
  }

  n <- N - k

  S <- crossprod(X) / n

  m <- sum(diag(S)) / p
  target <- m * diag(p)

  d2 <- sum((S - target)^2) / p

  b2_overline <- (sum(rowSums(X^2)^2) - n * sum(S^2)) / (p * n^2)
  b2 <- min(max(b2_overline, 0), d2)

  a2 <- d2 - b2

  if (d2 == 0) {
    Sstar <- target
  } else {
    Sstar <- (b2 / d2) * target + (a2 / d2) * S
  }

  return(Sstar)
}
