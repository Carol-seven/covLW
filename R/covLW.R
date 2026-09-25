#' Ledoit-Wolf Covariance Matrix Estimation
#'
#' Estimates a covariance matrix using a selected Ledoit-Wolf linear or
#' nonlinear shrinkage method.
#'
#' @param X A numeric matrix containing the data.
#' Rows represent observations and columns represent variables.
#' \code{X} must contain at least two rows and one column, and all entries
#' must be finite.
#'
#' @param k A numeric scalar (default = -1) controlling centering and
#' the effective sample size.
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
#' @param method A character string (default = "linear") specifying the method
#' for covariance matrix estimation:
#' \enumerate{
#' \item \code{"linear"}: Linear shrinkage toward a scaled identity matrix
#' \insertCite{ledoit2004well}{covLW}.
#' \item \code{"lis"}: Linear-inverse shrinkage, nonlinear shrinkage derived
#' under Stein's loss \insertCite{ledoit2022quadratic}{covLW}.
#' }
#'
#' @return A symmetric numeric matrix containing the selected Ledoit-Wolf
#' covariance matrix estimate.
#'
#' @references
#' \insertAllCited{}
#'
#' @example
#' inst/example/ex-covLW.R
#'
#' @importFrom Rdpack reprompt
#'
#' @export

covLW <- function(X, k = -1, method = "linear") {

  methods <- c("linear", "lis")
  method <- match.arg(method, choices = methods)

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }

  N <- nrow(X)
  p <- ncol(X)

  if (N < 2L || p < 1L) {
    stop("`X` must have at least two rows and one column.", call. = FALSE)
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop("`X` must contain only finite, non-missing values.", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k)) {
    stop("`k` must be a finite numeric scalar.", call. = FALSE)
  }

  if (k < 0) {
    X <- scale(X, center = TRUE, scale = FALSE)
    k <- 1L
  } else if (k != floor(k) || k >= N) {
    stop("`k` must be a nonnegative integer smaller than `nrow(X)`.", call. = FALSE)
  }

  n <- N - k
  S <- crossprod(X) / n

  if (method == "linear") {

    m <- sum(diag(S)) / p
    target <- m * diag(p)

    d2 <- sum((S - target)^2) / p

    if (d2 == 0) {
      dimnames(target) <- dimnames(S)
      return(target)
    }

    b2_overline <- (sum(rowSums(X^2)^2) - n * sum(S^2)) / (p * n^2)
    b2 <- min(max(b2_overline, 0), d2)
    a2 <- d2 - b2

    result <- (b2 / d2) * target + (a2 / d2) * S

  } else if (method == "lis") {

    r <- p / n

    if (r > 1) {
      stop("This estimator requires concentration ratio <= 1.", call. = FALSE)
    }

    eig <- eigen(S, symmetric = TRUE)
    lambda <- eig$values
    U <- eig$vectors

    tolerance <- max(lambda) * max(n, p) * .Machine$double.eps
    if (min(lambda) <= tolerance) {
      stop(paste0(
        "The sample covariance matrix is singular or numerically singular.\n",
        "Check for constant, duplicate, or linearly dependent columns."),
        call. = FALSE)
    }

    h <- min(r^2, 1 / r^2)^0.35 / p^0.35

    theta <- rowMeans(outer(lambda, lambda, function(li, lj) {
      li * (li - lj) / ((li - lj)^2 + (h * li)^2)
    }))

    d <- pmax((1 - r) / lambda + 2 * r / lambda * theta,
              min(1 / lambda))

    result <- tcrossprod(sweep(U, 2L, sqrt(1 / d), `*`))
  }

  dimnames(result) <- dimnames(S)
  return(result)
}
