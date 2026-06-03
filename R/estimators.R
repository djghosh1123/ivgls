#' IV-LASSO: Two-stage LASSO without graph structure
#'
#' @param Y Numeric vector of length n. Outcome.
#' @param X Numeric n x p matrix of endogenous exposures.
#' @param Z Numeric n x q matrix of instruments.
#' @return Numeric vector of length p of estimated causal effects.
#' @export
iv_lasso <- function(Y, X, Z) {
  n <- nrow(X)
  p <- ncol(X)

  X_hat <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) {
    cv_fit     <- glmnet::cv.glmnet(Z, X[, j], alpha = 1)
    X_hat[, j] <- stats::predict(cv_fit, newx = Z, s = "lambda.min")
  }

  cv_beta <- glmnet::cv.glmnet(X_hat, Y, alpha = 1)
  as.vector(stats::coef(cv_beta, s = "lambda.min"))[-1]
}


#' IVGL: IV regression with graph-fused Lasso
#'
#' @param Y Numeric vector of length n. Outcome.
#' @param X Numeric n x p matrix of endogenous exposures.
#' @param Z Numeric n x q matrix of instruments.
#' @param L Numeric p x p graph Laplacian (see \code{\link{get_laplacian}}).
#' @return Numeric vector of length p of estimated causal effects.
#' @export
ivgl <- function(Y, X, Z, L) {
  n <- nrow(X)
  p <- ncol(X)

  X_hat <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) {
    cv_fit     <- glmnet::cv.glmnet(Z, X[, j], alpha = 1)
    X_hat[, j] <- stats::predict(cv_fit, newx = Z, s = "lambda.min")
  }

  cv_fit <- glmgraph::cv.glmgraph(X_hat, Y, L = L,
                                  family = "gaussian", trace = FALSE)
  as.vector(cv_fit$beta.min)[-1]
}


#' IVGL-S: IV regression with graph Lasso and invalid-IV correction
#'
#' Extends IVGL with an alternating sisVIVE-style update to handle
#' partially invalid instruments that violate the exclusion restriction.
#'
#' @param Y Numeric vector of length n. Outcome.
#' @param X Numeric n x p matrix of endogenous exposures.
#' @param Z Numeric n x q matrix of instruments.
#' @param L Numeric p x p graph Laplacian.
#' @param max_iter Maximum number of alternating iterations. Default 20.
#' @param verbose Print CV loss at each iteration. Default FALSE.
#' @return A list with \code{beta} (length p causal effects) and
#'   \code{alpha} (length q direct IV-outcome effects).
#' @export
ivgl_s <- function(Y, X, Z, L, max_iter = 20, verbose = FALSE) {
  p     <- ncol(X)
  q     <- ncol(Z)
  beta  <- rep(0, p)
  alpha <- rep(0, q)

  for (iter in seq_len(max_iter)) {

    Y_tilde_a <- PZ(Y - X %*% beta, Z)
    cv_a      <- glmnet::cv.glmnet(Z, Y_tilde_a, alpha = 1,
                                   intercept = FALSE, standardize = TRUE)
    alpha     <- as.vector(stats::coef(cv_a, s = "lambda.1se"))[-1]

    Y_tilde_b <- PZ(Y - Z %*% alpha, Z)
    cv_b      <- glmgraph::cv.glmgraph(X, Y_tilde_b, L = L,
                                       family = "gaussian", trace = FALSE)
    beta      <- as.vector(cv_b$beta.min)[-1]

    if (verbose) {
      message(sprintf("Iter %d | CV loss = %.6f", iter, cv_b$cvmin))
    }
  }

  list(beta = beta, alpha = alpha)
}
