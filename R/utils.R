#' Project a vector onto the column space of Z
#' @keywords internal
PZ <- function(v, Z) {
  Z %*% MASS::ginv(crossprod(Z)) %*% t(Z) %*% v
}


#' Matthews Correlation Coefficient for support recovery
#'
#' @param true_support Integer vector of truly active indices.
#' @param estimated_support Integer vector of estimated active indices.
#' @param p Total number of predictors.
#' @return A scalar between -1 and 1.
#' @export
get_mcc <- function(true_support, estimated_support, p) {
  sel      <- integer(p); sel[estimated_support] <- 1L
  true_vec <- integer(p); true_vec[true_support]  <- 1L

  TP <- sum(sel == 1L & true_vec == 1L)
  TN <- sum(sel == 0L & true_vec == 0L)
  FP <- sum(sel == 1L & true_vec == 0L)
  FN <- sum(sel == 0L & true_vec == 1L)

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (denom == 0) return(0)
  (TP * TN - FP * FN) / denom
}


#' Compute the unnormalised graph Laplacian
#'
#' @param A Symmetric p x p binary adjacency matrix.
#' @return A p x p Laplacian matrix.
#' @export
get_laplacian <- function(A) {
  diag(rowSums(A)) - A
}
