#' @keywords internal
grow_cluster <- function(A, seed, max_size) {
  cluster  <- seed
  frontier <- which(A[seed, ] == 1)
  while (length(cluster) < max_size && length(frontier) > 0) {
    nxt      <- frontier[1]
    cluster  <- unique(c(cluster, nxt))
    frontier <- unique(setdiff(c(frontier[-1], which(A[nxt, ] == 1)), cluster))
  }
  cluster
}


#' Generate a sparse true coefficient vector on a graph
#'
#' @param A Symmetric p x p adjacency matrix.
#' @param s2 Number of active nodes.
#' @param signal Causal effect magnitude.
#' @param pattern One of \code{"smooth"}, \code{"nonsmooth"}, or \code{"community"}.
#' @param smooth_noise SD of noise added around the base signal.
#' @return A list with \code{beta_true} and \code{active_set}.
#' @export
generate_beta <- function(A, s2 = 5, signal = 3,
                          pattern = c("smooth", "nonsmooth", "community"),
                          smooth_noise = 0.20) {
  pattern   <- match.arg(pattern)
  p         <- nrow(A)
  beta_true <- rep(0, p)

  if (pattern == "smooth") {
    seed       <- sample.int(p, 1)
    active_set <- grow_cluster(A, seed, s2)
    if (length(active_set) < s2) {
      extra      <- sample(setdiff(seq_len(p), active_set), s2 - length(active_set))
      active_set <- c(active_set, extra)
    }
    active_set <- active_set[seq_len(s2)]
    base_effect <- sample(c(-1, 1), 1) * signal
    beta_true[active_set] <- base_effect + smooth_noise * stats::rnorm(s2)

  } else if (pattern == "nonsmooth") {
    active_set <- sample.int(p, s2)
    beta_true[active_set] <- stats::runif(s2, 0.5, 1) *
      sample(c(-1, 1), s2, TRUE) * signal

  } else {
    comps      <- igraph::components(
      igraph::graph_from_adjacency_matrix(A, mode = "undirected"))
    biggest    <- which.max(comps$csize)
    candidates <- which(comps$membership == biggest)
    active_set <- sample(candidates, min(s2, length(candidates)))
    beta_true[active_set] <- signal + smooth_noise * stats::rnorm(length(active_set))
  }

  list(beta_true = beta_true, active_set = active_set)
}


#' Simulate data for graph-IV regression
#'
#' @param n Sample size.
#' @param p Number of exposures.
#' @param q Number of instruments.
#' @param s1 Fraction of instruments relevant for each exposure.
#' @param s_alpha Number of invalid instruments.
#' @param alpha_strength Direct-effect magnitude of invalid instruments.
#' @param beta_true Numeric vector of length p of true causal effects.
#' @return A list with \code{Y}, \code{X}, \code{Z}, \code{A_true},
#'   and \code{alpha_true}.
#' @export
generate_data <- function(n = 100, p = 70, q = 500,
                          s1 = 0.10,
                          s_alpha = 10,
                          alpha_strength = 5,
                          beta_true) {
  Z      <- matrix(stats::rnorm(n * q), n, q)
  A_true <- matrix(0, q, p)
  for (j in seq_len(p)) {
    iv_idx <- sample.int(q, ceiling(s1 * q))
    A_true[iv_idx, j] <- stats::rnorm(length(iv_idx))
  }
  X          <- Z %*% A_true + matrix(stats::rnorm(n * p), n, p)
  alpha_true <- rep(0, q)
  alpha_true[seq_len(s_alpha)] <- alpha_strength
  Y          <- as.vector(X %*% beta_true + Z %*% alpha_true + stats::rnorm(n))

  list(Y = Y, X = X, Z = Z, A_true = A_true, alpha_true = alpha_true)
}
