#' Construct a graph adjacency matrix
#'
#' @param p Number of nodes.
#' @param type One of \code{"proximity"}, \code{"ring"}, \code{"chain"},
#'   \code{"community"}, or \code{"disconnected"}.
#' @return A symmetric p x p binary adjacency matrix.
#' @export
make_graph <- function(p = 70, type = c("proximity", "ring", "chain",
                                        "community", "disconnected")) {
  type <- match.arg(type)

  A <- switch(type,

              proximity = {
                coords <- matrix(stats::runif(p * 3, 0, 100), ncol = 3)
                Dmat   <- as.matrix(stats::dist(coords))
                1 * ((Dmat < 30) & (Dmat > 0))
              },

              ring = {
                M <- matrix(0L, p, p)
                for (j in seq_len(p)) {
                  nxt <- ifelse(j == p, 1L, j + 1L)
                  M[j, nxt] <- M[nxt, j] <- 1L
                }
                M
              },

              chain = {
                M <- matrix(0L, p, p)
                for (j in seq_len(p - 1L)) {
                  M[j, j + 1L] <- M[j + 1L, j] <- 1L
                }
                M
              },

              community = {
                blk <- c(floor(p / 3), floor(p / 3), p - 2L * floor(p / 3))
                pm  <- matrix(c(0.35, 0.03, 0.03,
                                0.03, 0.35, 0.03,
                                0.03, 0.03, 0.35), 3, 3)
                g   <- igraph::sample_sbm(n = p, pref.matrix = pm,
                                          block.sizes = blk,
                                          directed = FALSE, loops = FALSE)
                as.matrix(igraph::as_adjacency_matrix(g))
              },

              disconnected = {
                M   <- matrix(0L, p, p)
                blk <- c(floor(p / 3), floor(p / 3), p - 2L * floor(p / 3))
                start <- 1L
                for (b in blk) {
                  idx <- start:(start + b - 1L)
                  for (j in idx[-length(idx)]) {
                    M[j, j + 1L] <- M[j + 1L, j] <- 1L
                  }
                  start <- start + b
                }
                M
              }
  )

  diag(A) <- 0
  A
}


#' Corrupt a graph by random edge swaps
#'
#' @param A Symmetric p x p binary adjacency matrix.
#' @param corruption_rate Proportion of edges to remove and replace.
#' @return A corrupted adjacency matrix.
#' @export
corrupt_graph <- function(A, corruption_rate = 0.30) {
  A_new      <- A
  exist_idx  <- which(A == 1 & upper.tri(A), arr.ind = TRUE)
  noedge_idx <- which(A == 0 & upper.tri(A), arr.ind = TRUE)
  n_swap     <- floor(corruption_rate * nrow(exist_idx))
  if (n_swap == 0) return(A_new)

  rem <- exist_idx[sample.int(nrow(exist_idx), n_swap), , drop = FALSE]
  add <- noedge_idx[sample.int(nrow(noedge_idx), n_swap), , drop = FALSE]

  for (k in seq_len(n_swap)) {
    A_new[rem[k, 1], rem[k, 2]] <- A_new[rem[k, 2], rem[k, 1]] <- 0
    A_new[add[k, 1], add[k, 2]] <- A_new[add[k, 2], add[k, 1]] <- 1
  }
  A_new
}
