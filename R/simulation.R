#' Compute performance metrics for support recovery
#'
#' @param true_support Integer vector of true active indices.
#' @param estimated_support Integer vector of estimated active indices.
#' @param p Total number of predictors.
#' @return Named numeric vector with MCC, TPR, FPR, and Selected.
#' @export
eval_support <- function(true_support, estimated_support, p) {
  sel      <- integer(p); sel[estimated_support] <- 1L
  true_vec <- integer(p); true_vec[true_support]  <- 1L

  TP <- sum(sel == 1L & true_vec == 1L)
  TN <- sum(sel == 0L & true_vec == 0L)
  FP <- sum(sel == 1L & true_vec == 0L)
  FN <- sum(sel == 0L & true_vec == 1L)

  denom_mcc <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (denom_mcc == 0) 0 else (TP * TN - FP * FN) / denom_mcc
  tpr <- if ((TP + FN) == 0) 0 else TP / (TP + FN)
  fpr <- if ((FP + TN) == 0) 0 else FP / (FP + TN)

  c(MCC = mcc, TPR = tpr, FPR = fpr, Selected = sum(sel))
}


#' Run a single simulation replicate
#'
#' @param n Sample size.
#' @param p Number of exposures.
#' @param q Number of instruments.
#' @param graph_type Graph topology passed to \code{\link{make_graph}}.
#' @param signal_pattern One of \code{"smooth"}, \code{"nonsmooth"},
#'   or \code{"community"}.
#' @param fit_graph_type Graph supplied to estimators. If NULL uses
#'   \code{graph_type}.
#' @param graph_corruption Proportion of edges to corrupt. Default 0.
#' @param s2 Number of active nodes.
#' @param signal Causal effect magnitude.
#' @param s_alpha Number of invalid instruments.
#' @param alpha_strength Invalid-IV direct-effect magnitude.
#' @param smooth_noise Noise on the smooth signal pattern.
#' @param threshold Coefficients below this are treated as zero.
#' @return A data.frame with one row per method and columns Method,
#'   MSE, MCC, TPR, FPR, Selected.
#' @export
run_one_replicate <- function(n = 100, p = 70, q = 500,
                              graph_type       = "proximity",
                              signal_pattern   = "smooth",
                              fit_graph_type   = NULL,
                              graph_corruption = 0,
                              s2             = 5,
                              signal         = 3,
                              s_alpha        = 10,
                              alpha_strength = 5,
                              smooth_noise   = 0.20,
                              threshold      = 1e-4) {

  A_true_graph <- make_graph(p, graph_type)
  beta_obj     <- generate_beta(A_true_graph, s2 = s2, signal = signal,
                                pattern = signal_pattern,
                                smooth_noise = smooth_noise)
  beta_true  <- beta_obj$beta_true
  active_set <- beta_obj$active_set

  dat <- generate_data(n = n, p = p, q = q,
                       s_alpha = s_alpha,
                       alpha_strength = alpha_strength,
                       beta_true = beta_true)

  A_fit <- if (is.null(fit_graph_type)) A_true_graph else make_graph(p, fit_graph_type)
  if (graph_corruption > 0) A_fit <- corrupt_graph(A_fit, graph_corruption)
  L_fit <- get_laplacian(A_fit)

  beta_ivl  <- iv_lasso(dat$Y, dat$X, dat$Z)
  beta_ivgl <- ivgl(dat$Y, dat$X, dat$Z, L_fit)
  fit_s     <- ivgl_s(dat$Y, dat$X, dat$Z, L_fit)

  make_row <- function(method, beta_est) {
    supp <- which(abs(beta_est) > threshold)
    mets <- eval_support(active_set, supp, p)
    data.frame(Method = method,
               MSE    = mean((beta_est - beta_true)^2),
               MCC    = mets["MCC"],
               TPR    = mets["TPR"],
               FPR    = mets["FPR"],
               Selected = mets["Selected"],
               row.names = NULL)
  }

  rbind(
    make_row("IV_LASSO", beta_ivl),
    make_row("IVGL",     beta_ivgl),
    make_row("IVGL_S",   fit_s$beta)
  )
}


#' Run a simulation study with multiple replicates
#'
#' @param B Number of Monte Carlo replicates.
#' @inheritParams run_one_replicate
#' @return A data.frame with 3*B rows, one per method per replicate.
#' @export
run_simulation <- function(B = 100, n = 100, p = 70, q = 500,
                           graph_type       = "proximity",
                           signal_pattern   = "smooth",
                           fit_graph_type   = NULL,
                           graph_corruption = 0,
                           s2             = 5,
                           signal         = 3,
                           s_alpha        = 10,
                           alpha_strength = 5,
                           smooth_noise   = 0.20,
                           threshold      = 1e-4) {
  results <- vector("list", B)
  for (b in seq_len(B)) {
    results[[b]] <- run_one_replicate(
      n = n, p = p, q = q,
      graph_type       = graph_type,
      signal_pattern   = signal_pattern,
      fit_graph_type   = fit_graph_type,
      graph_corruption = graph_corruption,
      s2             = s2,
      signal         = signal,
      s_alpha        = s_alpha,
      alpha_strength = alpha_strength,
      smooth_noise   = smooth_noise,
      threshold      = threshold
    )
  }
  do.call(rbind, results)
}
