#' Re-weighted estimator
#'
#' @param task ...
#' @param delta_shift ...
#' @param lrnr_stack ...
#
estim_reweighted <- function(task, delta_shift, lrnr_stack) {
  # fit propensity score
  g_est <- fit_g_mech(task = task, delta_shift = delta_shift,
                      lrnr_stack = lrnr_stack)
  g_pred_shifted <- g_est$g_shifted

  # compute the clever intervention density for the mediator
  e_est <- fit_e_mech(task = task, lrnr_stack = lrnr_stack)
  e_pred <- e_est$e_pred

  # build the re-weighted estimator and return as output
  theta_re <- mean((g_pred_shifted / e_pred) * task$data$Y)
  return(theta_re)
}

