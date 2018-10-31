#' Re-weighted estimator
#'
#' @param task ...
#' @param delta_shift ...
#' @param lrnr_stack_g ...
#' @param lrnr_stack_e ...
#
estim_re <- function(task, delta_shift, lrnr_stack_g, lrnr_stack_e) {
  # fit propensity score
  g_est <- fit_g_mech(
    task = task, delta_shift = delta_shift,
    lrnr_stack = lrnr_stack_g
  )
  g_pred_shifted <- g_est$g_shifted

  # compute the clever intervention density for the mediator
  e_est <- fit_e_mech(task = task, lrnr_stack = lrnr_stack_e)
  e_pred <- e_est$e_pred

  # build the re-weighted estimator and return as output
  theta_re <- mean((g_pred_shifted / e_pred) * task$data$Y)
  return(theta_re)
}
