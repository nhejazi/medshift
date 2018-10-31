#' Substitution estimator
#'
#' @param task ...
#' @param delta_shift ...
#' @param lrnr_stack_g ...
#' @param lrnr_stack_m ...
#
estim_sub <- function(task, delta_shift, lrnr_stack_g, lrnr_stack_m) {
  # fit propensity score
  g_est <- fit_g_mech(
    task = task, delta_shift = delta_shift,
    lrnr_stack = lrnr_stack_g
  )
  g_pred_shifted <- g_est$g_shifted

  # compute estimator for the outcome regression
  m_out <- fit_m_mech(task = task, lrnr_stack = lrnr_stack_m)
  m_pred <- m_out$m_pred

  # loop over levels of intervention
  # TODO: FIX THIS, VERY BAD
  theta_sub <- 0
  for (a in unique(task$data$A)) {
    m_this_a <- m_pred[task$data$A == a]
    g_this_a <- g_pred_shifted[task$data$A == a]
    theta_sub_this_a <- mean(m_this_a * g_this_a)
    theta_sub <- theta_sub + theta_sub_this_a
  }

  # return value of substitution estimator
  return(theta_sub)
}
