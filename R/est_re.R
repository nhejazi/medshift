#' Re-weighted estimator
#'
#' @param data ...
#' @param delta ...
#' @param g_lrnrs ...
#' @param e_lrnrs ...
#' @param w_names ...
#' @param z_names ...
#'
#' @keywords internal
#
est_re <- function(data,
                   delta,
                   g_lrnrs,
                   e_lrnrs,
                   w_names,
                   z_names) {
  # fit regression for incremental propensity score intervention
  g_out <- fit_g_mech(
    data = data, delta = delta,
    lrnr_stack = g_lrnrs, w_names = w_names
  )

  # fit clever regression for treatment, conditional on mediators
  e_out <- fit_e_mech(
    data = data, lrnr_stack = e_lrnrs,
    z_names = z_names, w_names = w_names
  )

  # stabilize weights in AIPW by dividing by sample average since E[g/e] = 1
  g_shifted <- g_out$g_est$g_pred_shifted
  e_pred <- e_out$e_est$e_pred
  mean_aipw <- mean(g_shifted / e_pred)

  # build estimate
  estim_re <- mean(((g_shifted / e_pred) / mean_aipw) * data$Y)

  # output
  return(estim_re)
}
