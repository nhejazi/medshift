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

  # get indices of treated and control units in validation data
  idx_A1 <- which(data$A == 1)
  idx_A0 <- which(data$A == 0)

  # compute IPW  estimator components from estimates of nuisance parameters
  ipw_out <- compute_ipw(g_output = g_out, e_output = e_out,
                         idx_treat = idx_A1, idx_cntrl = idx_A0)
  g_shifted <- ipw_out$g_shifted
  e_pred <- ipw_out$e_pred
  mean_aipw <- ipw_out$mean_aipw

  # compute estimator
  estim_re <- mean(((g_shifted / e_pred) / mean_aipw) * data$Y)

  # output
  return(estim_re)
}

