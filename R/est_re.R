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

  # compute components for A = 0 based on symmetry with A = 1 case
  g_shifted_A1 <- g_out$g_est$g_pred_shifted
  g_shifted_A0 <- 1 - g_shifted_A1 
  e_pred_A1 <- e_out$e_est$e_pred
  e_pred_A0 <- 1 - e_pred_A1

  # get indices of treatment and control from observed data
  idx_A1 <- which(data$A == 1)
  idx_A0 <- which(data$A == 0)

  # subset computed components based on observed treatment status for g
  g_shifted_obs <- rep(NA, nrow(data))
  g_shifted_A1_obs <- g_shifted_A1[idx_A1]
  g_shifted_A0_obs <- g_shifted_A0[idx_A0]
  g_shifted_obs[idx_A1] <- g_shifted_A1_obs
  g_shifted_obs[idx_A0] <- g_shifted_A0_obs

  # subset computed components based on observed treatment status for e
  e_pred_obs <- rep(NA, nrow(data))
  e_pred_A1_obs <- e_out$e_est$e_pred[idx_A1]
  e_pred_A0_obs <- e_out$e_est$e_pred[idx_A0]
  e_pred_obs[idx_A1] <- e_pred_A1_obs
  e_pred_obs[idx_A0] <- e_pred_A0_obs

  # stabilize weights in A-IPW by dividing by sample average since E[g/e] = 1
  mean_aipw <- mean(g_shifted_obs / e_pred_obs)

  # compute estimator
  estim_re <- mean(((g_shifted_obs / e_pred_obs) / mean_aipw) * data$Y)

  # output
  return(estim_re)
}

