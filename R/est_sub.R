#' Substitution estimator
#'
#' @param data ...
#' @param delta ...
#' @param g_lrnrs ...
#' @param m_lrnrs ...
#' @param w_names ...
#' @param z_names ...
#'
#' @keywords internal
#
est_sub <- function(data,
                    delta,
                    g_lrnrs,
                    m_lrnrs,
                    w_names,
                    z_names) {
  # estimate propensity score
  g_out <- fit_g_mech(
    data = data, delta = delta,
    lrnr_stack = g_lrnrs, w_names = w_names
  )

  # fit regression for incremental propensity score intervention
  m_out <- fit_m_mech(
    data = data, lrnr_stack = m_lrnrs,
    z_names = z_names, w_names = w_names
  )

  # build estimate
  g_shifted_A1 <- g_out$g_est$g_pred_shifted
  g_shifted_A0 <- 1 - g_shifted_A1
  m_pred_A1 <- m_out$m_pred$m_pred_A1
  m_pred_A0 <- m_out$m_pred$m_pred_A0

  # find indices of observed treatment and control
  idx_A1 <- which(data$A == 1)
  idx_A0 <- which(data$A == 0)

  # subset predictions to evaluation under observed data
  g_shifted_A1_obs <- g_shifted_A1[idx_A1]
  g_shifted_A0_obs <- g_shifted_A0[idx_A0]
  m_pred_A1_obs <- m_pred_A1[idx_A1]
  m_pred_A0_obs <- m_pred_A0[idx_A0]

  # compute estimator
  estim_sub <- mean(m_pred_A0_obs * g_shifted_A0_obs) +
    mean(m_pred_A1_obs * g_shifted_A1_obs)

  # output
  return(estim_sub)
}

