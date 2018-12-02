#' Get Dzw component of efficient influence function from nuisance parameters
#'
#' @param g_output ...
#' @param m_output ...
#
compute_Dzw <- function(g_output, m_output) {
  # get g components from output for that nuisance parameter
  g_shifted_A1 <- g_output$g_est$g_pred_shifted_A1
  g_shifted_A0 <- g_output$g_est$g_pred_shifted_A0

  # get m components from output for that nuisance parameter
  m_pred_A1 <- m_output$m_pred$m_pred_A1
  m_pred_A0 <- m_output$m_pred$m_pred_A0

  # compute component Dzw from nuisance parameters
  Dzw_A1 <- g_shifted_A1 * m_pred_A1
  Dzw_A0 <- g_shifted_A0 * m_pred_A0

  # output as simple list
  return(list(dzw_cntrl = Dzw_A0,
              dzw_treat = Dzw_A1))
}

################################################################################

#' Get inverse probability weighted (IPW) estimate from nuisance parameters
#'
#' @param g_output ...
#' @param e_output ...
#' @param idx_treat ...
#' @param idx_cntrl ...
#
compute_ipw <- function(g_output, e_output, idx_treat, idx_cntrl) {
  # compute components for A = 0 based on symmetry with A = 1 case
  g_shifted_A1 <- g_output$g_est$g_pred_shifted_A1
  g_shifted_A0 <- g_output$g_est$g_pred_shifted_A0
  e_pred_A1 <- e_output$e_est$e_pred_A1
  e_pred_A0 <- e_output$e_est$e_pred_A0

  # subset computed components based on observed treatment status for g
  g_shifted_obs <- rep(NA, length(idx_treat) + length(idx_cntrl))
  g_shifted_A1_obs <- g_shifted_A1[idx_treat]
  g_shifted_A0_obs <- g_shifted_A0[idx_cntrl]
  g_shifted_obs[idx_treat] <- g_shifted_A1_obs
  g_shifted_obs[idx_cntrl] <- g_shifted_A0_obs

  # subset computed components based on observed treatment status for e
  e_pred_obs <- rep(NA, length(idx_treat) + length(idx_cntrl))
  e_pred_A1_obs <- e_pred_A1[idx_treat]
  e_pred_A0_obs <- e_pred_A0[idx_cntrl]
  e_pred_obs[idx_treat] <- e_pred_A1_obs
  e_pred_obs[idx_cntrl] <- e_pred_A0_obs

  # stabilize weights in A-IPW by dividing by sample average since E[g/e] = 1
  mean_aipw <- mean(g_shifted_obs / e_pred_obs)

  # output as simple list
  return(list(g_shifted = g_shifted_obs,
              e_pred = e_pred_obs,
              mean_aipw = mean_aipw))
}

