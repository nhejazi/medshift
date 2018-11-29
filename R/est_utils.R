#' Compute Dzw component of efficient influence funtion
#'
#' @param g_output ...
#' @param m_output ...
#
make_Dzw <- function(g_output, m_output) {
  # get m components from output for that nuisance parameter
  m_pred_A1 <- m_output$m_pred$m_pred_A1
  m_pred_A0 <- m_output$m_pred$m_pred_A0

  # compute component Dzw from nuisance parameters
  g_shifted_A1 <- g_output$g_est$g_pred_shifted_A1
  g_shifted_A0 <- g_output$g_est$g_pred_shifted_A0
  Dzw <- (m_pred_A0 * g_shifted_A0) + (m_pred_A1 * g_shifted_A1)
  return(Dzw)
}

################################################################################

