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

  # compute Dzw component of EIF using convenience function
  Dzw_groupwise <- make_Dzw(g_output = g_out, m_output = m_out)

  # compute estimator
  estim_sub <- mean(Dzw_groupwise$dzw_cntrl) + mean(Dzw_groupwise$dzw_treat)

  # output
  return(estim_sub)
}

