#' Nonparametric estimation of direct effects under mediation
#'
#' @param W ...
#' @param A ...
#' @param Z ...
#' @param Y ...
#' @param shift_value ...
#' @param g_lrnrs ...
#' @param e_lrnrs ...
#' @param Q_lrnrs ...
#' @param estimator ...
#'
#' @importFrom data.table as.data.table setnames
#
medshift <- function(W,
                     A,
                     Z,
                     Y,
                     shift_value = 0.5,
                     g_lrnrs = Lrnr_glm_fast$new(family = binomial()),
                     e_lrnrs = Lrnr_glm_fast$new(family = binomial()),
                     m_lrnrs = Lrnr_glm_fast$new(),
                     estimator = c("est_eqn", "sub", "reweighted")) {

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(W)[2]), sep = "_")
  z_names <- paste("Z", seq_len(dim(Z)[2]), sep = "_")
  data.table::setnames(data, c("Y", z_names, "A", w_names))

  # fit regression for incremental propensity score intervention
  g_out <- fit_g_mech(data, delta_shift = shift_value, lrnr_stack = g_lrnrs)


  # TODO: functions for estimating likelihood components
}
