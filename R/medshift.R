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
                     estimator = c("efficient", "substitution",
                                   "reweighted")) {
  # set defaults
  estimator <- match.arg(estimator)

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(W)[2]), sep = "_")
  z_names <- paste("Z", seq_len(dim(Z)[2]), sep = "_")
  data.table::setnames(data, c("Y", z_names, "A", w_names))

  # SUBSTITUTION ESTIMATOR TEMPLATE
  if (estimator == "substitution") {
    # fit regression for incremental propensity score intervention
    g_out <- fit_g_mech(data = data, delta_shift = shift_value,
                        lrnr_stack = g_lrnrs)

    # fit regression for incremental propensity score intervention
    m_out <- fit_m_mech(data = data, lrnr_stack = m_lrnrs,
                        z_names = z_names, w_names = w_names)

    # build estimate

  # REWEIGHTED ESTIMATOR TEMPLATE
  } else if (estimator == "reweighted") {
    # fit regression for incremental propensity score intervention
    g_out <- fit_g_mech(data = data, delta_shift = shift_value,
                        lrnr_stack = g_lrnrs)

    # fit clever regression for treatment, conditional on mediators
    e_out <- fit_e_mech(data = data, lrnr_stack = e_lrnrs,
                        z_names = z_names, w_names = w_names)

    # build estimate
    g_shifted <- g_out$g_est$g_pred_shifted
    e_pred <- e_out$e_pred
    estim <- mean((g_shifted / e_pred) * data$Y)

  # EFFICIENT ESTIMATOR TEMPLATE
  } else if (estimator == "efficient") {
    # TODO: use origami to perform CV-SL, fitting each EIF component per fold
  }

  # TODO: output

}
