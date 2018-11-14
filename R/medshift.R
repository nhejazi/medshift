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
#' @importFrom origami make_folds cross_validate
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
                        lrnr_stack = g_lrnrs, w_names = w_names)

    # fit regression for incremental propensity score intervention
    m_out <- fit_m_mech(data = data, lrnr_stack = m_lrnrs,
                        z_names = z_names, w_names = w_names)

    # build estimate
    g_shifted_A1 <- g_out$g_est$g_pred_shifted
    g_shifted_A0 <- 1 - g_shifted_A1
    m_pred_A1 <- m_out$m_pred$m_pred_A1
    m_pred_A0 <- m_out$m_pred$m_pred_A0
    estim_sub <- mean(m_pred_A0 * g_shifted_A0) +
      mean(m_pred_A1 * g_shifted_A1)

    # output
    est_out <- estim_sub

  # REWEIGHTED ESTIMATOR TEMPLATE
  } else if (estimator == "reweighted") {
    # fit regression for incremental propensity score intervention
    g_out <- fit_g_mech(data = data, delta_shift = shift_value,
                        lrnr_stack = g_lrnrs, w_names = w_names)

    # fit clever regression for treatment, conditional on mediators
    e_out <- fit_e_mech(data = data, lrnr_stack = e_lrnrs,
                        z_names = z_names, w_names = w_names)

    # build estimate
    g_shifted <- g_out$g_est$g_pred_shifted
    e_pred <- e_out$e_pred
    estim_re <- mean((g_shifted / e_pred) * data$Y)

    # output
    est_out <- estim_re

  # EFFICIENT ESTIMATOR TEMPLATE
  } else if (estimator == "efficient") {
    # use origami to perform CV-SL, fitting each EIF component per fold
    # and evaluating only on the holdouts
    eif_component_names <- c("Dy", "Da", "Dzw")
    folds <- origami::make_folds(data)
    cv_eif_results <- origami::cross_validate(cv_fun = cv_eif,
                                              folds = folds,
                                              data = data,
                                              delta_shift = shift_value,
                                              lrnr_stack_g = g_lrnrs,
                                              lrnr_stack_e = e_lrnrs,
                                              lrnr_stack_m = m_lrnrs,
                                              use_future = FALSE,
                                              .combine = FALSE)
    D_obs <- lapply(cv_eif_results[[1]], function(x) {
                      D_obs_fold <- rowSums(x[, ..eif_component_names])
                    })
    estim_eff <- mean(do.call(c, D_obs))

    # output
    est_out <- estim_eff
  }
  # final output
  return(est_out)
}
