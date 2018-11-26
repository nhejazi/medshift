utils::globalVariables(c("..eif_component_names"))

#' Efficient estimator
#'
#' @param data ...
#' @param delta ...
#' @param g_lrnrs ...
#' @param e_lrnrs ...
#' @param m_lrnrs ...
#' @param phi_lrnrs ...
#' @param w_names ...
#' @param z_names ...
#'
#' @importFrom origami make_folds cross_validate
#'
#' @keywords internal
#
est_eff <- function(data,
                    delta,
                    g_lrnrs,
                    e_lrnrs,
                    m_lrnrs,
                    phi_lrnrs,
                    w_names,
                    z_names) {
  # use origami to perform CV-SL, fitting/evaluating EIF components per fold
  eif_component_names <- c("Dy", "Da", "Dzw")

  # create folds for use with origami::cross_validate
  folds <- origami::make_folds(data)

  # perform the cv_eif procedure on a per-fold basis
  cv_eif_results <- origami::cross_validate(
    cv_fun = cv_eif,
    folds = folds,
    data = data,
    delta = delta,
    lrnr_stack_g = g_lrnrs,
    lrnr_stack_e = e_lrnrs,
    lrnr_stack_m = m_lrnrs,
    lrnr_stack_phi = phi_lrnrs,
    z_names = z_names,
    w_names = w_names,
    use_future = FALSE,
    .combine = FALSE
  )

  # combine results of EIF components for full EIF
  D_obs <- lapply(cv_eif_results[[1]], function(x) {
    D_obs_fold <- rowSums(x[, ..eif_component_names])
  })

  # compute parameter estimate as mean of EIF
  estim_eff <- mean(do.call(c, D_obs))

  # output
  return(estim_eff)
}
