#' Fit propensity score or intervention density
#'
#' @param task ...
#' @param delta_shift ...
#' @param lrnr_stack ...
#'
#' @importFrom data.table as.data.table
#
fit_g_mech <- function(task, delta_shift, lrnr_stack) {

  # make task for propensity score regression
  g_task <- task$get_regression_task("A")

  # fit and predict
  g_fit_stack <- lrnr_stack$train(g_task)
  g_pred_natural <- g_fit_stack$predict()

  # directly computed the shifted propensity score
  g_pred_shifted <- (delta_shift * g_pred_natural) /
      (delta_shift * g_pred_natural + (1 - g_pred_natural))

  # output
  out <- data.table::data.table(cbind(g_natural = g_pred_natural,
                                      g_shifted = g_pred_shifted))
  return(out)
}

