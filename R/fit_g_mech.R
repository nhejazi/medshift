#' Fit propensity score with incremental shift
#'
#' @param data ...
#' @param delta_shift ...
#' @param lrnr_stack ...
#'
#' @importFrom data.table as.data.table
#
fit_g_mech <- function(data, delta_shift, lrnr_stack) {
  #  construct task for propensity score fit
  g_task <- sl3_Task$new(data = data,
                         covariates = w_names,
                         outcome = "A")
  # fit and predict
  g_fit_stack <- lrnr_stack$train(g_task)

  # copy data and set intervention A = 1
  data_A1 <- copy(data)
  data_A1[, A := rep(1, .N)]
  g_task_A1 <- sl3_Task$new(data = data_A1,
                            covariates = w_names,
                            outcome = "A")
  g_pred_A1 <- g_fit_stack$predict(g_task_A1)
  g_pred_A0 <- 1 - g_pred_A1

  # directly computed the shifted propensity score
  g_pred_shifted <- (delta_shift * g_pred_A1) /
    (delta_shift * g_pred_A1 + (1 - g_pred_A1))

  # output
  out <- list(
    g_est = data.table::data.table(cbind(
    g_pred_A1 = g_pred_A1,
    g_pred_A0 = g_pred_A0,
    g_pred_shifted = g_pred_shifted
    )),
    g_fit = g_fit_stack)
  return(out)
}
