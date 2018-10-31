#' Fit outcome regression
#'
#' @param task ...
#' @param lrnr_stack ...
#'
#' @importFrom data.table as.data.table
#
fit_m_mech <- function(task, lrnr_stack) {

  # make task for propensity score regression
  m_task <- task$get_regression_task("Y")

  # fit and predict
  m_fit_stack <- lrnr_stack$train(m_task)
  m_pred <- m_fit_stack$predict()

  # output
  out <- list(m_pred = data.table::data.table(cbind(m_pred = m_pred)),
              m_fit_sl = m_fit_stack)
  return(out)
}
