#' Fit outcome regression
#'
#' @param data ...
#' @param valid_data A validation/holdout data set (optional).
#' @param lrnr_stack ...
#' @param z_names ...
#' @param w_names ...
#'
#' @importFrom data.table as.data.table
#' @importFrom sl3 sl3_Task
#
fit_m_mech <- function(data, valid_data = NULL,
                       lrnr_stack, z_names, w_names) {
  #  construct task for propensity score fit
  m_task <- sl3::sl3_Task$new(data = data,
                              covariates = c("A", z_names, w_names),
                              outcome = "Y")

  # fit and predict
  m_fit_stack <- lrnr_stack$train(m_task)

  if (is.null(valid_data)) {
    # copy full data
    data_A1 <- copy(data)
    data_A0 <- copy(data)
  } else {
    # copy only validation data
    data_A1 <- copy(valid_data)
    data_A0 <- copy(valid_data)
  }

  # copy data and set intervention A = 1
  data_A1[, A := rep(1, .N)]
  m_task_A1 <- sl3::sl3_Task$new(data = data_A1,
                                 covariates = c("A", z_names, w_names),
                                 outcome = "Y")
  m_pred_A1 <- m_fit_stack$predict(m_task_A1)

  # copy data and set intervention A = 0
  data_A0[, A := rep(0, .N)]
  m_task_A0 <- sl3::sl3_Task$new(data = data_A0,
                                 covariates = c("A", z_names, w_names),
                                 outcome = "Y")
  m_pred_A0 <- m_fit_stack$predict(m_task_A0)

  # output
  out <- list(m_pred = data.table::data.table(cbind(m_pred_A1 = m_pred_A1,
                                                    m_pred_A0 = m_pred_A0)),
              m_fit_sl = m_fit_stack)
  return(out)
}
