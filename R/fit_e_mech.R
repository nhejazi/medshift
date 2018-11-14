#' Fit clever propensity score regression for mediators
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
fit_e_mech <- function(data, valid_data = NULL,
                       lrnr_stack, z_names, w_names) {
  # construct task for nuisance parameter fit
  e_task <- sl3::sl3_Task$new(data = data,
                              covariates = c(z_names, w_names),
                              outcome = "A")

  # fit and predict
  e_fit_stack <- lrnr_stack$train(e_task)

  if (is.null(valid_data)) {
    # predict on full data
    e_pred <- e_fit_stack$predict()
  } else {
    # construct task for nuisance parameter fit for validation data
    e_task_valid <- sl3::sl3_Task$new(data = valid_data,
                                      covariates = c(z_names, w_names),
                                      outcome = "A")
    # predict on validation data
    e_pred <- e_fit_stack$predict(e_task_valid)
  }

  # bounding for positivity
  e_pred[e_pred == 0] <- .Machine$double.neg.eps
  e_pred[e_pred == 1] <- 1 - .Machine$double.neg.eps

  # output
  out <- data.table::data.table(cbind(e_pred = e_pred))
  return(out)
}
