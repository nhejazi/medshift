#' Fit clever propensity score or intervention density for mediators
#'
#' @param task ...
#' @param lrnr_stack ...
#'
#' @importFrom data.table as.data.table
#' @importFrom sl3 sl3_Task
#
fit_e_mech <- function(task, lrnr_stack) {
  # create clever task for mediation regression
  e_task <- task$get_regression_task("E")

  # fit and predict
  e_fit_stack <- lrnr_stack$train(e_task)
  e_pred <- e_fit_stack$predict()

  # bounding for positivity
  e_pred[e_pred == 0] <- .Machine$double.neg.eps
  e_pred[e_pred == 1] <- 1 - .Machine$double.neg.eps

  # output
  out <- data.table::data.table(cbind(e_pred = e_pred))
  return(out)
}
