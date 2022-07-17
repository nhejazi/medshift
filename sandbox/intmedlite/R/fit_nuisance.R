#' Fit Nuisance Parameters
#'
#' @import data.table origami sl3
#'
#' @param data ...
#' @param g_stack ...
#' @param e_stack ...
#' @param m_stack ...
#' @param b_stack ...
#' @param d_stack ...
#' @param cv_folds ...
fit_nuisance <- function(data,
                         g_stack,
                         e_stack,
                         m_stack,
                         b_stack,
                         d_stack,
                         cv_folds = 5L) {
  ## make fold structure
  folds <- origami::make_folds(data,
    fold_fun = origami::folds_vfold,
    V = cv_folds
  )

  ## fit nuisance functions with V-fold CV
  cv_fit <- origami::cross_validate(
    cv_fun = cv_fit_nuisance,
    folds = folds,
    data = data,
    g_stack = g_stack,
    e_stack = e_stack,
    m_stack = m_stack,
    b_stack = b_stack,
    d_stack = d_stack,
    use_future = FALSE,
    .combine = FALSE
  )
  cv_preds <- data.table::rbindlist(cv_fit[["preds"]])
  cv_fit[["preds"]] <- NULL

  ## output
  out <- list(cv_fits = cv_fit, cv_preds = cv_preds, folds = folds)
  return(out)
}


#' Fit Nuisance Parameters in a Cross-Validation Fold
#'
#' @import data.table origami sl3 stringr
#'
#' @param fold ...
#' @param data ...
#' @param g_stack ...
#' @param e_stack ...
#' @param m_stack ...
#' @param b_stack ...
#' @param d_stack ...
cv_fit_nuisance <- function(fold,
                            data,
                            g_stack,
                            e_stack,
                            m_stack,
                            b_stack,
                            d_stack) {
  ## make folds
  train_data <- origami::training(data)
  valid_data <- origami::validation(data)

  ## extract data
  Y <- data[, "Y"]
  Z <- data[, stringr::str_detect(names(data), "Z")]
  L <- data[, "L"]
  A <- data[, "A"]
  W <- data[, stringr::str_detect(names(data), "W")]

  ## get names of data components
  names_z <- stringr::str_subset(colnames(data), "Z")
  names_w <- stringr::str_subset(colnames(data), "W")

  ## fit nuisance functions
  g_task_train <- sl3::sl3_Task$new(
    data = train_data,
    covariates = names_w,
    outcome = "A",
    outcome_type = "binomial"
  )
  g_task_valid <- sl3::sl3_Task$new(
    data = valid_data,
    covariates = names_w,
    outcome = "A",
    outcome_type = "binomial"
  )
  g_fit <- g_stack$train(g_task_train)
  g_preds <- g_fit$predict(g_task_valid)

  e_task_train <- sl3::sl3_Task$new(
    data = train_data,
    covariates = c(names_w, names_z),
    outcome = "A",
    outcome_type = "binomial"
  )
  e_task_valid <- sl3::sl3_Task$new(
    data = valid_data,
    covariates = c(names_w, names_z),
    outcome = "A",
    outcome_type = "binomial"
  )
  e_fit <- e_stack$train(e_task_train)
  e_preds <- e_fit$predict(e_task_valid)

  m_task_train <- sl3::sl3_Task$new(
    data = train_data,
    covariates = c(names_w, "A", "L", names_z),
    outcome = "Y",
    outcome_type = "binomial"
  )
  m_task_valid <- sl3::sl3_Task$new(
    data = valid_data,
    covariates = c(names_w, "A", "L", names_z),
    outcome = "Y",
    outcome_type = "binomial"
  )
  m_fit <- m_stack$train(m_task_train)
  m_preds <- m_fit$predict(m_task_valid)

  b_task_train <- sl3::sl3_Task$new(
    data = train_data,
    covariates = c(names_w, "A"),
    outcome = "L",
    outcome_type = "binomial"
  )
  b_task_valid <- sl3::sl3_Task$new(
    data = valid_data,
    covariates = c(names_w, "A"),
    outcome = "L",
    outcome_type = "binomial"
  )
  b_fit <- b_stack$train(b_task_train)
  b_preds <- b_fit$predict(b_task_valid)

  d_task_train <- sl3::sl3_Task$new(
    data = train_data,
    covariates = c(names_w, "A", names_z),
    outcome = "L",
    outcome_type = "binomial"
  )
  d_task_valid <- sl3::sl3_Task$new(
    data = valid_data,
    covariates = c(names_w, "A", names_z),
    outcome = "L",
    outcome_type = "binomial"
  )
  d_fit <- d_stack$train(d_task_train)
  d_preds <- d_fit$predict(d_task_valid)

  ## return output as list
  preds <- data.table::as.data.table(list(
    g_preds, e_preds, m_preds, b_preds,
    d_preds
  ))
  data.table::setnames(preds, c(
    "g_preds", "e_preds", "m_preds", "b_preds",
    "d_preds"
  ))
  return(list(
    g_fit = g_fit, e_fit = e_fit, m_fit = m_fit, b_fit = b_fit,
    d_fit = d_fit, preds = preds
  ))
}
