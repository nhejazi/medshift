#' Efficient estimator
#'
#' @param task ...
#' @param delta_shift ...
#' @param lrnr_stack_g ...
#' @param lrnr_stack_e ...
#' @param lrnr_stack_m ...
#
estim_eff <- function(task, delta_shift,
                      lrnr_stack_g,
                      lrnr_stack_e,
                      lrnr_stack_m) {
  # fit propensity score
  g_out <- fit_g_mech(
    task = task, delta_shift = delta_shift,
    lrnr_stack = lrnr_stack_g
  )
  g_pred_shifted <- g_out$g_est$g_shifted


  # compute the clever intervention density for the mediator
  e_est <- fit_e_mech(task = task, lrnr_stack = lrnr_stack_e)
  e_pred <- e_est$e_pred


  # compute estimator for the outcome regression
  m_out <- fit_m_mech(task = task, lrnr_stack = lrnr_stack_m)
  m_pred <- as.numeric(unlist(m_out$m_pred))


  # construct efficient influence function components
  ## 1) D_ZW:
  # TODO: FIX THIS, VERY BAD
  Dzw <- 0
  for (a in unique(task$data$A)) {
    m_this_a <- m_pred[task$data$A == a]
    g_this_a <- g_pred_shifted[task$data$A == a]
    Dzw_this_a <- mean(m_this_a * g_this_a)
    Dzw <- Dzw + Dzw_this_a
  }

  ## 2) D_A:
  ### get already fit m regression SL
  m_fit <- m_out$m_fit_sl 
  names_W <- names(task$data)[stringr::str_detect(names(task$data), "W")]
  names_Z <- names(task$data)[stringr::str_detect(names(task$data), "Z")]
  ### extract data and make sl3 Tasks
  m_data_A1 <- copy(m_fit$training_task$data)
  m_data_A1[, A := rep(1, .N)]
  m_task_A1 <- sl3::sl3_Task$new(data = m_data_A1,
                               covariates = c(names_Z, "A", names_W),
                               outcome = "Y")
  m_data_A0 <- copy(m_fit$training_task$data)
  m_data_A0[, A := rep(0, .N)]
  m_task_A0 <- sl3::sl3_Task$new(data = m_data_A0,
                               covariates = c(names_Z, "A", names_W),
                               outcome = "Y")
  ### fit m regressions on counterfactual data and compute nuisance parameter
  m_pred_A1 <- m_fit$predict(task_A1)
  m_pred_A0 <- m_fit$predict(task_A0)
  phi_w <- mean(m_pred_A1 - m_pred_A0)

  ### need g(1|W)
  g_fit <- g_out$g_fit
  g_data_A1 <- copy(g_fit$training_task$data)
  g_data_A1[, A := rep(1, .N)]
  g_task_A1 <- sl3::sl3_Task$new(data = g_data_A1,
                               covariates = names_W,
                               outcome = "A")
  g1w <- g_fit$predict(g_task_A1)

  ###
  Da_num <- delta_shift * phi_w * (task$data$A - g1w)
  Da_denom <- (delta_shift * g1w + 1 - g1w)^2
  Da <- Da_num / Da_denom

  ## 3) D_Y:
  Dy <- (g_pred_shifted / e_pred) * (task$data$Y - m_pred)


  # build the re-weighted estimator and return as output
  theta_eff <- mean(Dzw + Da + Dy)
  return(theta_eff)
}

