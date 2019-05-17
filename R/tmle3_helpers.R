#' Stochastic Mediation NPSEM
#'
#' @param node_list A \code{list} object specifying the different nodes in the
#'  nonparametric structural equation model.
#' @param variable_types Used to define how variables are handled. No need for
#'  user specification of this argument.
#'
#' @importFrom tmle3 define_node
#'
#' @export
#
stochastic_mediation_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task
  npsem <- list(
    tmle3::define_node("W", node_list$W, variable_type = variable_types$W),
    tmle3::define_node("A", node_list$A, c("W"),
                       variable_type = variable_types$A),
    tmle3::define_node("Z", node_list$Z, c("A", "W"),
                       variable_type = variable_types$Z),
    tmle3::define_node("Y", node_list$Y, c("Z", "A", "W"),
                       variable_type = variable_types$Y, scale = TRUE)
  )
  return(npsem)
}


#' Stochastic Mediation Likelihood Factors
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'  NPSEM for use in constructing elements of TML estimator.
#' @param learner_list A \code{list} specifying which learners are to be applied
#'  for each of the regression tasks required for the TML estimator.
#'
#' @importFrom tmle3 define_lf LF_emp LF_fit Likelihood
#'
#' @export
#
stochastic_mediation_likelihood <- function(tmle_task, learner_list) {
  # covariates
  W_factor <- tmle3::define_lf(tmle3::LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continuous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  # treatment
  A_factor <- tmle3::define_lf(tmle3::LF_fit, "A",
    learner = learner_list[["A"]],
    bound = A_bound
  )

  # outcome
  Y_factor <- tmle3::define_lf(tmle3::LF_fit, "Y",
    learner = learner_list[["Y"]],
    type = "mean"
  )

  # construct and train likelihood
  factor_list <- list(W_factor, A_factor, Y_factor)

  likelihood_def <- tmle3::Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}


#' Make task for derived likelihood factor e(A,W)
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'   NPSEM for use in constructing elements of TML estimator.
#' @param likelihood A trained \code{Likelihood} object from \code{tmle3},
#'  constructed via the helper function \code{stochastic_mediation_likelihood}.
#'
#' @importFrom sl3 sl3_Task
#'
#' @export
#
make_e_task <- function(tmle_task, likelihood) {
  e_data <- tmle_task$internal_data
  e_task <- sl3::sl3_Task$new(
    data = e_data,
    outcome = tmle_task$npsem[["A"]]$variables,
    covariates = c(
      tmle_task$npsem[["Z"]]$variables,
      tmle_task$npsem[["W"]]$variables
    )
  )
  return(e_task)
}

#' Make task for derived likelihood factor phi(W)
#'
#' @param tmle_task A \code{tmle3_Task} object specifying the data and the
#'  NPSEM for use in constructing elements of TML estimator.
#' @param likelihood A trained \code{Likelihood} object from \code{tmle3},
#'  constructed via the helper function \code{stochastic_mediation_likelihood}.
#'
#' @importFrom data.table as.data.table data.table
#' @importFrom uuid UUIDgenerate
#' @importFrom sl3 sl3_Task
#'
#' @export
#
make_phi_task <- function(tmle_task, likelihood) {
  # create treatment and control tasks for intervention conditions
  treatment_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 1)
    )
  control_task <-
    tmle_task$generate_counterfactual_task(
      uuid = uuid::UUIDgenerate(),
      new_data = data.table::data.table(A = 0)
    )

  # create counterfactual outcomes and construct pseudo-outcome
  m1 <- likelihood$get_likelihood(treatment_task, "Y")
  m0 <- likelihood$get_likelihood(control_task, "Y")
  m_diff <- m1 - m0

  # create regression task for pseudo-outcome and baseline covariates
  phi_data <- data.table::as.data.table(list(
    m_diff = m_diff,
    tmle_task$get_tmle_node("W")
  ))
  phi_task <- sl3::sl3_Task$new(
    data = phi_data,
    outcome = "m_diff",
    covariates = tmle_task$npsem[["W"]]$variables,
    outcome_type = "continuous"
  )
  return(phi_task)
}
