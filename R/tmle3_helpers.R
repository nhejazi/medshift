#' Stochastic Mediation NPSEM
#'
#' @param node_list ...
#' @param variable_types ...
#'
#' @export
#
stochastic_mediation_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task
  npsem <- list(
    define_node("W", node_list$W, variable_type = variable_types$W),
    define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
    define_node("Z", node_list$Z, c("A", "W"), variable_type = variable_types$Z),
    define_node("Y", node_list$Y, c("Z", "A", "W"), variable_type = variable_types$Y, scale = TRUE)
  )
  return(npsem)
}


#' Stochastic Mediation Likelihood Factors
#'
#' @param tmle_task ...
#' @param learner_list ...
#'
#' @export
#
stochastic_mediation_likelihood <- function(tmle_task, learner_list) {
  # covariates
  W_factor <- define_lf(LF_emp, "W")

  # treatment (bound likelihood away from 0 (and 1 if binary))
  A_type <- tmle_task$npsem[["A"]]$variable_type
  if (A_type$type == "continuous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial","categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  # treatment
  A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]],
                        bound = A_bound)

  # mediator re-parametrization
  eA_factor <- define_lf(LF_fit, "eA", learner = learner_list[["eA"]],
                         bound = "mean")

  # outcome
  Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]],
                        type = "mean")

  # construct and train likelihood
  factor_list <- list(W_factor, A_factor, eA_factor, Y_factor)

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}
