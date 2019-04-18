#' Parameter for the pure mediated effect under binary stochastic interventions
#'
#' Parameter definition...
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @importFrom tmle3 Param_base
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_esteqn_medshift, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_medshift <- R6::R6Class(
  classname = "Param_medshift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::Param_base,
  public = list(
    initialize = function(observed_likelihood,
                          intervention_list_treatment,
                          intervention_list_control,
                          outcome_node = "Y") {
      # copied from parameter definition for ATE
      super$initialize(observed_likelihood, list(),
                       outcome_node = outcome_node)
      private$.cf_likelihood_treatment <-
        CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <-
        CF_Likelihood$new(observed_likelihood, intervention_list_control)

      # wrap up intervention nodes and counterfactuals for Phi task
      intervention_list <- list(treatment = intervention_list_treatment,
                                control = intervention_list_control)
      cf_likelihoods_list <- list(treatment = self$cf_likelihood_treatment,
                                  control = self$cf_likelihood_control)

      # compute extra regressions
      tmle_task <- observed_likelihood$training_task
      
      # build task for and compute e(A|Z,W) re-parameterization
      e_task <- build_e_task(observed_likelihood = observed_likelihood)
      # NOTE: needs something about revere sl3 tasks for eventual CV-TMLE

      # build task for and compute Phi nuisance parameter
      phi_task <- build_phi_task(observed_likelihood = observed_likelihood,
                                 cf_likelihoods = cf_likelihoods_list,
                                 intervention_list = intervention_list,
                                 fold_number = "full")

    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment),
                                  names(self$intervention_list_control))

      pA <- self$observed_likelihood$get_likelihoods(tmle_task,
                                                     intervention_nodes,
                                                     fold_number)
      cf_pA_treatment <-
        self$cf_likelihood_treatment$get_likelihoods(tmle_task,
                                                     intervention_nodes,
                                                     fold_number)
      cf_pA_control <-
        self$cf_likelihood_control$get_likelihoods(tmle_task,
                                                   intervention_nodes,
                                                   fold_number)

      HA <- (cf_pA_treatment - cf_pA_control) / pA
      return(list(Y = HA))
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment),
                                  names(self$intervention_list_control))

      # clever_covariates happens here but this is repeated computation
      HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]

      # todo: make sure we support updating these params
      pA <-
        self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes,
                                                 fold_number)
      cf_pA_treatment <-
        self$cf_likelihood_treatment$get_likelihoods(tmle_task,
                                                     intervention_nodes,
                                                     fold_number)
      cf_pA_control <-
        self$cf_likelihood_control$get_likelihoods(tmle_task,
                                                   intervention_nodes,
                                                   fold_number)

      # counterfactual tasks 
      cf_task_treatment <- self$cf_likelihood_treatment$cf_tasks[[1]]
      cf_task_control <- self$cf_likelihood_control$cf_tasks[[1]]

      Y <- tmle_task$get_tmle_node(self$outcome_node)

      EY <- self$observed_likelihood$get_likelihood(tmle_task,
                                                    self$outcome_node,
                                                    fold_number)
      EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment,
                                                     self$outcome_node,
                                                     fold_number)
      EY0 <- self$observed_likelihood$get_likelihood(cf_task_control,
                                                     self$outcome_node,
                                                     fold_number)

      # parameter and influence function
      psi <- mean(EY1 - EY0)
      IC <- HA * (Y - EY) + (EY1 - EY0) - psi

      # output
      result <- list(psi = psi, IC = IC)
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf(
        "E[%s_{%s}]", self$outcome_node,
        self$cf_likelihood$name
      )
      return(param_form)
    },
    cf_likelihood = function() {
      return(private$.cf_likelihood)
    },
    intervention_list = function() {
      return(self$cf_likelihood$intervention_list)
    },
    update_nodes = function() {
      return(self$outcome_node)
    }
  ),
  private = list(
    .type = "IPSI_binary",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    # function to build regression task for Phi nuisance parameter
    build_phi_task = function(observed_likelihood, cf_likelihoods,
                              intervention_list, fold_number = "full") {
      # get training task from observed data
      tmle_task <- observed_likelihood$training_task

      # build intervention nodes object as names of intervention list
      intervention_nodes <- union(names(intervention_list$treatment),
                                  names(intervention_list$control))

      # TODO: make sure we support updating these params
      pA <- observed_likelihood$get_likelihoods(tmle_task, intervention_nodes,
                                                fold_number)
      cf_pA_treatment <-
        cf_likelihoods$treatment$get_likelihoods(tmle_task, intervention_nodes,
                                                 fold_number)
      cf_pA_control <-
        cf_likelihoods$control$get_likelihoods(tmle_task, intervention_nodes,
                                               fold_number)

      # extract counterfactual tasks from input counterfactual likelihoods
      cf_task_treatment <- cf_likelihoods$treatment$cf_tasks[[1]]
      cf_task_control <- cf_likelihoods$control$cf_tasks[[1]]

      # counterfactual outcomes m(A=1,Z,W) and m(A=0,Z,W)
      EY1 <- observed_likelihood$get_likelihood(cf_task_treatment,
                                                self$outcome_node,
                                                fold_number)
      EY0 <- observed_likelihood$get_likelihood(cf_task_control,
                                                self$outcome_node,
                                                fold_number)

      # get difference for Phi reduced-dimension regression
      phi_outcome <- EY1 - EY0
      phi_data <- as.data.table(tmle_task$npsem[["W"]]$variables,
                                phi_outcome = phi_outcome)
      phi_task <- sl3_Task$new(data = phi_data,
                               outcome = "phi_outcome",
                               covariates =  tmle_task$npsem[["W"]]$variables)
      return(phi_task)
    },
    # function to build regression task for e(A|Z,W) nuisance parameter
    build_e_task = function(observed_likelihood) {
      # make task for e(A|Z,W) and output
      e_data <- observed_likelihood$training_task$internal_data
      e_task <- sl3_Task$new(data = e_data,
                             outcome = tmle_task$npsem[["A"]]$variables,
                             covariates = c(tmle_task$npsem[["Z"]]$variables,
                                            tmle_task$npsem[["W"]]$variables))
      return(e_task)
    }
  )
)
