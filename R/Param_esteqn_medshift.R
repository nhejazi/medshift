#' Parameter for efficient estimation equation-based estimator
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
#'     \item{\code{cf_likelihood}}{the counterfactual likelihood for this treatment
#'     }
#'     \item{\code{intervention_list}}{A list of objects inheriting from \code{\link{LF_base}}, representing the intervention
#'     }
#' }
#' @export
#
Param_esteqn_medshift <- R6::R6Class(
  classname = "Param_esteqn_medshift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list, ...,
                              outcome_node = "Y") {
      super$initialize(observed_likelihood, ..., outcome_node = outcome_node)
      private$.cf_likelihood <- make_CF_Likelihood(
        observed_likelihood,
        intervention_list
      )
    },
    clever_covariates = function(tmle_task = NULL, cv_fold = -1) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- names(self$intervention_list)

      # NOTE: empty list to avoid TMLE updates (trick)
      return(list())
    },
    estimates = function(tmle_task = NULL, cv_fold = -1) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      # todo: extend for stochastic
      cf_task <- self$cf_likelihood$cf_tasks[[1]]


      Y <- tmle_task$get_tmle_node(self$outcome_node)


      # clever_covariates happen here, but this is repeated computation
      HA <- self$clever_covariates(tmle_task, cv_fold)[[self$outcome_node]]

      # clever_covariates happen here, and this is repeated computation
      EYA <-
        unlist(self$observed_likelihood$get_likelihood(
          tmle_task,
          self$outcome_node,
          cv_fold
        ),
        use.names = FALSE
        )

      # clever_covariates happen here, and this is repeated computation
      EY1 <-
        unlist(self$cf_likelihood$get_likelihood(
          cf_task,
          self$outcome_node, cv_fold
        ),
        use.names = FALSE
        )

      # TODO: integrate unbounding logic into likelihood class
      variable_type <- tmle_task$npsem[[self$outcome_node]]$variable_type
      if ((variable_type$type == "continuous") &&
        (!is.na(variable_type$bounds))) {
        bounds <- variable_type$bounds
        scale <- bounds[2] - bounds[1]
        shift <- bounds[1]
        EYA <- EYA * scale + shift
        EY1 <- EY1 * scale + shift
      }

      # todo: separate out psi
      # todo: make this a function of f(W)
      psi <- mean(EY1)
      IC <- HA * (Y - EYA) + EY1 - psi

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
    .type = "TSM_mediation",
    .cf_likelihood = NULL
  )
)
