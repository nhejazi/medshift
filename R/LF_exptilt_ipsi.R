#' Exponential Tilt Likelihood Factor for Incremental Propensity Scores
#'
#' ...
#'
#' @references
#' \describe{
#'   \item{"Nonparametric Causal Effects Based on Incremental Propensity Score
#'         Interventions."}{Edward's JASA paper...}
#'   \item{"Causal Mediation Analysis for Stochastic Interventions"}{Díaz, Iván
#'         and Hejazi, Nima S (2019). submitted.}
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#'
#' @family Likelihood objects
#'
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_exptilt_ipsi, name, type = "density", likelihood_base,
#'     shift_delta, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node
#'           name in the nodes specified by \code{\link{tmle3_Task}$npsem}.
#'     }
#'     \item{\code{likelihood_base}}{The trained \code{\link{likelihood}}
#'           object, for use in generating a re-scaled likelihood factor.
#'     }
#'     \item{\code{shift_delta}}{\code{numeric}, specification of the magnitude
#'           of the desired shift (a multiplier for the propensity score).
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'   }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{likelihood_base}}{The trained \code{\link{likelihood}}
#'           object, for use in generating a re-scaled likelihood factor.
#'     }
#'     \item{\code{shift_delta}}{\code{numeric}, specification of the magnitude
#'           of the desired shift (a multiplier for the propensity score).
#'     }
#'     \item{\code{...}}{Additional arguments passed to the base class.
#'     }
#'   }
#'
#' @export
#
LF_exptilt_ipsi <- R6::R6Class(
  classname = "LF_exptilt_ipsi",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::LF_base,
  public = list(
    initialize = function(name, likelihood_base, shift_delta, ...) {
      super$initialize(name, ..., type = "density")
      private$.likelihood_base <- likelihood_base
      private$.shift_delta <- shift_delta
    },
    get_mean = function(tmle_task, fold_number) {
      stop(paste("get_mean not supported for", class(self)[1]))
    },
    get_density = function(tmle_task, fold_number = "full") {
      # create treatment and control tasks for intervention conditions
      treatment_task <-
        tmle_task$generate_counterfactual_task(uuid = uuid::UUIDgenerate(),
                                               new_data = data.table(A = 1))
      control_task <-
        tmle_task$generate_counterfactual_task(uuid = uuid::UUIDgenerate(),
                                               new_data = data.table(A = 0))

      # get shifted likelihood values for g(A,W)
      g1 <- self$likelihood_base$get_likelihood(treatment_task, "A",
                                                fold_number)
      g0 <- self$likelihood_base$get_likelihood(control_task, "A",
                                                fold_number)
      g_delta <- (exp(self$shift_delta * tmle_task$get_tmle_node("A")) *
                  self$likelihood_base$get_likelihood(tmle_task, "A",
                                                      fold_number)) /
        (self$shift_delta * g1 + g0)

      # return counterfactual likelihood for shifted propensity score
      cf_likelihood <- g_delta
      return(cf_likelihood)
    },
    cf_values = function(tmle_task) {
      cf_values <- rep(NA, tmle_task$nrow)
      return(cf_values)
      #stop(paste("cf_values is undefined for", class(self)[1]))
    }
  ),
  active = list(
    likelihood_base = function() {
      return(private$.likelihood_base)
    },
    shift_delta = function() {
      return(private$.shift_delta)
    }
  ),
  private = list(
    .name = NULL,
    .likelihood_base = NULL,
    .shift_delta = NULL
  )
)
