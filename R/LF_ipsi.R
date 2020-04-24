#' Likelihood Factor for Incremental Propensity Score Interventions
#'
#' @references
#' \describe{
#'   \item{"Nonparametric Causal Effects Based on Incremental Propensity Score
#'         Interventions."}{Kennedy, Edward H (2019). Journal of the American
#'         Statistical Association.
#'         <https://doi.org/10.1080/01621459.2017.1422737>}
#'   \item{"Causal Mediation Analysis for Stochastic Interventions"}{Díaz,
#'         Iván and Hejazi, Nima S (2020). Journal of the Royal Statistical
#'         Society, Series B. <https://doi.org/10.1111/rssb.12362>}
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
#' @return \code{\link[tmle3]{LF_base}} object.
#'
#' @format \code{\link[R6]{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_ipsi, name, type = "density", likelihood_base,
#'                   shift_param, treatment_task, control_task, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{A \code{character}, giving the name of the likelihood
#'           factor. Should match a node name in the nodes specified by the
#'           \code{npsem} slot of \code{\link[tmle3]{tmle3_Task}}.}
#'     \item{\code{likelihood_base}}{A trained \code{\link[tmle3]{Likelihood}}
#'           object, for use in generating a re-scaled likelihood factor.}
#'     \item{\code{shift_param}}{A \code{numeric}, specifying the magnitude of
#'           the desired incremental propensity score shift (a multiplier of
#'           the odds of receiving treatment).}
#'     \item{\code{treatment_task}}{A \code{\link[tmle3]{tmle3_Task}} object
#'           created by setting the intervention to the treatment condition:
#'           do(A = 1).}
#'     \item{\code{control_task}}{A \code{\link[tmle3]{tmle3_Task}} object
#'           created by setting the intervention to the control condition:
#'           do(A = 0).}
#'     \item{\code{...}}{Not currently used.}
#'   }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{likelihood_base}}{A trained \code{\link[tmle3]{Likelihood}}
#'           object, for use in generating a re-scaled likelihood factor.
#'     }
#'     \item{\code{shift_param}}{A \code{numeric}, specifying the magnitude of
#'           the desired incremental propensity score shift (a multiplier of
#'           the odds of receiving treatment).}
#'     \item{\code{treatment_task}}{A \code{\link[tmle3]{tmle3_Task}} object
#'           created by setting the intervention to the treatment condition:
#'           do(A = 1).}
#'     \item{\code{control_task}}{A \code{\link[tmle3]{tmle3_Task}} object
#'           created by setting the intervention to the control condition:
#'           do(A = 0).}
#'     \item{\code{...}}{Additional arguments passed to the base class.}
#'   }
#'
#' @export
LF_ipsi <- R6::R6Class(
  classname = "LF_ipsi",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3::LF_base,
  public = list(
    initialize = function(name, likelihood_base, shift_param,
                              treatment_task, control_task, ...) {
      super$initialize(name, ..., type = "density")
      private$.likelihood_base <- likelihood_base
      private$.shift_param <- shift_param
      private$.treatment_task <- treatment_task
      private$.control_task <- control_task
    },
    get_mean = function(tmle_task, fold_number) {
      stop(paste("get_mean not supported for", class(self)[1]))
    },
    get_density = function(tmle_task, fold_number = "full") {
      # treatment and control tasks for intervention conditions
      treatment_task <- self$treatment_task
      control_task <- self$control_task
      shift_param <- self$shift_param
      likelihood <- self$likelihood_base

      # get likelihood values for counterfactual g(A,W)
      g1 <- likelihood$get_likelihood(treatment_task, "A", fold_number)
      g0 <- likelihood$get_likelihood(control_task, "A", fold_number)

      # compute values for counterfactual (shifted) treatment mechanism
      shift_conditional_treatment <- ifelse(tmle_task$get_tmle_node("A") == 1,
        shift_param, 1
      )
      g_delta <- (shift_conditional_treatment *
        likelihood$get_likelihood(tmle_task, "A", fold_number)) /
        ((shift_param * g1) + g0)

      # return counterfactual likelihood for shifted propensity score
      cf_likelihood <- g_delta
      return(cf_likelihood)
    },
    cf_values = function(tmle_task) {
      cf_values <- rep(NA, tmle_task$nrow)
      return(cf_values)
      # stop(paste("cf_values is undefined for", class(self)[1]))
    }
  ),
  active = list(
    likelihood_base = function() {
      return(private$.likelihood_base)
    },
    shift_param = function() {
      return(private$.shift_param)
    },
    treatment_task = function() {
      return(private$.treatment_task)
    },
    control_task = function() {
      return(private$.control_task)
    }
  ),
  private = list(
    .name = NULL,
    .likelihood_base = NULL,
    .shift_param = NULL,
    .treatment_task = NULL,
    .control_task = NULL
  )
)
