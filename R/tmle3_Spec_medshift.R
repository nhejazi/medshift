#' Defines a TML Estimator for Outcome under Joint Stochastic Intervention on
#' Treatment and Mediator
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'
#' @export
#
tmle3_Spec_medshift <- R6::R6Class(
  classname = "tmle3_Spec_medshift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(shift_type = "exptilt", delta = 0, ...) {
      options <- list(
        shift_type = shift_type,
        delta_shift = delta,
        ...
      )
      do.call(super$initialize, options)
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      npsem <- stochastic_mediation_npsem(node_list)
      tmle_task <- tmle3_Task$new(data, npsem, variable_types)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      likelihood <- stochastic_mediation_likelihood(tmle_task, learner_list)
      return(likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      # add derived likelihood factors to targeted likelihood object
      lf_e <- tmle3::define_lf(
        tmle3::LF_derived, "E", cv_hal_binary_lrnr,
        targeted_likelihood, make_e_task
      )
      lf_phi <- tmle3::define_lf(
        tmle3::LF_derived, "phi", cv_hal_contin_lrnr,
        targeted_likelihood, make_phi_task
      )
      targeted_likelihood$add_factors(lf_e)
      targeted_likelihood$add_factors(lf_phi)

      # compute a tmle3 "by hand"
      tmle_params <- tmle3::define_param(Param_medshift, targeted_likelihood,
        shift_param = self$options$delta_shift
      )
      tmle_params <- list(tmle_params)
      return(tmle_params)
    },
    make_updater = function() {
      updater <- tmle3_Update$new(one_dimensional = TRUE,
                                  constrain_step = TRUE,
                                  maxit = 1e4,
                                  delta_epsilon = 1e-4,
                                  cvtmle = TRUE)
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Outcome under Joint Stochastic Intervention on Treatment and Mediator
#'
#' O = (W, A, Z, Y)
#' W = Covariates (possibly multivariate)
#' A = Treatment (binary or categorical)
#' Z = Mediators (binary or categorical; possibly multivariate)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param shift_type A \code{character} defining the type of shift to be applied
#'  to the treatment. By default, this is an exponential tilt intervention.
#' @param delta A \code{numeric}, specification of the magnitude of the
#'  desired shift.
##' @param ... Additional arguments (currently unused).
#'
#' @export
#
tmle_medshift <- function(shift_type = "exptilt",
                          delta = 1, ...) {
  # this is a factory function
  tmle3_Spec_medshift$new(
    shift_type, delta,
    ...
  )
}
