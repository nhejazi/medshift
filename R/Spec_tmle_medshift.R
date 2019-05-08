#' Defines TML Estimator for mediation with incremental propensity interventions
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}
#' See TODO notes for places generalization can be added
#'
#' @importFrom R6 R6Class
#' @importFrom tmle3 tmle3_Spec define_lf tmle3_Update Targeted_Likelihood
#'  Param_TSM
#'
#' @export
#
tmle3_Spec_medshift <- R6::R6Class(
  classname = "tmle3_Spec_medshift",
  portable = TRUE,
  class = TRUE,
  inherit = tmle3_Spec,
  public = list(
    initialize = function(shift_fxn = shift_ipsi,
                              shift_val = 0,
                              e_lrnrs,
                              phi_lrnrs,
                              ...) {
      options <- list(
        shift_fxn = shift_fxn,
        delta_ipsi = shift_val,
        e_lrnrs = e_lrnrs,
        phi_lrnrs = phi_lrnrs,
        ...
      )
      do.call(super$initialize, options)
    },
    make_params = function(tmle_task, likelihood) {
      # TODO: export and use sl3:::get_levels
      A_vals <- tmle_task$get_tmle_node("A")

      # unwrap internalized arguments
      shift_fxn <- self$options$shift_fxn
      delta_ipsi <- self$options$delta_ipsi
      e_lrnrs <- self$options$e_lrnrs
      phi_lrnrs <- self$options$phi_lrnrs

      # derived likelihood factors: e(A,W) and phi(W)
      lf_e <- tmle3::define_lf(
        tmle3::LF_derived, "E", e_lrnrs, likelihood,
        make_e_task
      )
      lf_phi <- tmle3::define_lf(
        tmle3::LF_derived, "phi", phi_lrnrs,
        likelihood, make_phi_task
      )
      likelihood$add_factors(lf_e)
      likelihood$add_factors(lf_phi)

      # define incremental propensity score intervention
      intervention <- tmle3::define_lf(LF_exptilt_ipsi,
        name = "A",
        original_lf = likelihood$factor_list[["A"]],
        likelihood_base = likelihood, # initialized likelihood
        shift_fxn, # shift function (from user)
        shift_factor = delta_ipsi # magnitude of shift multiplier
      )
      shifted_mean <- tmle3::Param_TSM$new(likelihood, intervention)

      # output should be a list
      tmle_params <- list(shifted_mean)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

################################################################################

#' Outcome under Shifted Treatment
#'
#' O = (W, A, Y)
#' W = Covariates
#' A = Treatment (binary or categorical)
#' Y = Outcome (binary or bounded continuous)
#'
#' @param shift_fxn A \code{function} defining the type of shift to be applied
#'  to the treatment. For an example, see \code{shift_additive}.
#' @param shift_fxn_inv A \code{function} defining the inverse of the type of
#'  shift to be applied to the treatment. For an example, see
#'  \code{shift_additive_inv}.
#' @param shift_val A \code{numeric}, specification of the magnitude of the
#'  desired shift (on the level of the treatment). This is a value passed to
#'  the \code{function}s above for modulating the treatment.
#' @param max_shifted_ratio A \code{numeric} value indicating the maximum
#'  tolerance for the ratio of the counterfactual and observed intervention
#'  densities. In particular, the shifted value of the intervention is assigned
#'  to a given observational unit when the ratio of counterfactual intervention
#'  density to the observed intervention density is below this value.
##' @param ... Additional arguments (currently unused).
#'
#' @importFrom sl3 make_learner Lrnr_mean
#'
#' @export
#
tmle_shift <- function(shift_fxn = shift_additive_bounded,
                       shift_fxn_inv = shift_additive_bounded_inv,
                       shift_val = 1, max_shifted_ratio = 2, ...) {
  # TODO: unclear why this has to be in a factory function
  tmle3_Spec_shift$new(
    shift_fxn, shift_fxn_inv,
    shift_val, max_shifted_ratio,
    ...
  )
}
