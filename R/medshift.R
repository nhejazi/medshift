#' Nonparametric estimation of decomposition term for causal mediation analysis
#' with stochastic interventions
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar corresponding to a
#'  set of baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Z A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
#'  similar corresponding to a set of mediators (on the causal pathway between
#'  the intervention A and the outcome Y).
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param delta A \code{numeric} value indicating the degree of shift in the
#'  intervention to be used in defining the causal quantity of interest. In the
#'  case of binary interventions, this takes the form of an incremental
#'  propensity score shift, acting as a multiplier of the probability with which
#'  a given observational unit receives the intervention (EH Kennedy, 2018,
#'  JASA; <doi:10.1080/01621459.2017.1422737>).
#' @param g_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting a model for the propensity
#'  score, i.e., g = P(A | W).
#' @param e_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting a cleverly parameterized
#'  propensity score that includes the mediators, i.e., e = P(A | Z, W).
#' @param m_lrnrs A \code{Stack} object, or other learner class (inheriting from
#'  \code{Lrnr_base}), containing a single or set of instantiated learners from
#'  the \code{sl3} package, to be used in fitting the outcome regression, i.e.,
#'  m(A, Z, W).
#' @param phi_lrnrs A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a reduced regression
#'  useful for computing the efficient one-step estimator, i.e., phi(W) =
#'  E[m(A = 1, Z, W) - m(A = 0, Z, W) | W).
#' @param shift_type A choice of the type of stochastic treatment regime to use
#'  -- either \code{"mtp"} for a modified treatment policy that shifts the
#'  center of the observed intervention distribution by the scalar \code{delta}
#'  or \code{"ipsi"} for an incremental propensity score shift that multiples
#'  the odds of receiving the intervention by the scalar \code{delta}.
#' @param estimator The desired estimator of the natural direct effect to be
#'  computed. Currently, choices are limited to a substitution estimator, a
#'  re-weighted estimator, and an efficient one-step estimator. The interested
#'  user should consider consulting Díaz & Hejazi (2019+) for a comparative
#'  investigation of each of these estimators.
#' @param estimator_args A \code{list} of extra arguments to be passed (via
#'  \code{...}) to the function call for the specified estimator. The default
#'  is so chosen as to allow the number of folds used in computing the AIPW
#'  estimator to be easily tweaked. Refer to the documentation for functions
#'  \code{\link{est_onestep}}, \code{\link{est_ipw}}, and
#'  \code{\link{est_substitution}} for details on what other arguments may be
#'  specified through this mechanism.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table setnames
#' @importFrom origami make_folds cross_validate
#' @importFrom sl3 Lrnr_glm_fast
#' @importFrom stats binomial
#'
#' @export
#
medshift <- function(W,
                     A,
                     Z,
                     Y,
                     delta,
                     g_lrnrs =
                       sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                     e_lrnrs =
                       sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                     m_lrnrs = sl3::Lrnr_glm_fast$new(),
                     phi_lrnrs = sl3::Lrnr_glm_fast$new(),
                     shift_type = c("ipsi", "mtp"),
                     estimator = c(
                       "onestep",
                       "tmle",
                       "sub",
                       "ipw"
                     ),
                     estimator_args = list(cv_folds = 5)) {
  # set defaults
  estimator <- match.arg(estimator)
  shift_type <- match.arg(shift_type)
  estimator_args <- unlist(estimator_args, recursive = FALSE)

  # check whether type of shift is appropriate for intervention node
  if (shift_type == "ipsi") {
    assertthat::assert_that(length(unique(A)) == 2)
    message("Intervention on binary A: incremental propensity score shift")
  } else {
    message("Intervention on continuous A: modified treatment policy")
  }

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste("Z", seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", z_names, "A", w_names))

  if (estimator == "sub") {
    # SUBSTITUTION/G-COMPUTATION ESTIMATOR
    sub_est_args <- list(
      data = data, delta = delta, g_lrnrs = g_lrnrs,
      m_lrnrs = m_lrnrs, w_names = w_names, z_names = z_names,
      shift_type = shift_type, estimator_args
    )
    est_out <- do.call(est_substitution, sub_est_args)
  } else if (estimator == "ipw") {
    # INVERSE PROBABILITY WEIGHTED ESTIMATOR
    ipw_est_args <- list(
      data = data, delta = delta, g_lrnrs = g_lrnrs,
      e_lrnrs = e_lrnrs, w_names = w_names, z_names = z_names,
      shift_type = shift_type, estimator_args
    )
    est_out <- do.call(est_ipw, ipw_est_args)
  } else if (estimator == "onestep") {
    # ONE-STEP ESTIMATOR WITH CROSS-FITTING
    aipw_est_args <- list(
      data = data, delta = delta, g_lrnrs = g_lrnrs,
      e_lrnrs = e_lrnrs, m_lrnrs = m_lrnrs,
      phi_lrnrs = phi_lrnrs, w_names = w_names,
      z_names = z_names, shift_type = shift_type, estimator_args
    )
    est_out <- do.call(est_onestep, aipw_est_args)
  } else if (estimator == "tmle") {
    # TODO: TARGETED MAXIMUM LIKELIHOOD ESTIMATOR
    # NOTE: should be a call to tmle3::fit_tmle
    stop("The TML estimator is currently under development")
  }

  # output as classed list
  class(est_out) <- "medshift"
  return(est_out)
}
