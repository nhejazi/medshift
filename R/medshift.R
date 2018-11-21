utils::globalVariables(c("..eif_component_names"))

#' Nonparametric estimation of direct effects under mediation
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
#' @param estimator The desired estimator of the natural direct effect to be
#'  computed. Currently, choices are limited to a substitution estimator, a
#'  re-weighted estimator, and an efficient one-step estimator. The interested
#'  user should consider consulting Díaz & Hejazi (2018+) for a comparative
#'  investigation of each of these estimators.
#'
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
                     delta = 0.5,
                     g_lrnrs =
                       sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                     e_lrnrs =
                       sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                     m_lrnrs = sl3::Lrnr_glm_fast$new(),
                     phi_lrnrs = sl3::Lrnr_glm_fast$new(),
                     estimator = c(
                       "efficient", "substitution",
                       "reweighted"
                     )) {
  # set defaults
  estimator <- match.arg(estimator)

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(W)[2]), sep = "_")
  z_names <- paste("Z", seq_len(dim(Z)[2]), sep = "_")
  data.table::setnames(data, c("Y", z_names, "A", w_names))

  # SUBSTITUTION ESTIMATOR
  if (estimator == "substitution") {
    # fit regression for incremental propensity score intervention
    g_out <- fit_g_mech(
      data = data, delta = delta,
      lrnr_stack = g_lrnrs, w_names = w_names
    )

    # fit regression for incremental propensity score intervention
    m_out <- fit_m_mech(
      data = data, lrnr_stack = m_lrnrs,
      z_names = z_names, w_names = w_names
    )

    # build estimate
    g_shifted_A1 <- g_out$g_est$g_pred_shifted
    g_shifted_A0 <- 1 - g_shifted_A1
    m_pred_A1 <- m_out$m_pred$m_pred_A1
    m_pred_A0 <- m_out$m_pred$m_pred_A0
    estim_sub <- mean(m_pred_A0 * g_shifted_A0) +
      mean(m_pred_A1 * g_shifted_A1)

    # output
    est_out <- estim_sub

    # REWEIGHTED ESTIMATOR
  } else if (estimator == "reweighted") {
    # fit regression for incremental propensity score intervention
    g_out <- fit_g_mech(
      data = data, delta = delta,
      lrnr_stack = g_lrnrs, w_names = w_names
    )

    # fit clever regression for treatment, conditional on mediators
    e_out <- fit_e_mech(
      data = data, lrnr_stack = e_lrnrs,
      z_names = z_names, w_names = w_names
    )

    # stabilize weights in AIPW by dividing by sample average since E[g/e] = 1
    g_shifted <- g_out$g_est$g_pred_shifted
    e_pred <- e_out$e_est$e_pred
    mean_aipw <- mean(g_shifted / e_pred)

    # build estimate
    estim_re <- mean(((g_shifted / e_pred) / mean_aipw) * data$Y)

    # output
    est_out <- estim_re

    # EFFICIENT ESTIMATOR
  } else if (estimator == "efficient") {
    # use origami to perform CV-SL, fitting each EIF component per fold
    # and evaluating only on the holdouts
    eif_component_names <- c("Dy", "Da", "Dzw")
    folds <- origami::make_folds(data)
    cv_eif_results <- origami::cross_validate(
      cv_fun = cv_eif,
      folds = folds,
      data = data,
      delta = delta,
      lrnr_stack_g = g_lrnrs,
      lrnr_stack_e = e_lrnrs,
      lrnr_stack_m = m_lrnrs,
      lrnr_stack_phi = phi_lrnrs,
      z_names = z_names,
      w_names = w_names,
      use_future = FALSE,
      .combine = FALSE
    )
    D_obs <- lapply(cv_eif_results[[1]], function(x) {
      D_obs_fold <- rowSums(x[, ..eif_component_names])
    })
    estim_eff <- mean(do.call(c, D_obs))

    # output
    est_out <- estim_eff
  }
  # final output
  return(est_out)
}
