#' Nonparametric estimation of the population intervention (in)direct effects
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
#'  propensity score shift, acting as a multiplier of the odds with which a
#'  unit receives the intervention (EH Kennedy, 2018, JASA;
#'  <doi:10.1080/01621459.2017.1422737>).
#' @param g_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'   inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'   instantiated learners from \pkg{sl3}, to be used in fitting the propensity
#'   score, i.e., g = P(A | W).
#' @param e_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'   inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'   instantiated learners from \pkg{sl3}, to be used in fitting a propensity
#'   score that conditions on the mediators, i.e., e = P(A | Z, W).
#' @param m_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'   inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'   instantiated learners from \pkg{sl3}, to be used in fitting the outcome
#'   regression, i.e., m(A, Z, W).
#' @param phi_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, to be used in a regression of a
#'  pseudo-outcome on the baseline covariates, i.e.,
#'  phi(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) | W).
#' @param estimator The desired estimator of the natural direct effect to be
#'  computed. Currently, choices are limited to a substitution estimator, a
#'  re-weighted estimator, a one-step estimator, and a targeted minimum loss
#'  estimator.
#' @param estimator_args A \code{list} of extra arguments to be passed (via
#'  \code{...}) to the function call for the specified estimator. The default
#'  is so chosen as to allow the number of folds used in computing the one-step
#'  estimator to be easily tweaked. Refer to the documentation for functions
#'  \code{\link{est_onestep}}, \code{\link{est_ipw}}, and
#'  \code{\link{est_substitution}} for details on what other arguments may be
#'  specified through this mechanism. For the option \code{"tmle"}, there is
#'  heavy reliance on the architecture provided by \pkg{tmle3}.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom origami make_folds cross_validate
#' @importFrom sl3 Lrnr_glm_fast
#' @importFrom tmle3 tmle3
#' @importFrom stats binomial
#' @importFrom assertthat assert_that
#'
#' @export
medshift <- function(W,
                     A,
                     Z,
                     Y,
                     delta,
                     g_learners =
                       sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                     e_learners =
                       sl3::Lrnr_glm_fast$new(family = stats::binomial()),
                     m_learners = sl3::Lrnr_glm_fast$new(),
                     phi_learners = sl3::Lrnr_glm_fast$new(),
                     estimator = c(
                       "onestep",
                       "tmle",
                       "substitution",
                       "reweighted"
                     ),
                     estimator_args = list(
                       cv_folds = 10,
                       max_iter = 1e4,
                       step_size = 1e-6
                     )) {
  # set defaults
  estimator <- match.arg(estimator)
  estimator_args <- unlist(estimator_args, recursive = FALSE)

  # NOTE: procedure does _not_ support static interventions currently
  assertthat::assert_that(delta > 0 && delta < Inf)

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste("Z", seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", z_names, "A", w_names))

  if (estimator == "substitution") {
    # SUBSTITUTION ESTIMATOR
    sub_est_args <- list(
      data = data, delta = delta,
      g_learners = g_learners, m_learners = m_learners,
      w_names = w_names, z_names = z_names
    )
    est_out <- do.call(est_substitution, sub_est_args)
  } else if (estimator == "reweighted") {
    # INVERSE PROBABILITY RE-WEIGHTED ESTIMATOR
    ipw_est_args <- list(
      data = data, delta = delta,
      g_learners = g_learners, e_learners = e_learners,
      w_names = w_names, z_names = z_names
    )
    est_out <- do.call(est_ipw, ipw_est_args)
  } else if (estimator == "onestep") {
    # CROSS-FITTED ONE-STEP ESTIMATOR
    os_est_args <- list(
      data = data, delta = delta,
      g_learners = g_learners, e_learners = e_learners,
      m_learners = m_learners, phi_learners = phi_learners,
      w_names = w_names, z_names = z_names,
      cv_folds = estimator_args[["cv_folds"]]
    )
    est_out <- do.call(est_onestep, os_est_args)
  } else if (estimator == "tmle") {
    # CROSS-VALIDATED TARGETED MINIMUM LOSS ESTIMATOR
    node_list <- list(W = w_names, A = "A", Z = z_names, Y = "Y")
    learner_list <- list(Y = m_learners, A = g_learners)
    tmle_spec <- tmle_medshift(
      delta = delta,
      e_learners = e_learners,
      phi_learners = phi_learners,
      max_iter = estimator_args[["max_iter"]],
      step_size = estimator_args[["step_size"]]
    )
    est_out <- tmle3::tmle3(tmle_spec, data, node_list, learner_list)
  }

  # lazily create output as S3 class, except for tmle3 output
  if (estimator != "tmle") {
    class(est_out) <- "medshift"
  }
  return(est_out)
}
