#' Nonparametric estimation of the stochastic and stochastic-interventional
#' (in)direct effects
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar rectangular data
#'  object of covariates, corresponding to measured baseline confounders.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  causal parameter is defined by an intervention that shifts this quantity
#'  through an (incremental) exponential tilt or a modified treatment policy.
#' @param L A \code{numeric} vector corresponding to an intermediate confounder
#'  affected by treatment (on the causal pathway between the treatment A, the
#'  mediators Z, and outcome Y. When set to \code{NULL} (default), stochastic
#'  (in)direct effects are estimated; when specified, stochastic-interventional
#'  effects are estimated instead, as the former are unidentifiable.
#' @param Z A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
#'  similar rectanguilar data object corresponding to a set of mediators that
#'  sit on the causal pathway between the treatment A and the outcome Y.
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param ids A \code{numeric} vector of observation-level IDs, allowing for
#'  observational units to be related through a hierarchical structure. The
#'  default is to assume all units are independent. When repeated measures are
#'  included, both the cross-validation scheme for nuisance paramter estimation
#'  and inferential procedures are accordingly adjusted.
#' @param delta A \code{numeric} value indicating the magnitude of shift in the
#'  treatment that defines the counterfactual quantity of interest. In the case
#'  of binary treatments, this defines the incremental propensity score shift,
#'  which acts as a multiplier of the odds with which an observational unit
#'  receives treatment (EH Kennedy, 2018; <doi:10.1080/01621459.2017.1422737>).
#'  When the treatment is quantitative, the defines the degree of shift for a
#'  modified treatment policy but non-binary treatments are not  yet supported.
#' @param g_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, used to fit a propensity score model.
#' @param e_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, used to fit a reparametrized model
#'  for the mediator(s) that is equivalent to a propensity score conditioning
#'  on the mediators.
#' @param b_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, used to fit a propensity score model
#'  for the intermediate confounder, adjusting for the baseline covariates and
#'  the treatment both. This is irrelevant for stochastic (in)direct effects,
#'  and so is ignored when L is left to its default of \code{NULL}.
#' @param d_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, used to fit a propensity score model
#'  for the intermediate confounder, adjusting for the baseline covariates, the
#'  treatment, and the mediators. This is irrelevant for stochastic (in)direct
#'  effects, and so is ignored when L is left to its default of \code{NULL}.
#' @param m_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, used to fit an outcome regression.
#' @param phi_learners A \code{\link[sl3]{Stack}} (or other learner class that
#'  inherits from \code{\link[sl3]{Lrnr_base}}), containing a single or set of
#'  instantiated learners from \pkg{sl3}, used to fit a nuisance regression of
#'  a "derived" pseudo-outcome (contrast of outcome regressions under different
#'  treatment conditions), conditional on baseline covariates. This is of the
#'  following form phi(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) | W). This is
#'  irrelevant for the stochastic-interventional effects (when intermediate
#'  confounding is present) and, so, is ignored when L is not \code{NULL}.
#' @param estimator The desired estimator of the direct or indirect effect to
#'  For the stochastic (in)direct effects, parametric substitution and inverse
#'  probability weighted estimators, alongside both efficient one-step and TML
#'  estimators, are available. For the stochastic-interventional (in)direct
#'  effects, only the latter two, which are compatible with machine learning
#'  estimation of nuisance functionals and incorporate cross-fitting to do so,
#'  are available. The one-step estimator is the default in all cases.
#' @param estimator_args A \code{list} of extra arguments to be passed (via
#'  \code{...}) to the internal function for the desired estimator. The default
#'  is chosen so as to allow the number of folds used in computing the one-step
#'  estimator to be easily tweaked. Refer to the documentation for functions
#'  \code{\link{est_onestep}}, \code{\link{est_ipw}}, and \code{\link{est_sub}}
#'  for details on what other arguments may be specified.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom origami make_folds cross_validate
#' @importFrom sl3 Lrnr_glm
#' @importFrom tmle3 tmle3
#' @importFrom assertthat assert_that
#'
#' @export
medshift <- function(W,
                     A,
                     L = NULL,
                     Z,
                     Y,
                     ids = seq_along(Y),
                     delta,
                     g_learners = sl3::Lrnr_glm$new(),
                     e_learners = sl3::Lrnr_glm$new(),
                     m_learners = sl3::Lrnr_glm$new(),
                     phi_learners = sl3::Lrnr_glm$new(),
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

  # NOTE: no static interventions -- need a valid IPSI stochastic intervention
  assertthat::assert_that(delta > 0 && delta < Inf)

  #browser()
  # construct input data structure
  if (is.null(L)) {
    data <- data.table::as.data.table(cbind(Y, Z, A, W, ids))
  } else {
    data <- data.table::as.data.table(cbind(Y, Z, L, A, W, ids))
  }
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste("Z", seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", z_names, "A", w_names, "ids"))

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
    node_list <- list(W = w_names, A = "A", Z = z_names, Y = "Y", id = "ids")
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

  # lazily create output as ad-hoc S3 class, except for tmle3 output
  if (estimator != "tmle") {
    class(est_out) <- "medshift"
  }
  return(est_out)
}
