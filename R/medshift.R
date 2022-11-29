utils::globalVariables(c(".N", "est"))

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
#'  modified treatment policy but support for continuous-valued treatments is
#'  not yet implemented.
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
#'  \code{\link{stoch_est_onestep}}, \code{\link{stoch_est_ipw}}, and
#'  \code{\link{stoch_est_sub}} for details on other arguments to specify.
#' @param ci_level A \code{numeric} indicating the confidence interval level.
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr mutate relocate
#' @importFrom data.table as.data.table setnames
#' @importFrom origami make_folds cross_validate
#' @importFrom sl3 Lrnr_glm_fast
#' @importFrom tmle3 tmle3
#'
#' @export
medshift <- function(W,
                     A,
                     L = NULL,
                     Z,
                     Y,
                     ids = seq_along(Y),
                     delta,
                     g_learners = sl3::Lrnr_glm_fast$new(),
                     e_learners = sl3::Lrnr_glm_fast$new(),
                     b_learners = sl3::Lrnr_glm_fast$new(),
                     d_learners = sl3::Lrnr_glm_fast$new(),
                     m_learners = sl3::Lrnr_glm_fast$new(),
                     phi_learners = sl3::Lrnr_glm_fast$new(),
                     estimator = c( "onestep", "tmle", "sub", "ipw"),
                     estimator_args = list(
                       cv_folds = 10L,
                       max_iter = 1e4,
                       step_size = 1e-6
                     ),
                     ci_level = 0.95
                    ) {
  # set defaults
  estimator <- match.arg(estimator)
  estimator_args <- unlist(estimator_args, recursive = FALSE)

  # NOTE: no static interventions: need valid stochastic IPSI
  assertthat::assert_that(delta > 0 && delta < Inf)

  # construct input data structure
  w_names <- paste("W", seq_len(ncol(as.matrix(W))), sep = "_")
  z_names <- paste("Z", seq_len(ncol(as.matrix(Z))), sep = "_")
  if (is.null(L)) {
    data <- data.table::as.data.table(cbind(Y, Z, A, W, ids))
    data.table::setnames(data, c("Y", z_names, "A", w_names, "ids"))
  } else {
    data <- data.table::as.data.table(cbind(Y, Z, L, A, W, ids))
    data.table::setnames(data, c("Y", z_names, "L", "A", w_names, "ids"))
  }

  # NOTE: in this case (without L), we provide several estimators for only the
  #       decomposition term defining the stochastic direct/indirect effects
  if (is.null(L)) {
    # PARAMETRIC SUBSTITUTION ESTIMATOR
    if (estimator == "sub") {
      sub_est_args <- list(
        data = data, delta = delta,
        g_learners = g_learners,
        m_learners = m_learners,
        w_names = w_names, z_names = z_names
      )
      # TODO: add inference based on the bootstrap
      est_out <- do.call(stoch_est_sub, sub_est_args)

    # INVERSE PROBABILITY WEIGHTED ESTIMATOR
    } else if (estimator == "ipw") {
      ipw_est_args <- list(
        data = data, delta = delta,
        g_learners = g_learners, e_learners = e_learners,
        w_names = w_names, z_names = z_names
      )
      # TODO: add inference based on IPW estimating function
      est_out <- do.call(stoch_est_ipw, ipw_est_args)

    # CROSS-FITTED ONE-STEP ESTIMATOR
    } else if (estimator == "onestep") {
      os_est_args <- list(
        data = data, delta = delta,
        g_learners = g_learners, e_learners = e_learners,
        m_learners = m_learners, phi_learners = phi_learners,
        w_names = w_names, z_names = z_names,
        cv_folds = estimator_args[["cv_folds"]]
      )
      est_out <- do.call(stoch_est_onestep, os_est_args)

    # CROSS-VALIDATED TARGETED MINIMUM LOSS ESTIMATOR
    } else if (estimator == "tmle") {
      node_list <- list(
        W = w_names, A = "A", Z = z_names, Y = "Y", id = "ids"
      )
      learner_list <- list(Y = m_learners, A = g_learners)
      tmle_spec <- tmle_medshift(
        delta = delta,
        e_learners = e_learners,
        phi_learners = phi_learners,
        max_iter = estimator_args[["max_iter"]],
        step_size = estimator_args[["step_size"]]
      )
      tml_est <- tmle3::tmle3(tmle_spec, data, node_list, learner_list)

      # reorganize tmle3 output for compatibility across estimator types
      est_out <- tibble::as_tibble(list(
        param = "decomp",
        est = tml_est$estimates[[1]]$psi,
        stderr = sqrt(
          stats::var(as.numeric(tml_est$estimates[[1]]$IC)) / data[, .N]
        )
      ))
      attr(est_out, "type") <- "tmle"
      attr(est_out, "param") <- "stochastic"
      attr(est_out, "eif") <- as.numeric(tml_est$estimates[[1]]$IC)
    }

  # NOTE: next branch handles stochastic-interventional effects, for which
  #       we provide only efficient estimators of direct/indirect effects
  } else {
    # PARAMETRIC SUBSTITUTION ESTIMATOR
    if (estimator == "sub") {
      stop("Parametric substitution estimator not available for
            stochastic-interventional direct/indirect effects")

    # INVERSE PROBABILITY WEIGHTED ESTIMATOR
    } else if (estimator == "ipw") {
      stop("Inverse probability weighted estimator not available for
            stochastic-interventional direct/indirect effects")

    # CROSS-FITTED ONE-STEP ESTIMATOR
    } else if (estimator == "onestep") {
      os_est_args <- list(
        data = data, delta = delta,
        g_learners = g_learners, e_learners = e_learners,
        b_learners = b_learners, d_learner = d_learners,
        m_learners = m_learners,
        w_names = w_names, z_names = z_names,
        cv_folds = estimator_args[["cv_folds"]]
      )
      est_out <- do.call(interv_est_onestep, os_est_args)

    # CROSS-VALIDATED TARGETED MINIMUM LOSS ESTIMATOR
    } else if (estimator == "tmle") {
      # set TMLE arguments and compute estimates
      tmle_est_args <- list(
        data = data, delta = delta,
        g_learners = g_learners, e_learners = e_learners,
        b_learners = b_learners, d_learner = d_learners,
        m_learners = m_learners,
        w_names = w_names, z_names = z_names,
        cv_folds = estimator_args[["cv_folds"]]
      )
      est_out <- do.call(interv_est_tmle, tmle_est_args)
    }
  }

  # compute inference for estimates -- first, need Wald-style multiplier
  wald_mult <- abs(stats::qnorm(p = (1 - ci_level) / 2))

  # inference for stochastic direct/indirect effects
  if (attr(est_out, "param") == "stochastic") {
    est_out <- est_out %>%
      dplyr::mutate(
        ci_lwr = est - wald_mult * stderr,
        ci_upr = est + wald_mult * stderr
      ) %>%
      dplyr::relocate("param", "ci_lwr", "est", "ci_upr", "stderr")

  # inference for stochastic-interventional direct/indirect effects
  } else if (attr(est_out, "param") == "stochastic_interventional") {
    est_out <- est_out %>%
      dplyr::mutate(
        ci_lwr = est - wald_mult * stderr,
        ci_upr = est + wald_mult * stderr
      ) %>%
      dplyr::relocate("param", "ci_lwr", "est", "ci_upr", "stderr")
  }
  return(est_out)
}
