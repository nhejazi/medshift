#' Nonparametric estimation of (in)direct effects under stochastic interventions
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
                     delta,
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

  if (estimator == "substitution") {
    # SUBSTITUTION ESTIMATOR
    est_out <- est_sub(
      data = data, delta = delta, g_lrnrs = g_lrnrs,
      m_lrnrs = m_lrnrs, w_names = w_names, z_names = z_names
    )
  } else if (estimator == "reweighted") {
    # RE-WEIGHTED ESTIMATOR
    est_out <- est_re(
      data = data, delta = delta, g_lrnrs = g_lrnrs,
      e_lrnrs = e_lrnrs, w_names = w_names, z_names = z_names
    )

  } else if (estimator == "efficient") {
    # EFFICIENT ESTIMATOR
    est_out <- est_eff(
      data = data, delta = delta, g_lrnrs = g_lrnrs,
      e_lrnrs = e_lrnrs, m_lrnrs = m_lrnrs,
      phi_lrnrs = phi_lrnrs, w_names = w_names,
      z_names = z_names
    )
  }
  # output
  return(est_out)
}
