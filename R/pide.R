#' One-step or TML estimation of the population intervention direct effect
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar corresponding to a
#'  set of baseline covariates.
#' @param A A \code{numeric} vector corresponding to a treatment variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Z A \code{numeric} vector, \code{matrix}, \code{data.frame}, or
#'  similar corresponding to a set of mediators (on the causal pathway between
#'  the intervention A and the outcome Y).
#' @param Y A \code{numeric} vector corresponding to an outcome variable.
#' @param ids A \code{numeric} vector of observation-level IDs, allowing for
#'  observational units to be related through a hierarchical structure. The
#'  default is to assume all units are IID. When repeated IDs are included,
#'  both the cross-validation procedures used for estimation and inferential
#'  procedures respect these IDs.
#' @param delta A \code{numeric} value indicating the degree of shift in the
#'  intervention to be used in defining the causal quantity of interest. In the
#'  case of binary interventions, this takes the form of an incremental
#'  propensity score shift, acting as a multiplier of the odds with which a
#'  unit receives the intervention (EH Kennedy, 2018, JASA;
#'  <doi:10.1080/01621459.2017.1422737>).
#' @param estimator The desired estimator of the natural direct effect to be
#'  computed. Currently, choices are limited to a substitution estimator, a
#'  re-weighted estimator, a one-step estimator, and a targeted minimum loss
#'  estimator.
#' @param ci_level A \code{numeric} indicating the desired coverage level of
#'  the confidence interval to be computed.
#' @param ... Additional arguments passed to \code{\link{medshift}}. Consult
#'  the documentation of that function for details.
#'
#' @export
pide <- function(W,
                 A,
                 Z,
                 Y,
                 ids = seq(1, length(Y)),
                 delta,
                 estimator = c("onestep", "tmle"),
                 ci_level = 0.95,
                 ...) {
  # set default estimator
  estimator <- match.arg(estimator)

  # bounds for confidence interval
  ci_norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - ci_level) / 2))

  # compute mean of outcome and EIF
  EY <- mean(Y)
  eif_EY <- Y - EY

  # pass to medshift for common term in effect decomposition
  est <- medshift(W = W, A = A, Z = Z, Y = Y, ids = ids, delta = delta,
                  estimator = estimator, ...)
  eif_medshift <- est$eif

  # reduce EIFs if there are repeated measures
  if (length(unique(ids)) < length(Y)) {
    eif_EY_combined <- by(eif_EY, as.numeric(ids), mean, simplify = FALSE)
    eif_EY <- unname(do.call(c, eif_EY_combined))
    eif_medshift_combined <- by(eif_medshift, as.numeric(ids), mean,
                                simplify = FALSE)
    eif_medshift <- unname(do.call(c, eif_medshift_combined))
  }

  # compute EIF for direct effect
  eif_de <- eif_EY - eif_medshift
  se_eif_de <- sqrt(var(eif_de) / length(eif_de))

  # compute parameter estimate for direct effect
  param_de <- EY - est$theta

  # direct effect parameter and inference
  param_ci <- param_de + ci_norm_bounds * se_eif_de
  out <- c(param_ci[1], param_de, param_ci[2], se_eif_de)
  names(out) <- c("lwr_ci", "pide_est", "upr_ci", "std_err")
  return(out)
}
