#' Confidence Intervals for Stochastic Mediation Parameters
#'
#' Compute confidence intervals for objects of class \code{medshift}, which
#' contain estimates produced by \code{\link{medshift}}.
#'
#' @param object An object of class \code{medshift}, as produced by invoking
#'  \code{\link{medshift}}, for which a confidence interval is to be computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the confidence interval level.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm
#' @importFrom assertthat assert_that
#'
#' @method confint medshift
#'
#' @export
confint.medshift <- function(object,
                             parm = seq_len(object$psi),
                             level = 0.95,
                             ...) {
  # inference is currently limited to the one-step efficient estimator
  # TODO: allow use for TML estimators once impelemented
  assertthat::assert_that(object$type == "onestep")

  # first, let's get Z_(1 - alpha)
  ci_norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

  # compute the EIF variance multiplier for the CI
  # NOTE: the variance value is already scaled by length of observations
  se_eif <- sqrt(object$var)

  # compute the interval around the point estimate
  ci_theta <- ci_norm_bounds * se_eif + object$theta

  # set up output CI object
  ci_out <- c(ci_theta[1], object$theta, ci_theta[2])
  names(ci_out) <- c("lwr_ci", "est", "upr_ci")
  return(ci_out)
}

################################################################################

#' Summary for Stochastic Mediation Parameter Objects
#'
#' Print a convenient summary for objects of \code{S3} class \code{medshift}.
#'
#' @param object An object of class \code{medshift}, as produced by invoking
#'  the function \code{\link{medshift}}, for which a confidence interval is to
#'  be constructed.
#' @param ... Other arguments. Not currently used.
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @importFrom stats confint
#'
#' @method summary medshift
#'
#' @export
summary.medshift <- function(object,
                             ...,
                             ci_level = 0.95) {
  # inference is currently limited to the one-step efficient estimator
  # TODO: allow use for TML estimators once impelemented
  if (object$type == "onestep") {
    # compute confidence interval using the pre-defined method
    ci <- stats::confint(object, level = ci_level)

    # only print useful info about the mean of the efficient influence function
    eif_mean <- format(mean(object$eif), scientific = TRUE)

    # create output table from input object and confidence interval results
    out <- c(round(c(ci, object$var), digits = 6), eif_mean, object$type)
    names(out) <- c(
      "lwr_ci", "param_est", "upr_ci", "param_var", "eif_mean", "estimator"
    )
  } else {
    out <- c(round(object$theta, digits = 6), object$type)
    names(out) <- c(
      "param_est", "estimator"
    )
  }
  print(noquote(out))
}

################################################################################

#' Print Method for Class medshift
#'
#' The \code{print} method for objects of class \code{medshift}.
#'
#' @param x An object of class \code{medshift}.
#' @param ... Other options (not currently used).
#'
#' @method print medshift
#'
#' @export
print.medshift <- function(x, ...) {
  # inference is currently limited to the one-step estimator
  if (x$type == "onestep") {
    print(x[c("theta", "var", "type")])
  } else {
    print(x[c("theta", "type")])
  }
}

###############################################################################

#' Bounding to numerical precision
#'
#' Bounds extreme values to a specified tolerance level, for use with sensitive
#' quantities that must be transformed, e.g., via \code{\link[stats]{qlogis}}.
#'
#' @param vals A \code{numeric} vector of values in the interval [0, 1].
#' @param tol A \code{numeric} indicating the tolerance limit to which extreme
#'  values should be truncated. Realizations of \code{val} less than \code{tol}
#'  are truncated to \code{tol} while those greater than (1 - \code{tol}) are
#'  truncated to (1 - \code{tol}).
#'
#' @keywords internal
bound_precision <- function(vals, tol = 1e-6) {
  vals_bounded <- pmax(pmin(vals, 1 - tol), tol)
  return(vals_bounded)
}

###############################################################################

#' Bounding propensity scores
#'
#' Bounds estimated propensity score values to be within a specified range.
#'
#' @param vals A \code{numeric} vector of values in the interval [0, 1].
#' @param bounds A \code{numeric} vector containing two values, the first being
#'  the minimum allowable value and the second being the maximum allowable for
#'  values appearing in the vector \code{vals} (the previous argument).
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
bound_propensity <- function(vals, bounds = c(0.001, 0.999)) {
  assertthat::assert_that(!(max(vals) > 1 || min(vals) < 0))
  vals_bounded <- pmax(pmin(vals, bounds[2]), bounds[1])
  return(vals_bounded)
}

###############################################################################

#' Map values to the unit interval
#'
#' @param vals A \code{numeric} vector of values to be scaled into the closed
#'  interval [0, 1].
#'
#' @keywords internal
scale_to_unit <- function(vals) {
  # rescale to unit interval
  max_vals <- max(vals)
  min_vals <- min(vals)
  vals_scaled <- (vals - min_vals) / (max_vals - min_vals)

  # store original min/max as attributes
  attr(vals_scaled, "max") <- max_vals
  attr(vals_scaled, "min") <- min_vals
  return(vals_scaled)
}

###############################################################################

#' Map values from the unit interval to their original scale
#'
#' @param vals_in_unit A \code{numeric} vector of values scaled to fall in the
#'  closed unit interval [0, 1] by use of \code{\link{scale_to_unit}}.
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
scale_from_unit <- function(vals_in_unit) {
  # check that input falls in the unit interval
  assertthat::assert_that(min(vals_in_unit) == 0)
  assertthat::assert_that(max(vals_in_unit) == 1)

  # rescale back to the original interval
  max_orig <- attr(vals_in_unit, "max")
  min_orig <- attr(vals_in_unit, "min")
  vals_rescaled <- (vals_in_unit * (max_orig - min_orig)) + min_orig
  return(vals_rescaled)
}

###############################################################################

#' Check TMLE Fluctuation Model
#'
#' @param tilt_mod ...
#' @param tilt_tol ...
#'
#' @importFrom stats coef
check_fluc <- function(tilt_mod, tilt_tol = 10) {
  # number of clever covariates
  n_coefs <- length(stats::coef(tilt_mod))

  # set NA coefficients to zero
  tilt_mod[["coefficients"]][is.na(stats::coef(tilt_mod))] <- rep(0, n_coefs)

  # check sanity of remaining coefficients
  if (!tilt_mod[["converged"]] | abs(max(stats::coef(tilt_mod))) > tilt_tol) {
    tilt_mod[["coefficients"]] <- rep(0, n_coefs)
  }
  # output cleaned up fluctuation model
  return(tilt_mod)
}

###############################################################################

#' Mapping for Incremental Propensity Score Interventions
#'
#' @param g_est A \code{numeric} vector of estimated propensity scores, ie, the
#'  conditional probability of receiving treatment (A = 1) given covariates W.
#' @param delta A \code{numeric} value indicating the magnitude of shift in the
#'  treatment that defines the counterfactual quantity of interest. In the case
#'  of binary treatments, this defines the incremental propensity score shift,
#'  which acts as a multiplier of the odds with which an observational unit
#'  receives treatment (EH Kennedy, 2018; <doi:10.1080/01621459.2017.1422737>).
#'  When the treatment is quantitative, the defines the degree of shift for a
#'  modified treatment policy but non-binary treatments are not  yet supported.
ipsi_delta <- function(g_est, delta) {
  delta * g_est / (delta * g_est + 1 - g_est)
}
