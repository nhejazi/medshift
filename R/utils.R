#' Bounding to numerical precision
#'
#' Bounds extreme values to a specified tolerance level, for use with sensitive
#' quantities that must be transformed, e.g., via \code{\link[stats]{qlogis}}.
#'
#' @param vals A \code{numeric} vector of values in the unit interval.
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
#' @param vals A \code{numeric} vector of values in the unit interval.
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
#'  unit interval.
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
#'  closed unit interval by use of \code{\link{scale_to_unit}}.
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
scale_from_unit <- function(vals_in_unit) {
  # check that input falls in the unit interval
  assertthat::assert_that(min(vals_in_unit) >= 0)
  assertthat::assert_that(max(vals_in_unit) <= 1)

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

#' Incremental Propensity Score Intervention
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
