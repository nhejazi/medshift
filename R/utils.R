#' Confidence Intervals for Stochastic Mediation Parameters
#'
#' Compute confidence intervals for estimates produced by \code{tmle_medshift}
#'
#' @param object An object of class \code{medshift}, as produced by invoking
#'  the function \code{tmle_medshift}, for which a confidence interval is to be
#'  computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the level of the confidence interval
#'  to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats qnorm
#'
#' @method confint medshift
#'
#' @export
#
confint.medshift <- function(object,
                             parm = seq_len(object$psi),
                             level = 0.95,
                             ...) {

  # first, let's get Z_(1 - alpha)
  norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - level) / 2))

  # compute the EIF variance multiplier for the CI
  # NOTE: the variance value is already scaled by length of observations
  sd_eif <- sqrt(object$var)

  # compute the interval around the point estimate
  ci_psi <- norm_bounds * sd_eif + object$psi

  # set up output CI object
  ci_out <- c(ci_psi[1], object$psi, ci_psi[2])
  names(ci_out) <- c("lwr_ci", "est", "upr_ci")
  return(ci_out)
}

################################################################################

#' Summary for Stochastic Mediation Parameter Objects
#'
#' Print a convenient summary for objects computed using \code{tmle_medshift}.
#'
#' @param object An object of class \code{medshift}, as produced by invoking
#'  the function \code{tmle_medshift}, for which a confidence interval is to be
#'  computed.
#' @param ... Other arguments. Not currently used.
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @importFrom stats confint
#'
#' @method summary medshift
#'
#' @export
#
summary.medshift <- function(object,
                             ...,
                             ci_level = 0.95) {

  # compute confidence interval using the pre-defined method
  ci <- stats::confint(object, level = ci_level)

  # only print useful info about the mean of the efficient influence function
  eif_mean <- format(mean(object$eif), scientific = TRUE)

  # create output table from input object and confidence interval results
  out <- c(round(c(ci, object$var), digits = 6), eif_mean, object$n_iter)
  names(out) <- c(
    "lwr_ci", "param_est", "upr_ci", "param_var"
  )
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
#' @export
#'
#' @method print medshift
#'
#
print.medshift <- function(x, ...) {
  print(x[c("psi", "var", "msg", "n_iter")])
}

################################################################################

#' Bounding Numerical Precision
#'
#' Bounds extreme values to numerical (machine) precision, for use with
#' sensitive quantities like estimated propensity scores.
#'
#' @param vals A \code{numeric} vector of values in the interval [0, 1].
#'
#' @importFrom assertthat assert_that
#'
#' @keywords internal
#
bound_precision <- function(vals) {
  assertthat::assert_that(!(max(vals) >= 1 | min(vals) <= 0))
  vals[vals == 0] <- .Machine$double.neg.eps
  vals[vals == 1] <- 1 - .Machine$double.neg.eps
  return(vals)
}

################################################################################

#' Bounding Propensity Scores
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
#
bound_propensity <- function(vals, bounds = c(0.01, 0.99)) {
  assertthat::assert_that(!(max(vals) >= 1 | min(vals) <= 0))
  vals[vals < bounds[1]] <- bounds[1]
  vals[vals > bounds[2]] <- bounds[2]
  return(vals)
}