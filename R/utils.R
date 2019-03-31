#' Confidence Intervals for Stochastic Mediation Parameter Objects
#'
#' Compute confidence intervals for objects of class \code{medshift}, which
#' contain estimates produced by \code{medshift}.
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
#' @importFrom assertthat assert_that
#'
#' @method confint medshift
#'
#' @export
#
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
  # inference is currently limited to the one-step efficient estimator
  if (object$type %in% c("onestep", "tmle")) {
    # compute confidence interval using the pre-defined method
    ci <- stats::confint(object, level = ci_level)

    # only print useful info about the mean of the efficient influence function
    eif_mean <- formatC(mean(object$eif), digits = 4, format = "e")

    # create output table from input object and confidence interval results
    out <- c(round(c(ci, object$var), digits = 4), eif_mean, object$type)
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

#' Print Method for Stochastic Mediation Parameter Objects
#'
#' The \code{print} method for objects of class \code{medshift}.
#'
#' @param x An object of class \code{medshift}.
#' @param ... Other options (not currently used).
#'
#' @method print medshift
#'
#' @export
#
print.medshift <- function(x, ...) {
  # inference is currently limited to the one-step efficient estimator
  # TODO: allow use for TML estimators once impelemented
  if (x$type == "one-step efficient") {
    print(x[c("theta", "var", "type")])
  } else {
    print(x[c("theta", "type")])
  }
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

################################################################################

#' Scale Values to the Unit Interval [0, 1]
#'
#' @param vals A \code{numeric} vector of values to be scaled into the interval
#'  [0, 1].
#'
#' @keywords internal
#
scale_to_unit <- function(vals) {
  vals_scaled <- (vals - min(vals)) / (max(vals) - min(vals))
  return(vals_scaled)
}

################################################################################

#' Scale Values to the Original Scale
#'
#' @param scaled_vals A \code{numeric} vector of values scaled to lie in the
#'  interval [0, 1] by use of \code{\link{scale_to_unit}}.
#' @param max_orig The maximum of the values on the original scale.
#' @param min_orig The minimum of the values on the original scale.
#'
#' @keywords internal
#
scale_to_original <- function(scaled_vals, max_orig, min_orig) {
  vals_orig <- (scaled_vals * (max_orig - min_orig)) + min_orig
  return(vals_orig)
}

################################################################################

#' Numerical Integration for Weighted Treatment Mechanism
#'
#' In the case of modified treatment policies, it is necessary to numerically
#' evaluate an integral over the domain of the treatment mechanism. This is a
#' simple procedure to numerically compute such an integral based on Monte
#' Carlo importance sampling from a uniform distribution.
#'
#' @param data A \code{data.table} containing the observed data, with columns
#'  in the order specified by the NPSEM (Y, Z, A, W), with column names set
#'  appropriately based on the original input data. Such a structure is merely
#'  a convenience utility to passing data around to the various core estimation
#'  routines and is automatically generated as part of a call to the user-facing
#'  wrapper function \code{\link{medshift}}.
#' @param delta A \code{numeric} value indicating the degree of shift in the
#'  intervention to be used in defining the causal quantity of interest. In the
#'  case of binary interventions, this takes the form of an incremental
#'  propensity score shift, acting as a multiplier of the probability with which
#'  a given observational unit receives the intervention (EH Kennedy, 2018,
#'  JASA; <doi:10.1080/01621459.2017.1422737>). In the case of continuous
#'  interventions, this is a modified treatment policy that shifts the observed
#'  treatment by the given value.
#' @param mc_draws A \code{numeric} vector corresponding to draws from a uniform
#'  distribution based on the range of the observed treatment mechanism. The
#'  numerical integral is computed by using this grid as values of intervention
#'  mechanism.
#' @param dens_mech Estimated conditional density corresponding to the natural
#'  or shifted value of the treatment for modified treatment policies based on
#'  the observed values of the treatment. This is an object with inherited class
#'  \code{Lrnr_base}, automatically produced by \code{\link{fit_g_mech}}.
#' @param wts_mech Estimated weighting mechanism learned from the observed data
#'  For the substitution estimator, this is the outcome mechanism, while, for
#'  the intervention score, this is a nuisance parameter defined by the type
#'  of shift intervention to be evaluated. This is an object with inherited
#'  class \code{Lrnr_base}.
#'
#' @importFrom data.table setkey setorder
#' @importFrom stringr str_subset
#' @importFrom sl3 sl3_Task
#'
#' @keywords internal
#
mc_integrate_dens <- function(data, delta, mc_draws, dens_mech, wts_mech) {
  # create expanded data set with multiple records for MC draws
  data[, id := 1:.N]
  data.table::setkey(data, id)
  mc_records_data <- data[rep(1:.N, length(mc_draws)), ]
  data.table::setorder(mc_records_data, id)
  mc_records_data[, A := rep(mc_draws, nrow(data))]

  # get names of baseline covariates and mediators
  w_names <- stringr::str_subset(colnames(mc_records_data), "W")
  z_names <- stringr::str_subset(colnames(mc_records_data), "Z")

  # numerical integration over the domain of A via Monte Carlo
  min_mc_scaling <- min(mc_draws)
  max_mc_scaling <- max(mc_draws)
  range_mc_scaling <- max_mc_scaling - min_mc_scaling

  # make outcome task for Monte Carlo integration (never shifted)
  m_mc_task <- sl3_Task$new(
    data = mc_records_data,
    covariates = c(w_names, "A", z_names),
    outcome_type = "continuous",
    outcome = "Y"
  )

  # compute Monte Carlo integral by predicting from task on expanded data
  m_mc_pred <- wts_mech$predict(m_mc_task)

  # estimate shifted intervention distribution with records data
  mc_records_data[, A := mtp_shift(A = A, delta = delta)]
  g_mc_task <- sl3_Task$new(
    data = mc_records_data,
    covariates = w_names,
    outcome_type = "continuous",
    outcome = "A"
  )

  # get predictions from natural propensity score model for Monte Carlo data
  g_mc_pred <- dens_mech$predict(g_mc_task)

  # compute weighted combination of terms for numerical integral
  # NOTE: sum over product of predicted values based on the observation ID,
  #       then re-scale based on the interval width
  mc_integ <- as.numeric(by(g_mc_pred * m_mc_pred, mc_records_data$id, sum)) *
    (range_mc_scaling / length(mc_draws))
  return(mc_integ)
}
