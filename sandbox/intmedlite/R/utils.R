#' Incremental Propensity Score Intervention
#'
#' @param g ...
#' @param delta ...
gdelta1 <- function(g, delta) {
  delta * g / (delta * g + 1 - g)
}

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

#' Bound Predicted Probability
#'
#' @param x ...
#' @param p ...
bound <- function(x, p = 0.005) {
  pmax(pmin(x, 1 - p), p)
}

#' Truncate Predictions
#'
#' @param x ...
#' @param p ...
truncate <- function(x, p = 0.01) {
  pmin(pmax(x, p), 1 - p)
}
