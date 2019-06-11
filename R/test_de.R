#' Hypothesis test of direct effect with mediated stochastic interventions
#' using the multiplier bootstrap
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
#' @param n_mult A \code{numeric} scalar giving the number of repetitions of the
#'  multipliers (Rademacher or Gaussian) to be used in computing the multiplier
#'  bootstrap.
#' @param mult_type A \code{character} identifying the type of multipliers to be
#'  used in the multiplier bootstrap. Choices are limited to \code{"rademacher"}
#'  or \code{"gaussian"}, with the default being the former.
#' @param ... Other arguments to be passed to \code{\link{medshift}}.
#' @param estimator The desired estimator of the natural direct effect to be
#'  computed. Currently, choices are limited to a substitution estimator, a
#'  re-weighted estimator, and an efficient one-step estimator. The interested
#'  user should consider consulting Díaz & Hejazi (2019+) for a comparative
#'  investigation of each of these estimators.
#' @param estimator_args A \code{list} of extra arguments to be passed (via
#'  \code{...}) to the function call for the specified estimator. The default
#'  is so chosen as to allow the number of folds used in computing the AIPW
#'  estimator to be easily tweaked. Refer to the documentation for functions
#'  \code{\link{est_onestep}}, \code{\link{est_ipw}}, and
#'  \code{\link{est_substitution}} for details on what other arguments may be
#'  specified through this mechanism. For the option \code{"tmle"}, there is
#'  heavy reliance on the architecture provided by the \code{tmle3} package.
#'
#' @importFrom stats var rbinom rnorm
#'
#' @export
#
test_de <- function(W,
                    A,
                    Z,
                    Y,
                    delta_grid,
                    n_mult = 10000,
                    mult_type = c("rademacher", "gaussian"),
                    ...,
                    estimator = c("onestep", "tmle"),
                    estimator_args = list(
                       cv_folds = 10,
                       max_iter = 1e4,
                       step_size = 1e-6
                   )) {
  # set default arguments
  mult_type <- match.arg(mult_type)
  estimator <- match.arg(estimator)

  # compute E[Y] for half of direct effect
  EY <- mean(Y)

  # use `medshift` wrapper function for other half of direct effect
  de_est <- lapply(delta_grid, function(delta) {
    theta_est <- medshift(W = W, A = A, Z = Z, Y = Y,
                          delta = delta, ...,
                          estimator = estimator,
                          estimator_args = estimator_args)

    # compute estimate of the direct effect
    beta_est <- EY - theta_est$theta

    # compute corresponding un-centered influence function
    eif_est <- Y - (theta_est$eif + theta_est$theta)

    # difference in estimated influence function and direct effect
    eif_de_diff <- eif_est - beta_est

    # estimated variance for current shift
    se_eif <- sqrt(stats::var(eif_de_diff) / length(Y))

    # output: need influence function difference and estimated standard error
    out <- list(eif_de_diff = eif_de_diff, se_eif = se_eif)
    return(out)
  })

  # generate multipliers
  if (mult_type == "rademacher") {
    rade_mult <- stats::rbinom(n_mult, 1, 0.5)
    rade_mult[rade_mult == 0] <- -1
  } else if (mult_type == "gaussian") {
    norm_mult <- stats::rnorm(n_mult)
  }

  # ...

}
