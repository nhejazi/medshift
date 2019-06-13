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
#' @param error_level A \code{numeric} [TO FILL IN]
#' @param mult_type A \code{character} identifying the type of multipliers to be
#'  used in the multiplier bootstrap. Choices are limited to \code{"rademacher"}
#'  or \code{"gaussian"}, with the default being the former.
#' @param ... Other arguments to be passed to \code{\link{medshift}}.
#'
#' @importFrom stats var rbinom rnorm
#'
#' @export
#
test_de <- function(W,
                    A,
                    Z,
                    Y,
                    delta_grid = seq(0.1, 0.9, 0.2),
                    n_mult = 10000,
                    mult_type = c("rademacher", "gaussian"),
                    error_level = 0.05,
                    ...) {
  # set default arguments
  mult_type <- match.arg(mult_type)
  estimator <- match.arg(estimator)

  # significance cutoffs
  sig_cutoffs <- c(error_level/2, 1 - error_level/2)

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste("Z", seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )
  data.table::setnames(data, c("Y", z_names, "A", w_names))

  # compute E[Y] for half of direct effect and size of observed data
  EY <- mean(Y)
  n_obs <- length(Y)

  # generate multipliers
  if (mult_type == "rademacher") {
    mult_boot <- lapply(seq_len(n_mult), function(iter) {
      mults <- stats::rbinom(n_obs, 1, 0.5)
      mults[mults == 0] <- -1
      return(mults)
    })
  } else if (mult_type == "gaussian") {
    mult_boot <- lapply(seq_len(n_mult), function(iter) {
      mults <- stats::rnorm(n_obs)
      return(mults)
    })
  }
  mult_boot_mat <- do.call(cbind, mult_boot)

  # use estimate one-step for other half of direct effect
  theta_est <- est_onestep(data = data,
                           delta = delta_grid,
                           g_lrnrs = hal_binary_lrnr,
                           e_lrnrs = hal_binary_lrnr,
                           m_lrnrs = hal_contin_lrnr,
                           phi_lrnrs = hal_contin_lrnr,
                           w_names = w_names,
                           z_names = z_names,
                           cv_folds = 5)

  # construct process M(delta) by using multiplier bootstrap
  rho_est <- apply(mult_boot_mat, 2, function(mults) {
    # NOTE: there ought to be a better way then looping over the grid of delta
    m_process_out <- lapply(seq_along(delta_grid), function(iter) {
      # compute estimate of the direct effect
      beta_est <- EY - theta_est[[iter]]$theta

      # compute corresponding un-centered influence function
      eif_est <- Y - (theta_est[[iter]]$eif + theta_est[[iter]]$theta)

      # difference in estimated influence function and direct effect
      eif_de_diff <- eif_est - beta_est

      # estimated variance for current shift
      se_eif <- sqrt(stats::var(eif_de_diff) / n_obs)

      # compute process M(delta) using multipliers
      m_process <- mean((mults * eif_de_diff) / se_eif)
      return(m_process)
    })
    m_process <- do.call(c, m_process_out)
    sup_m_process <- max(abs(m_process))
  })

  # need Pr(sup_{delta} M(delta) <= t | O_1,...,O_n)
  t_est <- quantile(rho_est, sig_levels)

}
