utils::globalVariables(c("beta_est", "se_eif"))

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
#' @param delta_grid A \code{numeric} of values giving the varous degrees of
#'  shift in the intervention to be used in defining the causal quantity of
#'  interest. In the case of binary interventions, this takes the form of an
#'  incremental propensity score shift, acting as a multiplier of probabilities
#'  with which a given observational unit receives the intervention (EH Kennedy,
#'  2018, JASA; <doi:10.1080/01621459.2017.1422737>).
#' @param mult_type A \code{character} identifying the type of multipliers to be
#'  used in the multiplier bootstrap. Choices are limited to \code{"rademacher"}
#'  or \code{"gaussian"}, with the default being the former.
#' @param ci_level A \code{numeric} indicating the (1 - alpha) level of the
#'  simultaneous confidence band to be computed around the estimates of the
#'  direct effect. The error level of the test reported in the p-value returned
#'  is simply alpha, i.e., one less this quantity.
#' @param g_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, used in fitting a model for the propensity
#'  score, i.e., g = P(A | W).
#' @param e_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting a cleverly parameterized
#'  propensity score that includes the mediators, i.e., e = P(A | Z, W).
#' @param m_learners A \code{Stack} object, or other learner class (inheriting
#'  from \code{Lrnr_base}), containing a single or set of instantiated learners
#'  from the \code{sl3} package, to be used in fitting the outcome regression,
#'  i.e., m(A, Z, W).
#' @param phi_learners A \code{Stack} object, or other learner class
#'  (inheriting from \code{Lrnr_base}), containing a single or set of
#'  instantiated learners from the \code{sl3} package, to be used in fitting a
#'  reduced regression useful for computing the efficient one-step estimator,
#'  i.e., phi(W) = E[m(A = 1, Z, W) - m(A = 0, Z, W) | W).
#' @param cv_folds A \code{numeric} integer value specifying the number of folds
#'  to be created for cross-validation. Use of cross-validation / cross-fitting
#'  allows for entropy conditions on the AIPW estimator to be relaxed. Note: for
#'  compatibility with \code{origami::make_folds}, this value specified here
#'  must be greater than or equal to 2; the default is to create 10 folds.
#' @param n_mult A \code{numeric} scalar giving the number of repetitions of the
#'  multipliers to be used in computing the multiplier bootstrap.
#'
#' @importFrom stats var rbinom rnorm quantile
#' @importFrom data.table as.data.table setnames rbindlist
#' @importFrom dplyr "%>%" transmute between
#' @importFrom tibble as_tibble
#' @importFrom assertthat assert_that
#'
#' @export
#
test_de <- function(W,
                    A,
                    Z,
                    Y,
                    delta_grid = seq(from = 0.5, to = 5.0, by = 0.9),
                    mult_type = c("rademacher", "gaussian"),
                    ci_level = 0.95,
                    g_learners,
                    e_learners,
                    m_learners,
                    phi_learners,
                    cv_folds = 10,
                    n_mult = 10000) {
  # set default arguments
  mult_type <- match.arg(mult_type)

  # NOTE: procedure does _not_ support static interventions currently
  assertthat::assert_that(all(delta_grid > 0) && all(delta_grid < Inf))

  # construct input data structure
  data_in <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- paste("W", seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste("Z", seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )
  data.table::setnames(data_in, c("Y", z_names, "A", w_names))

  # compute E[Y] for half of direct effect and size of observed data
  EY <- mean(Y)
  n_obs <- length(Y)

  # generate multipliers
  if (mult_type == "rademacher") {
    mult_boot <- lapply(seq_len(n_mult), function(iter) {
      mults <- (2 * stats::rbinom(n_obs, 1, 0.5)) - 1
      return(mults)
    })
  } else if (mult_type == "gaussian") {
    mult_boot <- lapply(seq_len(n_mult), function(iter) {
      mults <- stats::rnorm(n_obs)
      return(mults)
    })
  }
  mult_boot_mat <- do.call(cbind, mult_boot)

  # estimate via one-step for other half of direct effect
  theta_est <- est_onestep(
    data = data_in,
    delta = delta_grid,
    g_learners = g_learners,
    e_learners = e_learners,
    m_learners = m_learners,
    phi_learners = phi_learners,
    w_names = w_names,
    z_names = z_names,
    cv_folds = cv_folds
  )

  # perform multiplier bootstrap to construct ingredients for simultaneous CI
  rho_est <- apply(mult_boot_mat, 2, function(mults) {
    # NOTE: there ought to be a better way than looping over the grid of deltas
    m_process_out <- lapply(seq_along(delta_grid), function(iter) {
      # compute estimate of the direct effect
      beta_est <- EY - theta_est[[iter]]$theta

      # compute corresponding uncentered influence function
      s_eif_est <- Y - (theta_est[[iter]]$eif + theta_est[[iter]]$theta)

      # difference in estimated influence function and direct effect
      eif_de_diff <- s_eif_est - beta_est

      # estimated variance for current shift
      se_eif <- sqrt(stats::var(eif_de_diff) / n_obs)

      # compute process M(delta) using multipliers
      m_process <- mean((mults * eif_de_diff) / se_eif)

      # construct output
      return(list(beta_est = beta_est, se_eif = se_eif, m_est = m_process))
    })
    m_process_out <- data.table::rbindlist(m_process_out)
    sup_m_over_delta <- max(abs(m_process_out$m_est))
    return(list(est = m_process_out[, -3], sup_m = sup_m_over_delta))
  })
  # NOTE: this is _very_ inefficient as we store the parameter and standard
  #       error estimates for each of the bootstrap iterations even though they
  #       are the same, i.e., we save n_mult data frames when we need just 1
  param_est <- lapply(rho_est, `[[`, "est")[[1]]

  # construct simultaneous CI using quantiles from supremum of M(delta)
  sup_m_process <- do.call(c, lapply(rho_est, `[[`, "sup_m"))
  c_alpha <- unname(stats::quantile(sup_m_process, ci_level))
  est_with_cis <- param_est %>%
    dplyr::transmute(
      lwr_ci = beta_est - c_alpha * se_eif,
      beta_est = I(beta_est),
      upr_ci = beta_est + c_alpha * se_eif,
      se_eif = I(se_eif),
      test_stat = beta_est / se_eif,
      delta = delta_grid
    ) %>%
    tibble::as_tibble()

  # for p-value, evaluate 1-rho(t) at observed supremum test statistic
  pval_rho <- 1 - mean(sup_m_process < max(abs(est_with_cis$test_stat)))

  # output
  out <- list(est_de = est_with_cis, pval_de = pval_rho)
  return(out)
}
