context("Estimators agree for incremental propensity score interventions")

library(data.table)
library(stringr)
library(hal9001)
library(sl3)
library(tmle3)
set.seed(7128816)
delta <- 0.5

################################################################################
# setup learners for the nuisance parameters
################################################################################

# instantiate some learners
mean_lrnr <- Lrnr_mean$new()
fglm_contin_lrnr <- Lrnr_glm_fast$new()
fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
hal_contin_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5
)
hal_binary_lrnr <- Lrnr_hal9001$new(
  fit_type = "glmnet", n_folds = 5,
  family = "binomial"
)
cv_hal_contin_lrnr <- Lrnr_cv$new(hal_contin_lrnr, full_fit = TRUE)
cv_hal_binary_lrnr <- Lrnr_cv$new(hal_binary_lrnr, full_fit = TRUE)

################################################################################
# setup data and simulate to test with estimators
################################################################################
make_simulated_data <- function(n_obs = 1000, # no. observations
                                n_w = 3, # no. baseline covariates
                                delta = 0.5) { # shift parameter value

  # baseline covariate -- simple, binary
  W_1 <- rbinom(n_obs, 1, prob = 0.50)
  W_2 <- rbinom(n_obs, 1, prob = 0.65)
  W_3 <- rbinom(n_obs, 1, prob = 0.35)
  W <- cbind(W_1, W_2, W_3)

  # create treatment based on baseline W
  A <- as.numeric(rbinom(n_obs, 1, prob = (rowSums(W) / 4 + 0.1)))

  # mediators to affect the outcome
  ## 1st mediator (binary)
  z1_prob <- 1 - plogis((A^2 + W[, 1]) / (A + W[, 1]^3 + 0.5))
  Z_1 <- rbinom(n_obs, 1, prob = z1_prob)
  ## 2nd mediator (binary)
  z2_prob <- plogis((A - 1)^3 + W[, 2] / (W[, 3] + 3))
  Z_2 <- rbinom(n_obs, 1, prob = z2_prob)
  ## 3rd mediator (binary)
  z3_prob <- plogis((A - 1)^2 + 2 * W[, 1]^3 - 1 / (2 * W[, 1] + 0.5))
  Z_3 <- rbinom(n_obs, 1, prob = z3_prob)
  ## build matrix of mediators
  Z <- cbind(Z_1, Z_2, Z_3)

  # create outcome as a linear function of A, W + white noise
  Y <- Z_1 + Z_2 - Z_3 + A - 0.1 * rowSums(W)^2 +
    rnorm(n_obs, mean = 0, sd = 0.5)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W))
  setnames(data, c(
    "Y", paste("Z", 1:3, sep = "_"), "A",
    paste("W", seq_len(dim(W)[2]), sep = "_")
  ))
  return(data)
}

# get data and column names for sl3 tasks (for convenience)
data <- make_simulated_data()
z_names <- colnames(data)[str_detect(colnames(data), "Z")]
w_names <- colnames(data)[str_detect(colnames(data), "W")]


################################################################################
# test different estimators
################################################################################
set.seed(7128816)
theta_sub <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta,
  g_lrnrs = hal_binary_lrnr,
  e_lrnrs = hal_binary_lrnr,
  m_lrnrs = hal_contin_lrnr,
  phi_lrnrs = hal_contin_lrnr,
  estimator = "substitution"
)
theta_sub

set.seed(7128816)
theta_re <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta,
  g_lrnrs = hal_binary_lrnr,
  e_lrnrs = hal_binary_lrnr,
  m_lrnrs = hal_contin_lrnr,
  phi_lrnrs = hal_contin_lrnr,
  estimator = "reweighted"
)
theta_re

set.seed(7128816)
theta_os <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta,
  g_lrnrs = hal_binary_lrnr,
  e_lrnrs = hal_binary_lrnr,
  m_lrnrs = hal_contin_lrnr,
  phi_lrnrs = hal_contin_lrnr,
  estimator = "onestep",
  estimator_args = list(cv_folds = 10)
)
theta_os

set.seed(7128816)
theta_tmle <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta,
  g_lrnrs = cv_hal_binary_lrnr,
  e_lrnrs = cv_hal_binary_lrnr,
  m_lrnrs = cv_hal_contin_lrnr,
  phi_lrnrs = cv_hal_contin_lrnr,
  estimator = "tmle",
  estimator_args = list(
    max_iter = 100,
    step_size = 1e-5
  )
)
theta_tmle

test_that("Substitution and re-weighted estimator agree", {
  expect_equal(theta_sub$theta, theta_re$theta, tol = 1e-2)
})

test_that("Substitution and one-step estimator agree", {
  expect_equal(theta_sub$theta, theta_os$theta, tol = 1e-2)
})

test_that("Re-weighted and one-step estimator agree", {
  expect_equal(theta_re$theta, theta_os$theta, tol = 1e-2)
})

test_that("Substitution and TML estimator agree", {
  expect_equal(theta_sub$theta, theta_tmle$summary$psi, tol = 1e-2)
})

test_that("Re-weighted and TML estimator agree", {
  expect_equal(theta_re$theta, theta_tmle$summary$psi, tol = 1e-2)
})

test_that("One-step and TML estimator agree", {
  expect_equal(theta_re$theta, theta_tmle$summary$psi, tol = 1e-2)
})
