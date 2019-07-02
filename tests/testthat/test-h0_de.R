context("Uniform test of no direct effect")

library(data.table)
library(stringr)
library(hal9001)
library(sl3)
set.seed(7128816)
alpha_level <- 0.05

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
  Y <- Z_1 + Z_2 - Z_3 - 0.1 * rowSums(W)^2 + rnorm(n_obs, mean = 0, sd = 0.5)

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
# test of hypothesis testing procedure
################################################################################
de_test <- test_de(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta_grid = seq(from = 0.5, to = 5, by = 0.9),
  mult_type = "rademacher",
  ci_level = (1 - alpha_level),
  g_learners = hal_binary_lrnr,
  e_learners = hal_binary_lrnr,
  m_learners = hal_contin_lrnr,
  phi_learners = hal_contin_lrnr,
  cv_folds = 5
)

test_that("Uniform test of no direct effect rejects H0 with fixed p-value", {
  expect_equal(de_test$pval_de, 0.822, tol = 1e-3)
})

test_that("Simultaneous confidence band uniformly covers zero under H0", {
  covers_zero <- apply(
    de_test$est_de[, c("lwr_ci", "upr_ci")], 1,
    function(ci) {
      check_zero <- dplyr::between(0, ci[1], ci[2])
      return(check_zero)
    }
  )
  expect_true(all(covers_zero))
})
