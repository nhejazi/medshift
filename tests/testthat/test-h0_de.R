context("Uniform test of no direct effect")

library(data.table)
library(stringr)
library(hal9001)
library(sl3)
set.seed(7128816)
alpha_level <- 0.05

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

  # subject-level IDs
  ids <- seq_along(Y)

  # full data structure
  data <- as.data.table(cbind(Y, Z, A, W, ids))
  setnames(data, c(
    "Y", paste("Z", 1:3, sep = "_"), "A",
    paste("W", seq_len(dim(W)[2]), sep = "_"),
    "ids"
  ))
  return(data)
}

# get data and column names for sl3 tasks (for convenience)
data <- make_simulated_data()
z_names <- colnames(data)[str_detect(colnames(data), "Z")]
w_names <- colnames(data)[str_detect(colnames(data), "W")]

###############################################################################
# test of hypothesis testing procedure
###############################################################################
ranger_lrnr <- Lrnr_ranger$new(
  min.node.size = 5000, oob.error = FALSE, num.trees = 1000,
  sample.fraction = 0.7
)

de_test <- test_de(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  ids = data$ids, delta_grid = seq(from = 0.5, to = 5, by = 0.9),
  mult_type = "rademacher",
  ci_level = (1 - alpha_level),
  g_learners = ranger_lrnr,
  e_learners = ranger_lrnr,
  m_learners = ranger_lrnr,
  phi_learners = ranger_lrnr,
  cv_folds = 5
)

test_that("Uniform test of no direct effect rejects H0", {
  expect_gte(de_test$pval_de, alpha_level)
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
