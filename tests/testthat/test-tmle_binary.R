context("TML estimator for incremental propensity score interventions")

library(data.table)
library(stringr)
library(hal9001)
library(sl3)
library(tmle3)
set.seed(7128816)
delta_ipsi <- 0.5

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

# create node list and learner list
node_list <- list(
  W = c("W_1", "W_2", "W_3"),
  A = "A",
  Z = c("Z_1", "Z_2", "Z_3"),
  Y = "Y"
)
learner_list <- list(
  Y = cv_hal_contin_lrnr,
  A = cv_hal_binary_lrnr
)

## instantiate tmle3 spec for stochastic mediation
tmle_spec <- tmle_medshift(
  delta = delta_ipsi,
  e_learners = cv_hal_binary_lrnr,
  phi_learners = cv_hal_contin_lrnr
)

## define data (from tmle3_Spec base class)
tmle_task <- tmle_spec$make_tmle_task(data, node_list)

## define likelihood (from tmle3_Spec base class)
likelihood_init <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

## define update method (submodel and loss function)
updater <- tmle_spec$make_updater()
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

## define param
tmle_params <- tmle_spec$make_params(tmle_task, likelihood_targeted)
updater$tmle_params <- tmle_params

## fit tmle update
tmle_fit <- fit_tmle3(tmle_task, likelihood_targeted, tmle_params, updater)


## one-line call with faster with tmle3 wrapper
set.seed(71281)
tmle_fit <- tmle3(tmle_spec, data, node_list, learner_list)
tmle_fit


# fit one-step estimator
set.seed(71281)
os_fit <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta_ipsi,
  g_learners = hal_binary_lrnr,
  e_learners = hal_binary_lrnr,
  m_learners = hal_contin_lrnr,
  phi_learners = hal_contin_lrnr,
  estimator = "onestep",
)
summary(os_fit)

# fit substitution estimator
set.seed(71281)
sub_fit <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta_ipsi,
  g_learners = hal_binary_lrnr,
  e_learners = hal_binary_lrnr,
  m_learners = hal_contin_lrnr,
  phi_learners = hal_contin_lrnr,
  estimator = "substitution",
)
summary(sub_fit)

# test --- what exactly?
test_that("TML estimate...", {
})

################################################################################
if (FALSE) {
  get_sim_truth <- function(n_obs = 1e7, # number of observations
                              n_w = 3, # number of baseline covariates
                              delta = 0.5) { # value of shift parameter

    # compute large data set for true values
    data <- make_simulated_data(
      n_obs = n_obs,
      n_w = n_w,
      delta = delta
    )
    w_names <- str_subset(colnames(data), "W")
    z_names <- str_subset(colnames(data), "Z")
    W <- data[, ..w_names]
    Z <- data[, ..z_names]
    Y <- data$Y

    # compute TRUE G under counterfactual regimes
    g_Ais1 <- rowSums(W) / 4 + 0.1
    g_Ais0 <- 1 - g_Ais1

    # compute TRUE SHIFTED G under counterfactual regimes
    g_shifted_Ais1 <- (delta * g_Ais1) / (delta * g_Ais1 + g_Ais0)
    g_shifted_Ais0 <- 1 - g_shifted_Ais1

    # compute TRUE M under counterfactual regimes
    m_Ais1 <- Z$Z_1 + Z$Z_2 - Z$Z_3 + 1 - 0.1 * rowSums(W)^2
    m_Ais0 <- Z$Z_1 + Z$Z_2 - Z$Z_3 + 0 - 0.1 * rowSums(W)^2

    # output: true values of nuisance parameters
    return(list(
      g_obs_true = data.table(
        A1 = g_Ais1,
        A0 = g_Ais0
      ),
      g_shifted_true = data.table(
        A1 = g_shifted_Ais1,
        A0 = g_shifted_Ais0
      ),
      m_true = data.table(
        A1 = m_Ais1,
        A0 = m_Ais0
      ),
      EY_true = mean(Y)
    ))
  }

  # simulate data and extract components for computing true parameter value
  sim_truth <- get_sim_truth()
  m_A1 <- sim_truth$m_true$A1
  m_A0 <- sim_truth$m_true$A0
  g_shifted_A1 <- sim_truth$g_shifted_true$A1
  g_shifted_A0 <- sim_truth$g_shifted_true$A0
  EY <- sim_truth$EY_true

  # compute true parameter value based on the substitution estimator
  true_param <- mean(m_A1 * g_shifted_A1) + mean(m_A0 * g_shifted_A0)
}
