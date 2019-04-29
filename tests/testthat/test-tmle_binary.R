context("TML estimator for incremental propensity score interventions")

library(data.table)
library(stringr)
library(future)
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

# set up TMLE components: NPSEM, likelihood
npsem <- list(
  define_node("W", c(
    "W_1", "W_2", "W_3"
  )),
  define_node("A", c("A"), c("W")),
  define_node("Z", c("Z_1", "Z_2", "Z_3"), c("A", "W")),
  define_node("Y", c("Y"), c("Z", "A", "W"))
)

factor_list <- list(
  define_lf(LF_emp, "W"),
  define_lf(LF_fit, "A", hal_binary_lrnr),
  define_lf(LF_fit, "Y", hal_contin_lrnr, type = "mean")
)

# create TMLE task
tmle_task <- tmle3_Task$new(data, npsem = npsem)
likelihood_def <- Likelihood$new(factor_list)
likelihood_init <- likelihood_def$train(tmle_task)
likelihood_init$get_likelihoods(tmle_task)

# NEXT, need targeted_likelihood constructor
updater <- tmle3_Update$new(cvtmle = FALSE)
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

# compute a tmle3 "by hand"
tsm <- define_param(Param_TSM, likelihood, intervention)
updater$tmle_params <- tsm

if (FALSE) {
  # define data (from tmle3_Spec base class)
  tmle_task <- tmle_spec$make_tmle_task(data, node_list)

  # define likelihood (from tmle3_Spec base class)
  likelihood_init <- tmle_spec$make_initial_likelihood(tmle_task, learner_list)

  # define update method (fluctuation submodel and loss function)
  updater <- tmle_spec$make_updater()
  likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

  # invoke params specified in spec
  tmle_params <- tmle_spec$make_params(tmle_task, likelihood_targeted)
  updater$tmle_params <- tmle_params

  # fit TML estimator update
  tmle_fit <- fit_tmle3(tmle_task, likelihood_targeted, tmle_params, updater)

  # extract results from tmle3_Fit object
  tmle_fit
}
