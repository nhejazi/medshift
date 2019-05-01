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

# create node list and learner list
node_list <- list(W = c("W_1", "W_2", "W_3"),
                  A = "A",
                  Z = c("Z_1", "Z_2", "Z_3"),
                  Y = "Y")
learner_list <- list(Y = hal_contin_lrnr,
                     A = hal_binary_lrnr)

# set up TMLE components: NPSEM, likelihood, TMLE task
npsem <- stochastic_mediation_npsem(node_list)
tmle_task <- tmle3_Task$new(data, npsem = npsem)
likelihood_init <- stochastic_mediation_likelihood(tmle_task, learner_list)
likelihood_init$get_likelihoods(tmle_task)

# NEXT, need targeted_likelihood constructor
updater <- tmle3_Update$new(cvtmle = FALSE)
likelihood_targeted <- Targeted_Likelihood$new(likelihood_init, updater)

# add derived likelihood factors to targeted likelihood object
lf_e <- tmle3::define_lf(tmle3::LF_derived, "E", hal_binary_lrnr,
                         likelihood_targeted, medshift::make_e_task)
lf_phi <- tmle3::define_lf(tmle3::LF_derived, "phi", hal_contin_lrnr,
                           likelihood_targeted, medshift::make_phi_task)
likelihood_targeted$add_factors(lf_e)
likelihood_targeted$add_factors(lf_phi)

# compute a tmle3 "by hand"
tmle_params <- define_param(Param_medshift, likelihood_targeted,
                            shift_param = delta)
updater$tmle_params <- tmle_params

# separately test param methods
tmle_params$clever_covariates(tmle_task)
theta_tmle <- tmle_params$estimates(tmle_task)

# fit TML estimator update
tmle_fit <- fit_tmle3(tmle_task, likelihood_targeted, tmle_params, updater)

# fit one-step estimator
os_fit <- medshift(
  W = data[, ..w_names], A = data$A, Z = data[, ..z_names], Y = data$Y,
  delta = delta,
  g_lrnrs = hal_binary_lrnr,
  e_lrnrs = hal_binary_lrnr,
  m_lrnrs = hal_contin_lrnr,
  phi_lrnrs = hal_contin_lrnr,
  estimator = "onestep",
)



# how Specs are used...
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

