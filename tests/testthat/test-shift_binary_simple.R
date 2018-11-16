context("Estimators work for simple incremental propensity score interventions")

library(data.table)
library(stringr)
library(future)
library(sl3)
set.seed(429153)

################################################################################
# setup learners for the nuisance parameters
################################################################################

# instantiate some learners
mean_lrnr <- Lrnr_mean$new()
fglm_contin_lrnr <- Lrnr_glm_fast$new()
fglm_binary_lrnr <- Lrnr_glm_fast$new(family = binomial())
rf_contin_lrnr <- Lrnr_ranger$new()
rf_binary_lrnr <- Lrnr_ranger$new(family = "binomial")
xgb_contin_lrnr <- Lrnr_xgboost$new(nrounds = 10)
xgb_binary_lrnr <- Lrnr_xgboost$new(nrounds = 10, objective = "reg:logistic")
hal_contin_lrnr <- Lrnr_hal9001$new(fit_type = "glmnet", n_folds = 5)
hal_binary_lrnr <- Lrnr_hal9001$new(fit_type = "glmnet", n_folds = 5,
                                    family = "binomial")

# learner stack for the propensity score
sl_lrn_g <- Lrnr_sl$new(
  learners = list(fglm_binary_lrnr, rf_binary_lrnr, hal_binary_lrnr),
  metalearner = Lrnr_nnls$new()
)

# learner stack for the clever conditional regression e
sl_lrn_e <- Lrnr_sl$new(
  learners = list(fglm_binary_lrnr, rf_binary_lrnr, hal_binary_lrnr),
  metalearner = Lrnr_nnls$new()
)

# learner stack for the outcome regression m
sl_lrn_m <- Lrnr_sl$new(
  learners = list(mean_lrnr, fglm_contin_lrnr, rf_contin_lrnr, hal_contin_lrnr),
  metalearner = Lrnr_nnls$new()
)

# learner stack for reduced-dimension regression phi
sl_lrn_phi <- Lrnr_sl$new(
  learners = list(mean_lrnr, fglm_contin_lrnr, rf_contin_lrnr, hal_contin_lrnr),
  metalearner = Lrnr_nnls$new()
)


################################################################################
# setup data and simulate to test with estimators
################################################################################

# simulate simple data for simple simulation
n_obs <- 1000  # number of observations
n_w <- 3  # number of baseline covariates
p_w <- 0.5  # probability of a success ("1") in the baseline variables
tx_mult <- 2  # multiplier for the effect of W = 1 on the treatment
delta_shift <- 0.5  # posited value of the shift parameter

# baseline covariate -- simple, binary
W <- as.matrix(replicate(n_w, rbinom(n_obs, 1, prob = p_w)))

# create treatment based on baseline W
A <- as.numeric(rbinom(n_obs, 1, prob = (rowSums(W)/5 + 0.1)))

# mediators to affect the outcome
## 1st mediator (binary)
z1_prob <-(A + rowSums(W) + runif(n_obs, 0, 0.1)) / max(A + rowSums(W) + 0.1)
Z_1 <- rbinom(n_obs, 1, prob = z1_prob)
## 2nd mediator (binary)
z2_form <- (A - 1)^3 + W[, 2] - W[, 3] + runif(n_obs, 0, 0.2)
z2_prob <- abs(z2_form) / max(abs(z2_form))
z2_prob <- z2_prob + runif(n_obs, 0, 1 - max(z2_prob))
Z_2 <- rbinom(n_obs, 1, prob = z2_prob)
## 3rd mediator (binary)
z3_form <- (A - 1) + rnorm(n_obs)
z3_prob <- abs(z3_form) / max(abs(z3_form))
Z_3 <- rbinom(n_obs, 1, prob = z3_prob)
## build matrix of mediators
Z <- cbind(Z_1, Z_2, Z_3)

# create outcome as a linear function of A, W + white noise
Y <- Z_1 + Z_2 - Z_3 + A - rowSums(W) + rnorm(n_obs, mean = 0, sd = 1)

# full data structure
data <- as.data.table(cbind(Y, Z, A, W))
setnames(data, c("Y", paste("Z", 1:3, sep = "_"), "A",
                paste("W", seq_len(dim(W)[2]), sep = "_")))
head(data)

# column names for sl3 tasks (for convenience)
z_names <- colnames(data)[str_detect(colnames(data), "Z")]
w_names <- colnames(data)[str_detect(colnames(data), "W")]



################################################################################
# test different estimators
################################################################################
theta_sub <- medshift(W = W, A = A, Z = Z, Y = Y,
                      shift_value = delta_shift,
                      #g_lrnrs = sl_lrn_g,
                      #e_lrnrs = sl_lrn_g,
                      #m_lrnrs = sl_lrn_m,
                      #phi_lrnrs = sl_lrn_phi,
                      estimator = "substitution")
theta_sub

theta_re <- medshift(W = W, A = A, Z = Z, Y = Y,
                     shift_value = delta_shift,
                     #g_lrnrs = sl_lrn_g,
                     #e_lrnrs = sl_lrn_g,
                     #m_lrnrs = sl_lrn_m,
                     #phi_lrnrs = sl_lrn_phi,
                     estimator = "reweighted")
theta_re

theta_eff <- medshift(W = W, A = A, Z = Z, Y = Y,
                      shift_value = delta_shift,
                      #g_lrnrs = sl_lrn_g,
                      #e_lrnrs = sl_lrn_g,
                      #m_lrnrs = sl_lrn_m,
                      #phi_lrnrs = sl_lrn_phi,
                      estimator = "efficient")
theta_eff

