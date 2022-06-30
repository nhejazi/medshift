###############################################################################
# fit estimators
###############################################################################
fit_estimators <- function(data, delta, cv_folds = 5) {
  # NOTE: specify HAL arguments based on Ivan's experience: "main change is
  #       that I am using 'mse' for the loss function (default is 'deviance'),
  #       I’m using 'use_min = FALSE' + I increased the value of 'nlambda'."
  ## NOTE: need to force `family = "gaussian"` to use `lambda.min.ratio`
  hal_lrnr <- Lrnr_hal9001$new(max_degree = NULL,
                               n_folds = 5,
                               fit_type = "glmnet",
                               use_min = FALSE,
                               type.measure = "mse",
                               standardize = FALSE,
                               family = "gaussian",
                               lambda.min.ratio = 1 / nrow(data),
                               nlambda = 1000,
                               yolo = FALSE)
  source(here::here("src", "sl_id.R"))
  myglmnet_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.myglmnet")
  myglm_lrnr <- make_learner(Lrnr_pkg_SuperLearner, "SL.myglm")
  mean_lrnr <- Lrnr_mean$new()

  # meta-learner to ensure predicted probabilities do not go outside [0,1]
  logistic_metalearner <- make_learner(Lrnr_solnp,
                                       metalearner_logistic_binomial,
                                       loss_loglik_binomial)

  # set nuisance regression learners based on ID's successful simulations
  sl <- Lrnr_sl$new(learners = list(myglmnet_lrnr, myglm_lrnr, hal_lrnr),
                    metalearner = Lrnr_nnls$new())
                    #metalearner = logistic_metalearner)

  # compute TML and one-step estimators
  ## 1) all nuisance functions correctly specified
  est_allc <- estimators(data = data,
                         delta = delta,
                         g_stack = sl,
                         e_stack = sl,
                         m_stack = sl,
                         b_stack = sl,
                         d_stack = sl,
                         cv_folds = cv_folds)

  ## 2) g misspecified
  est_misg <- estimators(data = data,
                         delta = delta,
                         g_stack = mean_lrnr,
                         e_stack = sl,
                         m_stack = sl,
                         b_stack = sl,
                         d_stack = sl,
                         cv_folds = cv_folds)

  ## 3) e misspecified
  est_mise <- estimators(data = data,
                         delta = delta,
                         g_stack = sl,
                         e_stack = mean_lrnr,
                         m_stack = sl,
                         b_stack = sl,
                         d_stack = sl,
                         cv_folds = cv_folds)

  ## 4) m misspecified
  est_mism <- estimators(data = data,
                         delta = delta,
                         g_stack = sl,
                         e_stack = sl,
                         m_stack = mean_lrnr,
                         b_stack = sl,
                         d_stack = sl,
                         cv_folds = cv_folds)

  ## 5) b misspecified
  est_misb <- estimators(data = data,
                         delta = delta,
                         g_stack = sl,
                         e_stack = sl,
                         m_stack = sl,
                         b_stack = mean_lrnr,
                         d_stack = sl,
                         cv_folds = cv_folds)

  ## 6) d misspecified
  est_misd <- estimators(data = data,
                         delta = delta,
                         g_stack = sl,
                         e_stack = sl,
                         m_stack = sl,
                         b_stack = sl,
                         d_stack = mean_lrnr,
                         cv_folds = cv_folds)

  # bundle estimators in list
  estimates <- list(allc = est_allc, misg = est_misg, mise = est_mise,
                    mism = est_mism, misb = est_misb, misd = est_misd)
  sim_out <- bind_rows(estimates, .id = "sim_type")
  return(sim_out)
}
