# packages
renv::activate(here::here())
library(here)
#library(foreach)
#library(future)
#library(doFuture)
#library(doRNG)
library(data.table)
library(tidyverse)
library(hal9001)
library(origami)
library(sl3)
library(SuperLearner)

# load package and ivan's code
devtools::load_all(here("intmedlite"))
source(here("orig_id", "slfunctions.r"))
source(here("orig_id", "utils.r"))
source(here("orig_id", "estim.r"))

# load scripts, parallelization, PRNG
source(here("src", "setup_data.R"))
source(here("src", "fit_estim.R"))
#registerDoFuture()
#plan(multiprocess, workers = 24)

# simulation parameters
set.seed(11249)
ipsi_delta <- 2
cv_folds <- 5L
n_obs <- 5000
n_sim <- 3

# perform simulation across sample sizes
#sim_results <- lapply(n_obs, function(n_obs) {
  ## get results in parallel
  #results <- foreach(this_iter = seq_len(n_sim),
                     #.options.multicore = list(preschedule = FALSE),
                     #.errorhandling = "remove") %dorng% {
    # generate simulated data
    data_sim <- sim_data(n_obs = n_obs)

    # run estimators from original implementation
    nuisance_allc_id <- fitnuisance(
      data = as.data.frame(data_sim), nfolds = cv_folds
    )
    nuisance_misb_id <- nuisance_misd_id <- nuisance_mise_id <-
      nuisance_misg_id <- nuisance_mism_id <- nuisance_allc_id

    # set nuisance parameter mis-configurations from ivan's code, based on
    # loop over difference nuisance configurations as in orig_id/code.R
    nuisance_misg_id$fitg <- mySL(
        Y = data_sim$A,
        X = data_sim[, .SD, .SDcols = patterns("W")],
        family = stats::binomial(),
        SL.library = c("SL.mean"),
        validRows = nuisance_allc_id$valSets
    )
    nuisance_mise_id$fite <- mySL(
        Y = data_sim$A,
        X = cbind(data_sim[, .SD, .SDcols = patterns("W")], data_sim$Z),
        family = stats::binomial(),
        SL.library = c("SL.mean"),
        validRows = nuisance_allc_id$valSets
    )
    nuisance_mism_id$fitm <- mySL(
        Y = data_sim$Y,
        X = cbind(data_sim[, .SD, .SDcols = patterns("W")], data_sim$Z,
                  data_sim$A, data_sim$L),
        family = stats::binomial(),
        SL.library = c("SL.mean"),
        validRows = nuisance_allc_id$valSets
    )
    nuisance_misb_id$fitb <- mySL(
        Y = data_sim$L,
        X = cbind(data_sim[, .SD, .SDcols = patterns("W")], data_sim$A),
        family = stats::binomial(),
        SL.library = c("SL.mean"),
        validRows = nuisance_allc_id$valSets
    )
    nuisance_misd_id$fitd <- mySL(
        Y = data_sim$L,
        X = cbind(data_sim[, .SD, .SDcols = patterns("W")], data_sim$Z,
                  data_sim$A),
        family = stats::binomial(),
        SL.library = c("SL.mean"),
        validRows = nuisance_allc_id$valSets
    )

    # fit estimators across grid of nuiance paramter configurations
    est_allc_id <- estimators(
      data = as.data.frame(data_sim), delta = ipsi_delta,
      nuisance = nuisance_allc_id, seed = NA_real_, type = "allc"
    )
    est_misg_id <- estimators(
      data = as.data.frame(data_sim), delta = ipsi_delta,
      nuisance = nuisance_misg_id, seed = NA_real_, type = "misg"
    )
    est_mise_id <- estimators(
      data = as.data.frame(data_sim), delta = ipsi_delta,
      nuisance = nuisance_mise_id, seed = NA_real_, type = "mise"
    )
    est_mism_id <- estimators(
      data = as.data.frame(data_sim), delta = ipsi_delta,
      nuisance = nuisance_mism_id, seed = NA_real_, type = "mism"
    )
    est_misb_id <- estimators(
      data = as.data.frame(data_sim), delta = ipsi_delta,
      nuisance = nuisance_misb_id, seed = NA_real_, type = "misb"
    )
    est_misd_id <- estimators(
      data = as.data.frame(data_sim), delta = ipsi_delta,
      nuisance = nuisance_misd_id, seed = NA_real_, type = "misd"
    )
    est_id <- rbind(est_allc_id, est_misg_id, est_mise_id, est_mism_id,
                    est_misb_id, est_misd_id) %>%
      setDT() %>%
      mutate(
        parameter = case_when(parameter == "DE" ~ "direct",
                              parameter == "IE" ~ "indirect"),
        estimator = tolower(estimator)
      ) %>%
      rename(sim_type = type, n_obs = n) %>%
      select(-seed) %>%
      relocate(sim_type, parameter, estimator, estimate, ses, n_obs)


    # run estimators from sl3-based re-implementation
    glmnet_learner <- Lrnr_glmnet$new()
    est_nh <- fit_estimators(
      data = data_sim, delta = ipsi_delta, cv_folds = cv_folds
    )
    est_nh <- setDT(est_nh)


    # organize results across corresponding implementations
    est_combined <- bind_rows(est_id, est_nh, .id = "impl") %>%
      mutate(impl = case_when(impl == 1 ~ "id", impl == 2 ~ "nh"))

    est_diff_summary <- est_combined %>%
      group_by(sim_type, parameter, estimator) %>%
      summarize(diff(estimate), diff(ses))
    print(setDT(est_diff_summary))

    # output simulation results
    #return(est_out)
  #}

  # concatenate iterations
  #results_out <- bind_rows(results, .id = "sim_iter")
  #return(results_out)
#})

# save results to file
#names(sim_results) <- paste("n", n_obs, sep = "_")
#timestamp <- str_replace_all(Sys.time(), " ", "_")
#saveRDS(object = sim_results,
        #file = here("data", paste0("compare_id_nh_", timestamp, ".rds")))
