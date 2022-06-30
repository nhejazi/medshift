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
n_obs <- 400
n_sim <- 3

# perform simulation across sample sizes
sim_results <- lapply(n_obs, function(n_obs) {
  # get results in parallel
  results <- foreach(this_iter = seq_len(n_sim),
                     .options.multicore = list(preschedule = FALSE),
                     .errorhandling = "remove") %dorng% {
    # generate simulated data
    data_sim <- sim_data(n_obs = n_obs)

    # run estimators from original implementation
    est_id <- NULL
    nuisance_id <- fitnuisance(data = as.data.frame(data_sim), nfolds = 5L)
    # TODO: add calls to ivan's code to loop over difference nuisance
    #       configurations as hard-coded in orig_id/code.R
    est_id <- estimators(data = as.data.frame(data_sim), delta = ipsi_delta,
                         nuisance = nuisance_id, seed = NA_real_,
                         type = "allc")
    est_id <- setDT(est_id)

    # run estimators from sl3-based re-implementation
    est_nh <- fit_estimators(data = data_sim, delta = ipsi_delta)
    est_nh <- setDT(est_nh)

    # TODO: organize results across corresponding implementations
    #return(est_out)
  }

  # concatenate iterations
  #results_out <- bind_rows(results, .id = "sim_iter")
  #return(results_out)
})

# save results to file
#names(sim_results) <- paste("n", n_obs, sep = "_")
#timestamp <- str_replace_all(Sys.time(), " ", "_")
#saveRDS(object = sim_results,
        #file = here("data", paste0("compare_id_nh_", timestamp, ".rds")))
