library(here)
library(tidyverse)
library(data.table)
source(here("src", "dgp_funs.R"))

# generate data based on DGP components
sim_data <- function(n_obs = 1000) {
    ## baseline covariate -- simple, binary
    w_1 <- rbinom(n_obs, 1, prob = 0.6)
    w_2 <- rbinom(n_obs, 1, prob = 0.3)
    w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
    w <- cbind(w_1, w_2, w_3)
    w_names <- paste("W", seq_len(ncol(w)), sep = "")
    colnames(w) <- w_names

    ## exposure/treatment
    a <- as.numeric(rbinom(n_obs, 1, prob = gfun(1, w)))

    ## mediator-outcome confounder affected by treatment
    l <- rbinom(n_obs, 1, pl(1, a, w))

    ## mediator (possibly multivariate)
    z <- rbinom(n_obs, 1, pz(1, l, a, w))
    z_names <- "Z"

    ## outcome
    y <- rbinom(n_obs, 1, my(z, l, a, w))

    ## construct data for output
    dat <- as.data.table(cbind(W = w, A = a, Z = z, L = l, Y = y))
    return(dat)
}
