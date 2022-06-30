library(here)
library(tidyverse)
library(data.table)
source(here("R", "dgp_funs.R"))
source(here("R", "01_setup_data.R"))

# compute truth based on DGP components
get_sim_truth <- function(n_obs = 1e7, delta) {
    ## simulate data at very large sample size
    data_large <- sim_data(n_obs)

    ## extract components from large data set
    w_names <- str_subset(colnames(data_large), "W")
    z_names <- str_subset(colnames(data_large), "Z")
    w <- data_large[, ..w_names]
    a <- data_large[["A"]]
    l <- data_large[["L"]]
    z <- data_large[[z_names]]
    y <- data_large[["Y"]]

    ## compute EIF components for direct and indirect effects
    eif10 <- pzaw(z, a, w) / pz(z, l, a, w) * (y - my(z, l, a, w)) +
        v(l, a, w) - vbar(a, w) + u(z, a, w)

    eif1d <- gdeltafun(a, w, delta) / gfun(a, w) * pzaw(z, a, w) /
        pz(z, l, a, w) * (y - my(z, l, a, w)) +
        gdeltafun(a, w, delta) / gfun(a, w) *
        (v(l, a, w) - vbar(a, w) + u(z, a, w)) +
        (1 - gdeltafun(a, w, delta) / gfun(a, w)) *
        (ubar(1, w) * gdeltafun(1, w, delta) + ubar(0, w) *
         gdeltafun(0, w, delta))

    eif2d <- gdeltafun(a, w, delta) / gfun(a, w) * pzw(z, w) / pz(z, l, a, w) *
        (y - my(z, l, a, w)) + gdeltafun(a, w, delta) / gfun(a, w) *
        (s(l, a, w) - sbar(a, w)) + u(z, 1, w) * gdeltafun(1, w, delta) +
        u(z, 0, w) * gdeltafun(0, w, delta) + gdeltafun(a, w, delta) /
        gfun(a, w) * (q(a, w) - (q(1, w) * gdeltafun(1, w, delta) + q(0, w) *
                                 gdeltafun(0, w, delta)))

    ## compute truth and efficiency bound approximations based on EIF
    de <- mean(eif10 - eif2d)
    ie <- mean(eif2d - eif1d)
    efbde <- var(eif10 - eif2d)
    efbie <- var(eif2d - eif1d)

    ## output
    out <- list(effect = c("direct", "indirect"),
                truth = c(de, ie),
                eff_bound = c(efbde, efbie))
    return(as_tibble(out))
}
