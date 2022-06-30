bound <- function(x, p = 0.0001) {
    pmax(pmin(x, 1 - p), p)
}

truncate <- function(x, p = 0.01) {
    pmin(pmax(x, p), 1 - p)
}

gfun <- function(a, w) {
    pscore <- plogis(2 - 1 / (- rowSums(w) / 4 + 1))
    return(a * pscore + (1 - a) * (1 - pscore))
}

gdeltafun <- function(a, w, delta){
    pscore <- delta * gfun(1, w) / (delta * gfun(1, w) + 1 - gfun(1, w))
    return(a * pscore + (1 - a) * (1 - pscore))
}

gdelta1 <- function(g, delta){
    delta * g / (delta * g + 1 - g)
}

pl <- function(l, a, w) {
    prob1 <- plogis(rowMeans(w) - a - log(2) + 0.2)
    return(l * prob1 + (1 - l) * (1 - prob1))
}

pz <- function(z, l, a, w) {
    prob1 <- plogis(rowSums(log(3) * w[, -3]) + a - l)
    return(z * prob1 + (1 - z) * (1 - prob1))
}

my <- function(z, l, a, w) {
    plogis(1 - 3 / (2 + rowSums(w) / 3 - l - 3 * a + z))
    ## plogis((rowSums(w) - l + 3 * a * z))
}

pzaw <- function(z, a, w) pz(z, 1, a, w) * pl(1, a, w) + pz(z, 0, a, w) * pl(0, a, w)
pzw  <- function(z, w) pzaw(z, 1, w) * gfun(1, w) + pzaw(z, 0, w) * gfun(0, w)

u <- function(z, a, w) my(z, 1, a, w) * pl(1, a, w) + my(z, 0, a, w) * pl(0, a, w)
v <- function(l, a, w) my(1, l, a, w) * pzaw(1, a, w) + my(0, l, a, w) * pzaw(0, a, w)
s <- function(l, a, w) my(1, l, a, w) * pzw(1, w) + my(0, l, a, w) * pzw(0, w)

ubar <- function(a, w) u(1, a, w) * pzaw(1, a, w) + u(0, a, w) * pzaw(0, a, w)
sbar <- function(a, w) s(1, a, w) * pl(1, a, w) + s(0, a, w) * pl(0, a, w)
vbar <- ubar
q <- function(a, w) u(1, a, w) * pzw(1, w) + u(0, a, w) * pzw(0, w)

simdata <- function(n_obs = 1000) {
    ## helper functions for nuisance parameter estimation
    ## baseline covariate -- simple, binary
    w_1 <- rbinom(n_obs, 1, prob = 0.6)
    w_2 <- rbinom(n_obs, 1, prob = 0.3)
    w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
    w <- cbind(w_1, w_2, w_3)
    w_names <- paste("W", seq_len(ncol(w)), sep = "")

    ## exposure/treatment
    a <- as.numeric(rbinom(n_obs, 1, prob = gfun(1, w)))

    ## mediator-outcome confounder affected by treatment
    l <- rbinom(n_obs, 1, pl(1, a, w))

    ## mediator (possibly multivariate)
    z <- rbinom(n_obs, 1, pz(1, l, a, w))
    z_names <- "Z"

    ## outcome
    y <- rbinom(n_obs, 1, my(z, l, a, w))

    colnames(w) <- w_names
    ## construct data for output
    dat <- as.data.frame(cbind(W = w, A = a, Z = z, L = l, Y = y))
    return(dat)
}


truth <- function(n_obs = 1e7, delta = 2){

    w_1 <- rbinom(n_obs, 1, prob = 0.6)
    w_2 <- rbinom(n_obs, 1, prob = 0.3)
    w_3 <- rbinom(n_obs, 1, prob = pmin(0.2 + (w_1 + w_2) / 3, 1))
    w <- cbind(w_1, w_2, w_3)
    w_names <- paste("W", seq_len(ncol(w)), sep = "")

    ## exposure/treatment
    a <- as.numeric(rbinom(n_obs, 1, prob = gfun(1, w)))

    ## mediator-outcome confounder affected by treatment
    l <- rbinom(n_obs, 1, pl(1, a, w))

    ## mediator (possibly multivariate)
    z <- rbinom(n_obs, 1, pz(1, l, a, w))
    z_names <- "Z"

    ## outcome
    y <- rbinom(n_obs, 1, my(z, l, a, w))

    eif10 <- pzaw(z, a, w) / pz(z, l, a, w) * (y - my(z, l, a, w)) +
        v(l, a, w) - vbar(a, w) + u(z, a, w)

    eif1d <- gdeltafun(a, w, delta) / gfun(a, w) * pzaw(z, a, w) / pz(z, l, a, w) *
        (y - my(z, l, a, w)) +
        gdeltafun(a, w, delta) / gfun(a, w) * (v(l, a, w) - vbar(a, w) + u(z, a, w)) +
        (1 - gdeltafun(a, w, delta) / gfun(a, w)) *
        (ubar(1, w) * gdeltafun(1, w, delta) + ubar(0, w) * gdeltafun(0, w, delta))

    eif2d <- gdeltafun(a, w, delta) / gfun(a, w) * pzw(z, w) / pz(z, l, a, w) *
        (y - my(z, l, a, w)) +
        gdeltafun(a, w, delta) / gfun(a, w) * (s(l, a, w) - sbar(a, w)) +
        u(z, 1, w) * gdeltafun(1, w, delta) + u(z, 0, w) * gdeltafun(0, w, delta) +
        gdeltafun(a, w, delta) / gfun(a, w) * (q(a, w) - (q(1, w) * gdeltafun(1, w, delta) +
                                                          q(0, w) * gdeltafun(0, w, delta)))

    efbde <- var(eif10 - eif2d)
    efbie <- var(eif2d - eif1d)

    de <- mean(eif10 - eif2d)
    ie <- mean(eif2d - eif1d)

    return(data.frame(parameter = c('DE', 'IE'),
                      truth = c(de, ie),
                      eff_bound = c(efbde, efbie)))

}


