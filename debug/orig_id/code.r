rm(list = ls())
library(SuperLearner)
source('utils.r')
source('slfunctions.r')
source('estim.r')

estima <- function(i){
    seed <- task[i, 'seed']
    n <- task[i, 'n']

    set.seed(seed)
    data <- simdata(n)

    Y <- data[, "Y"]
    Z <- data[, substr(names(data), 1, 1) == "Z"]
    L <- data[, "L"]
    A <- data[, "A"]
    W <- data[, substr(names(data), 1, 1) == "W"]

    nuisance <- fitnuisance(data, nfolds = 5)

    nuisance1 <- nuisance2 <- nuisance3 <- nuisance4 <- nuisance5 <- nuisance

    nuisance1$fitg <- mySL(
        Y = A,
        X = W,
        family = stats::binomial(),
        SL.library = c('SL.mean'),
        validRows = nuisance$valSets)
    nuisance2$fite <- mySL(
        Y = A,
        X = data.frame(Z, W),
        family = stats::binomial(),
        SL.library = c('SL.mean'),
        validRows = nuisance$valSets)
    nuisance3$fitm <- mySL(
        Y = Y,
        X = data.frame(Z, A, L, W),
        family = stats::binomial(),
        SL.library = c('SL.mean'),
        validRows = nuisance$valSets)
    nuisance4$fitb <- mySL(
        Y = L,
        X = data.frame(A, W),
        family = stats::binomial(),
        SL.library = c('SL.mean'),
        validRows = nuisance$valSets)
    nuisance5$fitd <- mySL(
        Y = L,
        X = data.frame(Z, A, W),
        family = stats::binomial(),
        SL.library = c('SL.mean'),
        validRows = nuisance$valSets)

    nuisances <- list(nuisance, nuisance1, nuisance2, nuisance3, nuisance4, nuisance5)

    out <- lapply(1:6, function(j) {
        out <- estimators(data, delta = 2, nuisances[[j]], seed = seed,
                          type = j - 1)
        return(out)
        })

    return(Reduce(rbind, out))

}

rep <- 1000
set.seed(48838)
seed <- sample(1e8, rep)
task <- expand.grid(seed = seed, n = cumsum(rep(sqrt(200), 9))^2)
## task <- expand.grid(seed = seed, n = rep(sqrt(200), 4)^2)

funslave <- function(j){

    index <- (1:nrow(task))[(1:nrow(task) - 1) %% 1200 + 1 == j]

    out <- lapply(index, function(k){
        ret <- estima(k)
        return(ret)
    })

    return(out)

}

## n    <- 1e5
## seed <- 29803457
## set.seed(seed)
## data <- simdata(n)

## nuisance <- fitnuisance(data, nfolds = 5)

## z <- Z; a <- A; l <- L; y <- Y; w <- W
## plot(pzaw(z, a, w) / pz(z, l, a, w), b / d)
## abline(0,1)

## plot(gdeltafun(a, w, delta) / gfun(a, w), gdelta / g)
## abline(0,1)

## plot(my(z, l, a, w), m)
## abline(0,1)

## eif10 <- pzaw(z, a, w) / pz(z, l, a, w) * (y - my(z, l, a, w)) +
##     v(l, a, w) - vbar(a, w) + u(z, a, w)

## eif1d <- gdeltafun(a, w, delta) / gfun(a, w) * pzaw(z, a, w) / pz(z, l, a, w) *
##     (y - my(z, l, a, w)) +
##     gdeltafun(a, w, delta) / gfun(a, w) * (v(l, a, w) - vbar(a, w) + u(z, a, w)) +
##     (1 - gdeltafun(a, w, delta) / gfun(a, w)) *
##     (ubar(1, w) * gdeltafun(1, w, delta) + ubar(0, w) * gdeltafun(0, w, delta))

## eif2d <- gdeltafun(a, w, delta) / gfun(a, w) * pzw(z, w) / pz(z, l, a, w) *
##     (y - my(z, l, a, w)) +
##     gdeltafun(a, w, delta) / gfun(a, w) * (s(l, a, w) - sbar(a, w)) +
##     u(z, 1, w) * gdeltafun(1, w, delta) + u(z, 0, w) * gdeltafun(0, w, delta) +
##     gdeltafun(a, w, delta) / gfun(a, w) * (q(a, w) - (q(1, w) * gdeltafun(1, w, delta) +
##                                                       q(0, w) * gdeltafun(0, w, delta)))


## eifde <- eif10 - eif2d
## eifie <- eif2d - eif1d

## plot(eifde, eifD)
## abline(0,1)

## plot(eifie, eifI, col = a+1)
## abline(0,1)

## tmpA <- gdeltafun(a, w, delta) / gfun(a, w) *
##     ((q(a, w) - (q(1, w) * gdeltafun(1, w, delta) + q(0, w) * gdeltafun(0, w, delta))) -
##      (ubar(a, w) - (ubar(1, w) * gdeltafun(1, w, delta) + ubar(0, w) * gdeltafun(0, w, delta))))
## plot(MI * (A - g1), tmp)
## abline(0,1)

## tmpY <- (gdeltafun(a, w, delta) / gfun(a, w) * pzw(z, w) / pz(z, l, a, w) -
##                     gdeltafun(a, w, delta) / gfun(a, w) * pzaw(z, a, w) / pz(z, l, a, w)) *
##                    (y - my(z, l, a, w))
## plot(HI * (Y - m), tmpY)
## abline(0,1)

## tmpZ <- gdeltafun(a, w, delta) / gfun(a, w) *  (u(z, a, w) - ubar(a, w))
## plot(HI * (Y - m), tmpY)
## abline(0,1)

## eifD <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1) + (uA - ubarA) + de
## eifI <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1) + J * (uA - ubarA) + ie

