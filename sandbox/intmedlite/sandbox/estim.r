library(SuperLearner)
source('utils.r')
source('slfunctions.r')

fitnuisance <- function(data, nfolds,
                        sl_lib = c('SL.myglmnet', 'SL.myglm', 'SL.mean')) {

    candidatesg <- candidatese <- candidatesm <- candidatesb <-
        candidatesd <- sl_lib
    ## extract data
    Y <- data[, "Y"]
    Z <- data[, substr(names(data), 1, 1) == "Z"]
    L <- data[, "L"]
    A <- data[, "A"]
    W <- data[, substr(names(data), 1, 1) == "W"]

    ## observations and cross-validation
    n <- length(A)
    valSets <- split(sample(seq_len(n)), rep(seq_len(nfolds), length = n))

    ## fit nuisance functions
    fitg <- mySL(
        Y = A,
        X = W,
        family = stats::binomial(),
        SL.library = candidatesg,
        validRows = valSets
    )
    fite <- mySL(
        Y = A,
        X = data.frame(Z, W),
        family = stats::binomial(),
        SL.library = candidatese,
        validRows = valSets
    )
    fitm <- mySL(
        Y = Y,
        X = data.frame(Z, L, A, W),
        family = binomial(),
        SL.library = candidatesm,
        validRows = valSets
    )
    fitb <- mySL(
        Y = L,
        X = data.frame(A, W),
        family = stats::binomial(),
        SL.library = candidatesb,
        validRows = valSets
    )
    fitd <- mySL(
        Y = L,
        X = data.frame(Z, A, W),
        family = stats::binomial(),
        SL.library = candidatesd,
        validRows = valSets
    )

    return(list(fitg = fitg, fite = fite, fitm = fitm, fitb = fitb,
                fitd = fitd, valSets = valSets))
}

estimators <- function(data, delta, nuisance, seed, type,
                       sl_lib = c('SL.myglmnet', 'SL.myglm', 'SL.mean')) {

    candidatesv  <- candidatess <- candidatesq <- candidatesubar <- sl_lib

    ## extract data
    Y <- data[, "Y"]
    Z <- data[, substr(names(data), 1, 1) == "Z"]
    L <- data[, "L"]
    A <- data[, "A"]
    W <- data[, substr(names(data), 1, 1) == "W"]

    n <- length(A)

    g1  <- truncate(nuisance$fitg$SL.predict[, 1])
    e1  <- truncate(nuisance$fite$SL.predict[, 1])
    b1A <- nuisance$fitb$SL.predict[, 1]
    d1A <- truncate(nuisance$fitd$SL.predict[, 1])

    g <- A * g1 + (1 - A) * (1 - g1)
    e <- A * e1 + (1 - A) * (1 - e1)
    b <- L * b1A + (1 - L) * (1 - b1A)
    d <- L * d1A + (1 - L) * (1 - d1A)

    b11  <- predict(nuisance$fitb, newdata = data.frame(A = 1, W))$pred[, 1]
    b10  <- predict(nuisance$fitb, newdata = data.frame(A = 0, W))$pred[, 1]
    d11  <- truncate(predict(nuisance$fitd, newdata = data.frame(Z, A = 1, W))$pred[, 1])
    d10  <- truncate(predict(nuisance$fitd, newdata = data.frame(Z, A = 0, W))$pred[, 1])
    m11  <- predict(nuisance$fitm, newdata = data.frame(Z, L = 1, A = 1, W))$pred[, 1]
    m10  <- predict(nuisance$fitm, newdata = data.frame(Z, L = 1, A = 0, W))$pred[, 1]
    m01  <- predict(nuisance$fitm, newdata = data.frame(Z, L = 0, A = 1, W))$pred[, 1]
    m00  <- predict(nuisance$fitm, newdata = data.frame(Z, L = 0, A = 0, W))$pred[, 1]

    m <- nuisance$fitm$SL.predict[, 1]
    m1A <- A * m11 + (1 - A) * m10
    m0A <- A * m01 + (1 - A) * m00

    uA <- m1A * b1A + m0A * (1 - b1A)
    u1 <- m11 * b11 + m01 * (1 - b11)
    u0 <- m10 * b10 + m00 * (1 - b10)

    vout <- m * b / d
    sout <- m * b / d * g / e

    if(sd(vout) < .Machine$double.eps * 10) {
        v1 <- v0 <- rep(mean(vout), length(vout))
    } else {
        fitv <- mySL(
            Y = vout,
            X = data.frame(L, A, W),
            family = stats::gaussian(),
            SL.library = candidatesv,
            validRows = nuisance$valSets
        )
        v1 <- predict(fitv, newdata = data.frame(L = 1, A, W))$pred[, 1]
        v0 <- predict(fitv, newdata = data.frame(L = 0, A, W))$pred[, 1]
    }

    if(sd(sout) < .Machine$double.eps * 10) {
        s1 <- s0 <- rep(mean(sout), length(sout))
    } else {
        fits <- mySL(
            Y = sout,
            X = data.frame(L, A, W),
            family = stats::gaussian(),
            SL.library = candidatess,
            validRows = nuisance$valSets
        )
        s1 <- predict(fits, newdata = data.frame(L = 1, A, W))$pred[, 1]
        s0 <- predict(fits, newdata = data.frame(L = 0, A, W))$pred[, 1]
    }

    if(sd(uA) < .Machine$double.eps * 10) {
        ubar1 <- ubar0 <- rep(mean(uA), length(uA))
    } else {
        fitubar <- mySL(
            Y = uA,
            X = data.frame(A, W),
            family = stats::gaussian(),
            SL.library = candidatesubar,
            validRows = nuisance$valSets
        )
        ubar1 <- predict(fitubar, newdata = data.frame(A = 1, W))$pred[, 1]
        ubar0 <- predict(fitubar, newdata = data.frame(A = 0, W))$pred[, 1]
    }

    ubarA  <- A * ubar1 + (1 - A) * ubar0
    q1    <- ubar1 - ubar0

    if(sd(u1 - u0) < .Machine$double.eps * 10) {
        q2  <- rep(mean(u1 - u0), length(u1))
    } else {
        fitq2 <- mySL(
            Y = u1 - u0,
            X = data.frame(W),
            family = stats::gaussian(),
            SL.library = candidatesq,
            validRows = nuisance$valSets
        )
        q2  <- fitq2$SL.predict[, 1]
    }

    g1delta <- gdelta1(g1, delta)
    gdelta <- A * g1delta + (1 - A) * (1 - g1delta)

    HD <- b / d * (1 - gdelta / e)
    HI <- b / d * gdelta * (1 / e -  1 / g)
    KD <- v1 - v0 - gdelta / g * (s1 - s0)
    KI <- gdelta / g * (s1 - s0 - v1 + v0)
    MD <- - g1delta * (1 - g1delta) / (g1 * (1 - g1)) * q2
    MI <- g1delta * (1 - g1delta) / (g1 * (1 - g1)) * (q2 - q1)
    J  <- - gdelta / g

    de <- ubar1 * g1 + ubar0 * (1 - g1) - (u1 * g1delta + u0 * (1 - g1delta))
    ie <- (u1 - ubar1) * g1delta + (u0 - ubar0) * (1 - g1delta)
    eifD <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1) + (uA - ubarA) + de
    eifI <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1) + J * (uA - ubarA) + ie

    osde <- mean(eifD)
    osie <- mean(eifI)
    seosde <- sd(eifD) / sqrt(n)
    seosie <- sd(eifI) / sqrt(n)

    ## now, compute the TMLE
    stopcrit <- FALSE
    iter <- 1

    ## iterative TMLE
    while (!stopcrit) {
        tiltm <- glm(Y ~ 0 + offset(qlogis(bound(m)))   + HD + HI, family = binomial())
        tiltb <- glm(L ~ 0 + offset(qlogis(bound(b1A))) + KD + KI, family = binomial())
        tiltg <- glm(A ~ 0 + offset(qlogis(bound(g1)))  + MD + MI, family = binomial())

        m   <- predict(tiltm, type = 'response')
        b1A <- predict(tiltb, type = 'response')
        g1  <- predict(tiltg, type = 'response')

        HD11 <- b11 / d11 * (1 - g1delta / e1)
        HI11 <- b11 / d11 * g1delta * (1 / e1 -  1 / g1)
        HD10 <- b10 / d10 * (1 - (1 - g1delta) / (1 - e1))
        HI10 <- b10 / d10 * (1 - g1delta) * (1 / (1 - e1) -  1 / (1 - g1))

        HD01 <- (1 - b11) / (1 - d11) * (1 - g1delta / e1)
        HI01 <- (1 - b11) / (1 - d11) * g1delta * (1 / e1 -  1 / g1)
        HD00 <- (1 - b10) / (1 - d10) * (1 - (1 - g1delta) / (1 - e1))
        HI00 <- (1 - b10) / (1 - d10) * (1 - g1delta) * (1 / (1 - e1) -  1 / (1 - g1))

        m11 <- plogis(qlogis(bound(m11)) + coef(tiltm)[1] * HD11 + coef(tiltm)[2] * HI11)
        m10 <- plogis(qlogis(bound(m10)) + coef(tiltm)[1] * HD10 + coef(tiltm)[2] * HI10)
        m01 <- plogis(qlogis(bound(m01)) + coef(tiltm)[1] * HD01 + coef(tiltm)[2] * HI01)
        m00 <- plogis(qlogis(bound(m00)) + coef(tiltm)[1] * HD00 + coef(tiltm)[2] * HI00)

        g <- A * g1 + (1 - A) * (1 - g1)
        b <- L * b1A + (1 - L) * (1 - b1A)
        g1delta <- gdelta1(g1, delta)
        gdelta <- A * g1delta + (1 - A) * (1 - g1delta)

        HD <- b / d * (1 - gdelta / e)
        HI <- b / d * gdelta * (1 / e -  1 / g)
        KD <- v1 - v0 - gdelta / g * (s1 - s0)
        KI <- gdelta / g * (s1 - s0 - v1 + v0)
        MD <- - g1delta * (1 - g1delta) / (g1 * (1 - g1)) * q2
        MI <- g1delta * (1 - g1delta) / (g1 * (1 - g1)) * (q2 - q1)

        de <- ubar1 * g1 + ubar0 * (1 - g1) - (u1 * g1delta + u0 * (1 - g1delta))
        ie <- (u1 - ubar1) * g1delta + (u0 - ubar0) * (1 - g1delta)

        eifDp <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1)
        eifIp <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1)

        ## iterate iterator
        iter <- iter + 1
        ## note: interesting stopping criterion
        stopcrit <- max(abs(c(mean(eifDp), mean(eifIp)))) <
            max(c(sd(eifDp), sd(eifIp))) * log(n) / n | iter > 6
    }

    J   <- - gdelta / g
    J1  <- - g1delta / g1
    J0  <- - (1 - g1delta) / (1 - g1)

    m1A <- A * m11 + (1 - A) * m10
    m0A <- A * m01 + (1 - A) * m00

    uA <- m1A * b1A + m0A * (1 - b1A)
    u1 <- m11 * b11 + m01 * (1 - b11)
    u0 <- m10 * b10 + m00 * (1 - b10)

    tiltu <- glm(uA ~ offset(qlogis(bound(ubarA))) + J, family = binomial())
    ubar1 <- plogis(coef(tiltu)[1] + qlogis(bound(ubar1)) + coef(tiltu)[2] * J1)
    ubar0 <- plogis(coef(tiltu)[1] + qlogis(bound(ubar0)) + coef(tiltu)[2] * J0)
    ubarA  <- A * ubar1 + (1 - A) * ubar0

    de <- ubar1 * g1 + ubar0 * (1 - g1) - (u1 * g1delta + u0 * (1 - g1delta))
    ie <- (u1 - ubar1) * g1delta + (u0 - ubar0) * (1 - g1delta)

    eifD <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1) + (uA - ubarA) + de
    eifI <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1) + J * (uA - ubarA) + ie

    tmlede <- mean(de)
    tmleie <- mean(ie)
    setmlede <- sd(eifD) / sqrt(n)
    setmleie <- sd(eifI) / sqrt(n)

    out <- data.frame(
        parameter = c('DE', 'IE', 'DE', 'IE'),
        estimator = c('OS', 'OS', 'TMLE', 'TMLE'),
        estimate  = c(osde, osie, tmlede, tmleie),
        ses       = c(seosde, seosie, setmlede, setmleie),
        n         = n,
        seed      = seed,
        type      = type)

    return(out)
}
