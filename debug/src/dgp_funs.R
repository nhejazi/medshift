# helper functions for nuisance parameter estimation
gfun <- function(a, w) {
    pscore <- plogis(2 - 1 / (- rowSums(w) / 4 + 1))
    return(a * pscore + (1 - a) * (1 - pscore))
}
gdeltafun <- function(a, w, delta) {
    pscore <- delta * gfun(1, w) / (delta * gfun(1, w) + 1 - gfun(1, w))
    return(a * pscore + (1 - a) * (1 - pscore))
}
gdelta1 <- function(g, delta) {
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

pzaw <- function(z, a, w) {
    pz(z, 1, a, w) * pl(1, a, w) + pz(z, 0, a, w) * pl(0, a, w)
}
pzw  <- function(z, w) {
    pzaw(z, 1, w) * gfun(1, w) + pzaw(z, 0, w) * gfun(0, w)
}

u <- function(z, a, w) {
    my(z, 1, a, w) * pl(1, a, w) + my(z, 0, a, w) * pl(0, a, w)
}
v <- function(l, a, w) {
    my(1, l, a, w) * pzaw(1, a, w) + my(0, l, a, w) * pzaw(0, a, w)
}
s <- function(l, a, w) {
    my(1, l, a, w) * pzw(1, w) + my(0, l, a, w) * pzw(0, w)
}

ubar <- function(a, w) {
    u(1, a, w) * pzaw(1, a, w) + u(0, a, w) * pzaw(0, a, w)
}
sbar <- function(a, w) {
    s(1, a, w) * pl(1, a, w) + s(0, a, w) * pl(0, a, w)
}
vbar <- ubar
q <- function(a, w) {
    u(1, a, w) * pzw(1, w) + u(0, a, w) * pzw(0, w)
}
