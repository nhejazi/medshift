#' Incremental Propensity Score Intervention
#'
#' @param g ...
#' @param delta ...
gdelta1 <- function(g, delta) {
  delta * g / (delta * g + 1 - g)
}

