#' Incremental Propensity Score Intervention
#'
#' compute shift where gn := g(A = 1 | W)
#'
#' @param gn_est ...
#' @param delta ...
#'
#' @export
#
ipsi_shift <- function(gn_est, delta) {
  gn_est_shifted <- (delta * gn_est) / (delta * gn_est + (1 - gn_est))
  return(gn_est_shifted)
}

#' Modified Treatment Policy Shift
#'
#' compute shift where d(A,W) = A + delta
#'
#' @param A ...
#' @param delta ...
#'
#' @export
#
mtp_shift <- function(A, delta) {
  # just assume unbounded support for A|W for now...
  A_shifted <- A + delta
  return(A_shifted)
}
