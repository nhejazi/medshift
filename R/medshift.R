#' Nonparametric estimation of direct effects under mediation
#'
#' @param W ...
#' @param A ...
#' @param Z ...
#' @param Y ...
#' @param g_lrnrs ...
#' @param e_lrnrs ...
#' @param Q_lrnrs ...
#'
#' @importFrom data.table as.data.table
#' @importFrom tmle3 define_node tmle3_Task
#
medshift <- function(W,
                     A,
                     Z,
                     Y,
                     g_lrnrs,
                     e_lrnrs,
                     m_lrnrs) {

  # construct input data structure
  data <- data.table::as.data.table(cbind(Y, Z, A, W))
  w_names <- names(W)
  z_names <- names(Z)

  # construct NPSEM and TMLE task
  npsem <- list(
    tmle3::define_node("W", w_names),
    tmle3::define_node("A", "A", c("W")),
    tmle3::define_node("Z", z_names, c("A", "W")),
    tmle3::define_node("E", "A", c(z_names, "W")),
    tmle3::define_node("Y", "Y", c("Z", "A", "W"))
  )
  tmle_task <- tmle3::tmle3_Task$new(data, npsem = npsem)

  # TODO: functions for estimating likelihood components
}
