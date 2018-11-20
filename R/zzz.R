.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "medshift v", utils::packageDescription("medshift")$Version,
    ": Causal Mediation Analysis for Stochastic Interventions"
  ))
}
