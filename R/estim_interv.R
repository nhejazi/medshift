utils::globalVariables(c("..w_names", "..z_names"))

#' Compute Cross-Fitted One-Step Estimator
#'
#' @param data ...
#' @param delta ...
#' @param g_learners ...
#' @param e_learners ...
#' @param b_learners ...
#' @param d_learners ...
#' @param m_learners ...
#' @param cv_folds ...
#'
#' @import data.table origami sl3 stats stringr tibble
interv_est_onestep <- function(data,
                               delta,
                               g_learners,
                               e_learners,
                               b_learners,
                               d_learners,
                               m_learners,
                               cv_folds = 5L) {

  # create folds for use with origami::cross_validate
  folds <- origami::make_folds(data,
    fold_fun = origami::folds_vfold,
    V = cv_folds,
    cluster_ids = data[["ids"]]
  )

  # perform the cv_eif procedure on a per-fold basis
  cv_eif_results <- origami::cross_validate(
    cv_fun = interv_cv_eif,
    folds = folds,
    data = data,
    delta = delta,
    g_learners = g_learners,
    e_learners = e_learners,
    b_learners = b_learners,
    d_learners = d_learners,
    m_learners = m_learners,
    w_names = w_names,
    z_names = z_names,
    use_future = FALSE,
    .combine = FALSE
  )

  # compute one-step estimate and inference
  osde <- mean(eifD)
  osie <- mean(eifI)
  seosde <- sd(eifD) / sqrt(n_obs)
  seosie <- sd(eifI) / sqrt(n_obs)

  # output
  out <- list(
    parameter = c("direct", "indirect", "direct", "indirect"),
    estimator = c("os", "os", "tmle", "tmle"),
    estimate = c(osde, osie, tmlede, tmleie),
    ses = c(seosde, seosie, setmlede, setmleie)
  )
  return(tibble::as_tibble(out))


}



#' Compute Cross-Validated TML estimator
#'
#' @param data ...
#' @param delta ...
#' @param g_learners ...
#' @param e_learners ...
#' @param b_learners ...
#' @param d_learners ...
#' @param m_learners ...
#' @param cv_folds ...
#' @param tiltmod_tol ...
#'
#' @import data.table origami sl3 stats stringr tibble
interv_est_tmle <- function(data,
                            delta,
                            g_learners,
                            e_learners,
                            b_learners,
                            d_learners,
                            m_learners,
                            cv_folds = 5L,
                            tiltmod_tol = 10) {

  # create folds for use with origami::cross_validate
  folds <- origami::make_folds(data,
    fold_fun = origami::folds_vfold,
    V = cv_folds,
    cluster_ids = data[["ids"]]
  )

  # perform the cv_eif procedure on a per-fold basis
  cv_eif_results <- origami::cross_validate(
    cv_fun = interv_cv_eif,
    folds = folds,
    data = data,
    delta = delta,
    g_learners = g_learners,
    e_learners = e_learners,
    b_learners = b_learners,
    d_learners = d_learners,
    m_learners = m_learners,
    w_names = w_names,
    z_names = z_names,
    use_future = FALSE,
    .combine = FALSE
  )



  # now, compute the TMLE
  stopcrit <- FALSE
  iter <- 1

  # iterative TMLE
  while (!stopcrit) {
    suppressWarnings(
      tiltm <- glm(Y ~ -1 + offset(qlogis(bound(m))) + HD + HI,
        family = "binomial", start = c(0, 0)
      )
    )
    tiltm <- check_fluc(tiltm, tiltmod_tol)
    m <- predict(tiltm, type = "response")

    suppressWarnings(
      tiltb <- glm(L ~ -1 + offset(qlogis(bound(b1A))) + KD + KI,
        family = "binomial", start = c(0, 0)
      )
    )
    tiltb <- check_fluc(tiltb, tiltmod_tol)
    b1A <- predict(tiltb, type = "response")

    suppressWarnings(
      tiltg <- glm(A ~ -1 + offset(qlogis(bound(g1))) + MD + MI,
        family = "binomial", start = c(0, 0)
      )
    )
    tiltg <- check_fluc(tiltg, tiltmod_tol)
    g1 <- predict(tiltg, type = "response")

    HD11 <- with(bdm_cv_cf_preds, b11 / d11 * (1 - g1delta / e1))
    HI11 <- with(bdm_cv_cf_preds, b11 / d11 * g1delta * (1 / e1 - 1 / g1))
    HD10 <- with(bdm_cv_cf_preds, b10 / d10 * (1 - (1 - g1delta) / (1 - e1)))
    HI10 <- with(bdm_cv_cf_preds, b10 / d10 * (1 - g1delta) *
                 (1 / (1 - e1) - 1 / (1 - g1)))

    HD01 <- with(bdm_cv_cf_preds, (1 - b11) / (1 - d11) * (1 - g1delta / e1))
    HI01 <- with(bdm_cv_cf_preds, (1 - b11) / (1 - d11) * g1delta *
                 (1 / e1 - 1 / g1))
    HD00 <- with(bdm_cv_cf_preds, (1 - b10) / (1 - d10) *
                 (1 - (1 - g1delta) / (1 - e1)))
    HI00 <- with(bdm_cv_cf_preds, (1 - b10) / (1 - d10) *
                 (1 - g1delta) * (1 / (1 - e1) - 1 / (1 - g1)))

    m11 <- with(bdm_cv_cf_preds, plogis(qlogis(bound(m11)) +
                                        coef(tiltm)[1] * HD11 +
                                        coef(tiltm)[2] * HI11))
    m10 <- with(bdm_cv_cf_preds, plogis(qlogis(bound(m10)) +
                                        coef(tiltm)[1] * HD10 +
                                        coef(tiltm)[2] * HI10))
    m01 <- with(bdm_cv_cf_preds, plogis(qlogis(bound(m01)) +
                                        coef(tiltm)[1] * HD01 +
                                        coef(tiltm)[2] * HI01))
    m00 <- with(bdm_cv_cf_preds, plogis(qlogis(bound(m00)) +
                                        coef(tiltm)[1] * HD00 +
                                        coef(tiltm)[2] * HI00))

    g <- A * g1 + (1 - A) * (1 - g1)
    b <- L * b1A + (1 - L) * (1 - b1A)
    g1delta <- ipsi_delta(g1, delta)
    gdelta <- A * g1delta + (1 - A) * (1 - g1delta)

    HD <- b / d * (1 - gdelta / e)
    HI <- b / d * gdelta * (1 / e - 1 / g)
    KD <- v1 - v0 - gdelta / g * (s1 - s0)
    KI <- gdelta / g * (s1 - s0 - v1 + v0)
    MD <- - g1delta * (1 - g1delta) / (g1 * (1 - g1)) * q2
    MI <- g1delta * (1 - g1delta) / (g1 * (1 - g1)) * (q2 - q1)

    de <- ubar1 * g1 + ubar0 * (1 - g1) - (u1 * g1delta + u0 *
      (1 - g1delta))
    ie <- (u1 - ubar1) * g1delta + (u0 - ubar0) * (1 - g1delta)

    eifDp <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1)
    eifIp <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1)

    ## iterate iterator
    iter <- iter + 1
    ## note: interesting stopping criterion
    stopcrit <- max(abs(c(mean(eifDp), mean(eifIp)))) <
      max(c(sd(eifDp), sd(eifIp))) * log(n_obs) / n_obs | iter > 6
  }

  J <- - gdelta / g
  J1 <- - g1delta / g1
  J0 <- - (1 - g1delta) / (1 - g1)

  m1A <- A * m11 + (1 - A) * m10
  m0A <- A * m01 + (1 - A) * m00

  uA <- m1A * b1A + m0A * (1 - b1A)
  u1 <- m11 * bdm_cv_cf_preds$b11 + m01 * (1 - bdm_cv_cf_preds$b11)
  u0 <- m10 * bdm_cv_cf_preds$b10 + m00 * (1 - bdm_cv_cf_preds$b10)

  suppressWarnings(
    tiltu <- glm(uA ~ offset(qlogis(bound(ubarA))) + J, family = "binomial")
  )
  tiltu <- check_fluc(tiltu, tiltmod_tol)
  ubar1 <- plogis(coef(tiltu)[1] + qlogis(bound(ubar1)) +
    coef(tiltu)[2] * J1)
  ubar0 <- plogis(coef(tiltu)[1] + qlogis(bound(ubar0)) +
    coef(tiltu)[2] * J0)
  ubarA <- A * ubar1 + (1 - A) * ubar0

  de <- ubar1 * g1 + ubar0 * (1 - g1) - (u1 * g1delta + u0 * (1 - g1delta))
  ie <- (u1 - ubar1) * g1delta + (u0 - ubar0) * (1 - g1delta)

  # updated EIFs from TML estimation
  eifD <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1) +
    (uA - ubarA) + de
  eifI <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1) + J *
    (uA - ubarA) + ie

  # compute TML estimate and inference
  tmlede <- mean(de)
  tmleie <- mean(ie)
  setmlede <- sd(eifD) / sqrt(n_obs)
  setmleie <- sd(eifI) / sqrt(n_obs)

  # output
  out <- list(
    parameter = c("direct", "indirect", "direct", "indirect"),
    estimator = c("os", "os", "tmle", "tmle"),
    estimate = c(osde, osie, tmlede, tmleie),
    ses = c(seosde, seosie, setmlede, setmleie)
  )
  return(tibble::as_tibble(out))
}
