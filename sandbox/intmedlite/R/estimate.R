utils::globalVariables(c("..w_names", "..z_names"))

#' Compute One-step and TML estimators
#'
#' @import data.table origami sl3 stats stringr tibble
#'
#' @param data ...
#' @param delta ...
#' @param g_stack ...
#' @param e_stack ...
#' @param m_stack ...
#' @param b_stack ...
#' @param d_stack ...
#' @param cv_folds ...
#' @param tiltmod_tol ...
estimators <- function(data,
                       delta,
                       g_stack,
                       e_stack,
                       m_stack,
                       b_stack,
                       d_stack,
                       cv_folds = 5L,
                       tiltmod_tol = 10) {
  ## extract data
  data <- data.table::as.data.table(data)
  w_names <- stringr::str_detect(names(data), "W")
  z_names <- stringr::str_detect(names(data), "Z")
  Y <- data[["Y"]]
  Z <- data[, ..z_names]
  L <- data[["L"]]
  A <- data[["A"]]
  W <- data[, ..w_names]
  n_obs <- nrow(data)

  ## get names of data components
  names_z <- stringr::str_subset(colnames(data), "Z")
  names_w <- stringr::str_subset(colnames(data), "W")

  ## instantiate HAL and GLM for pseudo-outcome regressions
  hal_learner <- Lrnr_hal9001$new(max_degree = 4L,
                                  n_folds = 10L,
                                  fit_type = "glmnet",
                                  family = "gaussian",
                                  cv_select = TRUE,
                                  reduce_basis = 0.01,
                                  return_lasso = FALSE,
                                  type.measure = "mse",
                                  standardize = FALSE,
                                  lambda.min.ratio = 1e-4,
                                  nlambda = 1000L,
                                  yolo = FALSE)

  ## fit nuisance function estimators
  nuisance_fits <- fit_nuisance(data, g_stack, e_stack, m_stack,
                                b_stack, d_stack, cv_folds)
  trn_idx <- lapply(nuisance_fits$folds, `[[`, "training_set")
  val_idx <- lapply(nuisance_fits$folds, `[[`, "validation_set")
  reorder_by_val_idx <- order(do.call(c, val_idx))

  ## truncate predictions
  g1 <- truncate(nuisance_fits$cv_preds$g_preds[reorder_by_val_idx])
  e1 <- truncate(nuisance_fits$cv_preds$e_preds[reorder_by_val_idx])
  b1A <- nuisance_fits$cv_preds$b_preds[reorder_by_val_idx]
  d1A <- truncate(nuisance_fits$cv_preds$d_preds[reorder_by_val_idx])

  ## compute "natural" estimates
  g <- unname(unlist(A * g1 + (1 - A) * (1 - g1)))
  e <- unname(unlist(A * e1 + (1 - A) * (1 - e1)))
  b <- unname(unlist(L * b1A + (1 - L) * (1 - b1A)))
  d <- unname(unlist(L * d1A + (1 - L) * (1 - d1A)))
  m <- nuisance_fits$cv_preds$m_preds[reorder_by_val_idx]

  ## compute cross-validated counterfactual predictions
  # NOTE: notation renamed, so...
  bdm_cv_cf_preds <- lapply(seq_along(val_idx), function(fold_idx) {
    ## predict for b
    b0_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, A := 0][],
      covariates = c(names_w, "A"),
      outcome = "L",
      outcome_type = "binomial"
    )
    b1_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, A := 1][],
      covariates = c(names_w, "A"),
      outcome = "L",
      outcome_type = "binomial"
    )
    b11 <- nuisance_fits$cv_fits$b_fit[[fold_idx]]$predict(b1_task)
    b10 <- nuisance_fits$cv_fits$b_fit[[fold_idx]]$predict(b0_task)

    ## predict for d
    d0_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, A := 0][],
      covariates = c(names_w, "A", names_z),
      outcome = "L",
      outcome_type = "binomial"
    )
    d1_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, A := 1][],
      covariates = c(names_w, "A", names_z),
      outcome = "L",
      outcome_type = "binomial"
    )
    d11 <- nuisance_fits$cv_fits$d_fit[[fold_idx]]$predict(d1_task)
    d10 <- nuisance_fits$cv_fits$d_fit[[fold_idx]]$predict(d0_task)

    ## predict for m
    m00_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, `:=`(A = 0, L = 0)][],
      covariates = c(names_w, "A", "L", names_z),
      outcome = "Y",
      outcome_type = "binomial"
    )
    m01_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, `:=`(A = 1, L = 0)][],
      covariates = c(names_w, "A", "L", names_z),
      outcome = "Y",
      outcome_type = "binomial"
    )
    m10_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, `:=`(A = 0, L = 1)][],
      covariates = c(names_w, "A", "L", names_z),
      outcome = "Y",
      outcome_type = "binomial"
    )
    m11_task <- sl3::sl3_Task$new(
      data = data[val_idx[[fold_idx]], ][, `:=`(A = 1, L = 1)][],
      covariates = c(names_w, "A", "L", names_z),
      outcome = "Y",
      outcome_type = "binomial"
    )
    m00 <- nuisance_fits$cv_fits$m_fit[[fold_idx]]$predict(m00_task)
    m01 <- nuisance_fits$cv_fits$m_fit[[fold_idx]]$predict(m01_task)
    m10 <- nuisance_fits$cv_fits$m_fit[[fold_idx]]$predict(m10_task)
    m11 <- nuisance_fits$cv_fits$m_fit[[fold_idx]]$predict(m11_task)

    ## return holdout-specific predictions
    cf_preds <- data.table::as.data.table(list(
      b10, b11, d10, d11, m00, m01,
      m10, m11
    ))
    setnames(cf_preds, c(
      "b10", "b11", "d10", "d11", "m00", "m01", "m10",
      "m11"
    ))
    return(cf_preds)
  })
  bdm_cv_cf_preds <- data.table::rbindlist(bdm_cv_cf_preds)
  bdm_cv_cf_preds <- bdm_cv_cf_preds[reorder_by_val_idx, ]

  m1A <- unname(unlist(A * bdm_cv_cf_preds$m11 +
    (1 - A) * bdm_cv_cf_preds$m10))
  m0A <- unname(unlist(A * bdm_cv_cf_preds$m01 +
    (1 - A) * bdm_cv_cf_preds$m00))

  uA <- unname(unlist(m1A * b1A + m0A * (1 - b1A)))
  u1 <- with(bdm_cv_cf_preds, m11 * b11 + m01 * (1 - b11))
  u0 <- with(bdm_cv_cf_preds, m10 * b10 + m00 * (1 - b10))

  vout <- m * b / d
  sout <- m * b / d * g / e

  if (sd(vout) < .Machine$double.eps * 10) {
    v1 <- v0 <- rep(mean(vout), length(vout))
  } else {
    # make data for estimating v
    v_data <- data.table::as.data.table(list(W, A, L, vout))
    data.table::setnames(v_data, c(names_w, "A", "L", "vout"))

    # cross-validated counterfactual predictions
    v_cv_cf_preds <- lapply(seq_along(val_idx), function(fold_idx) {
      v_task <- sl3::sl3_Task$new(
        data = v_data[trn_idx[[fold_idx]], ],
        covariates = c(names_w, "A", "L"),
        outcome = "vout",
        outcome_type = "gaussian"
      )

      # safely fit pseudo-outcome regression
      v_fit <- hal_learner$train(v_task)

      # counterfactual tasks
      v_task_L1 <- sl3_Task$new(
        data = v_data[val_idx[[fold_idx]], ][, L := 1][],
        covariates = c(names_w, "A", "L"),
        outcome = "vout",
        outcome_type = "gaussian"
      )
      v_task_L0 <- sl3_Task$new(
        data = v_data[val_idx[[fold_idx]], ][, L := 0][],
        covariates = c(names_w, "A", "L"),
        outcome = "vout",
        outcome_type = "gaussian"
      )

      # generate counterfactual predictions
      v1 <- v_fit$predict(v_task_L1)
      v0 <- v_fit$predict(v_task_L0)
      out <- data.table::as.data.table(list(v1, v0))
      data.table::setnames(out, c("v1", "v0"))
      return(out)
    })
    v_cv_cf_preds <- data.table::rbindlist(v_cv_cf_preds)
    v_cv_cf_preds <- v_cv_cf_preds[reorder_by_val_idx, ]
    v1 <- v_cv_cf_preds$v1
    v0 <- v_cv_cf_preds$v0
  }

  if (sd(sout) < .Machine$double.eps * 10) {
    s1 <- s0 <- rep(mean(sout), length(sout))
  } else {
    # make data for estimating s
    s_data <- data.table::as.data.table(list(W, A, L, sout))
    data.table::setnames(s_data, c(names_w, "A", "L", "sout"))

    # cross-validated counterfactual predictions
    s_cv_cf_preds <- lapply(seq_along(val_idx), function(fold_idx) {
      s_task <- sl3_Task$new(
        data = s_data[trn_idx[[fold_idx]], ],
        covariates = c(names_w, "A", "L"),
        outcome = "sout",
        outcome_type = "gaussian"
      )

      # safely fit pseudo-outcome regression
      s_fit <- hal_learner$train(s_task)

      # counterfactual tasks
      s_task_L1 <- sl3_Task$new(
        data = s_data[val_idx[[fold_idx]], ][, L := 1][],
        covariates = c(names_w, "A", "L"),
        outcome = "sout",
        outcome_type = "gaussian"
      )
      s_task_L0 <- sl3_Task$new(
        data = s_data[val_idx[[fold_idx]], ][, L := 0][],
        covariates = c(names_w, "A", "L"),
        outcome = "sout",
        outcome_type = "gaussian"
      )

      # generate counterfactual predictions
      s1 <- s_fit$predict(s_task_L1)
      s0 <- s_fit$predict(s_task_L0)
      out <- data.table::as.data.table(list(s1, s0))
      data.table::setnames(out, c("s1", "s0"))
      return(out)
    })
    s_cv_cf_preds <- data.table::rbindlist(s_cv_cf_preds)
    s_cv_cf_preds <- s_cv_cf_preds[reorder_by_val_idx, ]
    s1 <- s_cv_cf_preds$s1
    s0 <- s_cv_cf_preds$s0
  }

  if (sd(uA) < .Machine$double.eps * 10) {
    ubar1 <- ubar0 <- rep(mean(uA), length(uA))
  } else {
    # make data for estimating ubar
    ubar_data <- data.table::as.data.table(list(W, A, uA))
    data.table::setnames(ubar_data, c(names_w, "A", "uA"))

    # cross-validated counterfactual predictions
    ubar_cv_cf_preds <- lapply(seq_along(val_idx), function(fold_idx) {
      ubar_task <- sl3_Task$new(
        data = ubar_data[trn_idx[[fold_idx]], ],
        covariates = c(names_w, "A"),
        outcome = "uA",
        outcome_type = "gaussian"
      )

      # safely fit pseudo-outcome regression
      ubar_fit <- hal_learner$train(ubar_task)

      # counterfactual tasks
      ubar_task_A1 <- sl3_Task$new(
        data = ubar_data[val_idx[[fold_idx]], ][, A := 1][],
        covariates = c(names_w, "A"),
        outcome = "uA",
        outcome_type = "gaussian"
      )
      ubar_task_A0 <- sl3_Task$new(
        data = ubar_data[val_idx[[fold_idx]], ][, A := 0][],
        covariates = c(names_w, "A"),
        outcome = "uA",
        outcome_type = "gaussian"
      )

      # generate counterfactual predictions
      ubar1 <- ubar_fit$predict(ubar_task_A1)
      ubar0 <- ubar_fit$predict(ubar_task_A0)
      out <- data.table::as.data.table(list(ubar1, ubar0))
      data.table::setnames(out, c("ubar1", "ubar0"))
      return(out)
    })
    ubar_cv_cf_preds <- data.table::rbindlist(ubar_cv_cf_preds)
    ubar_cv_cf_preds <- ubar_cv_cf_preds[reorder_by_val_idx, ]
    ubar1 <- ubar_cv_cf_preds$ubar1
    ubar0 <- ubar_cv_cf_preds$ubar0
  }

  # NOTE: the next bit is only relevant to IPSI's
  ubarA <- A * ubar1 + (1 - A) * ubar0
  q1 <- ubar1 - ubar0

  if (sd(u1 - u0) < .Machine$double.eps * 10) {
    q2 <- rep(mean(u1 - u0), length(u1))
  } else {
    # make data for estimating q2
    q2_data <- data.table::as.data.table(list(W, u1 - u0))
    data.table::setnames(q2_data, c(names_w, "u_diff"))

    # cross-validated predictions
    q2_cv_preds <- lapply(seq_along(val_idx), function(fold_idx) {
      q2_train_task <- sl3_Task$new(
        data = q2_data[trn_idx[[fold_idx]], ],
        covariates = names_w,
        outcome = "u_diff",
        outcome_type = "gaussian"
      )
      q2_valid_task <- sl3_Task$new(
        data = q2_data[val_idx[[fold_idx]], ],
        covariates = names_w,
        outcome = "u_diff",
        outcome_type = "gaussian"
      )

      # safely fit pseudo-outcome regression
      q2_fit <- hal_learner$train(q2_train_task)
      q2 <- q2_fit$predict(q2_valid_task)
    })
    q2 <- do.call(c, q2_cv_preds)[reorder_by_val_idx]
  }

  # compute IPSI shifts
  g1delta <- gdelta1(g1, delta)
  gdelta <- A * g1delta + (1 - A) * (1 - g1delta)

  # compute clever covariates
  HD <- b / d * (1 - gdelta / e)
  HI <- b / d * gdelta * (1 / e - 1 / g)
  KD <- v1 - v0 - gdelta / g * (s1 - s0)
  KI <- gdelta / g * (s1 - s0 - v1 + v0)
  MD <- -g1delta * (1 - g1delta) / (g1 * (1 - g1)) * q2
  MI <- g1delta * (1 - g1delta) / (g1 * (1 - g1)) * (q2 - q1)
  J <- -gdelta / g

  # compute (in)direct effects and EIF
  de <- ubar1 * g1 + ubar0 * (1 - g1) - (u1 * g1delta + u0 * (1 - g1delta))
  ie <- (u1 - ubar1) * g1delta + (u0 - ubar0) * (1 - g1delta)
  eifD <- HD * (Y - m) + KD * (L - b1A) + MD * (A - g1) +
    (uA - ubarA) + de
  eifI <- HI * (Y - m) + KI * (L - b1A) + MI * (A - g1) + J *
    (uA - ubarA) + ie

  # compute one-step estimate and inference
  osde <- mean(eifD)
  osie <- mean(eifI)
  seosde <- sd(eifD) / sqrt(n_obs)
  seosie <- sd(eifI) / sqrt(n_obs)

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
    g1delta <- gdelta1(g1, delta)
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

  #browser()
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
    ses = c(seosde, seosie, setmlede, setmleie),
    n_obs = n_obs
  )
  return(tibble::as_tibble(out))
}
