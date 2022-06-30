library(SuperLearner)

SL.caretRF <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'rf',  tuneLength = 20,
             trControl =  caret::trainControl(method = "cv", number = 5,
                                              search = 'random',
                                              verboseIter = TRUE), ...)
}

SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree',
             tuneLength = 30,
             trControl =  caret::trainControl(method = "cv", number = 5,
                                              search = 'random',
                                              verboseIter = TRUE), ...)
}

SL.myglmnet <- function (Y, X, newX, family, obsWeights, id, alpha = 1,
                         nfolds = 10, nlambda = 500, useMin = TRUE,
                         loss = "mse", ...) {
    .SL.require("glmnet")
    if (!is.matrix(X)) {
        formula <- as.formula(paste('~ -1 + .^', eval(ncol(X))))
        X <- model.matrix(formula, X)
        newX <- model.matrix(formula, newX)
    }
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
        lambda = NULL, type.measure = loss, nfolds = nfolds,
        family = "gaussian", lambda.min.ratio = 1 / nrow(X),
        alpha = alpha, nlambda = nlambda,
        ...)
    pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
        "lambda.min", "lambda.1se"))
    fit <- list(object = fitCV, useMin = useMin)
    class(fit) <- "SL.myglmnet"
    out <- list(pred = pred, fit = fit)
    return(out)
}

SL.myglm <- function (Y, X, newX, family, obsWeights, model = TRUE, ...) {
    if (is.matrix(X)) {
        X = as.data.frame(X)
    }
    formula <- as.formula(paste('Y ~ .^', eval(ncol(X))))
    fit.glm <- glm(formula, data = X, family = family, weights = obsWeights,
        model = model)
    if (is.matrix(newX)) {
        newX = as.data.frame(newX)
    }
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred = pred, fit = fit)
    return(out)
}

SL.caretXGB <- function(Y, X, newX, family, obsWeights, ...) {
    SL.caret(Y, X, newX, family, obsWeights, method = 'xgbTree',
             tuneLength = 30,
             trControl =  caret::trainControl(method = "cv", number = 5,
                                              search = 'random',
                                              verboseIter = TRUE), ...)
}

predict.SL.myglmnet <- function(object, newdata, remove_extra_cols = FALSE,
                                add_missing_cols = FALSE, ...) {
    .SL.require("glmnet")
    if (!is.matrix(newdata)) {
        formula <- as.formula(paste('~ -1 + .^', eval(ncol(newdata))))
        newdata <- model.matrix(formula, newdata)
    }
    original_cols = rownames(object$object$glmnet.fit$beta)
    if (remove_extra_cols) {
        extra_cols = setdiff(colnames(newdata), original_cols)
        if (length(extra_cols) > 0) {
            warning(paste("Removing extra columns in prediction data:",
                paste(extra_cols, collapse = ", ")))
            newdata = newdata[, !colnames(newdata) %in% extra_cols,
                drop = F]
        }
    }
    if (add_missing_cols) {
        missing_cols = setdiff(original_cols, colnames(newdata))
        if (length(missing_cols) > 0) {
            warning(paste("Adding missing columns in prediction data:",
                paste(missing_cols, collapse = ", ")))
            new_cols = matrix(0, nrow = nrow(newdata),
                              ncol = length(missing_cols))
            colnames(new_cols) = missing_cols
            newdata = cbind(newdata, new_cols)
            newdata = newdata[, original_cols]
        }
    }
    pred <- predict(object$object, newx = newdata, type = "response",
        s = ifelse(object$useMin, "lambda.min", "lambda.1se"))
    return(pred)
}

environment(SL.myglm) <- environment(SuperLearner)
environment(SL.myglmnet) <- environment(SuperLearner)
environment(predict.SL.myglmnet) <- environment(SuperLearner)
