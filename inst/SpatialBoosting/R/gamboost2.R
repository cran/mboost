
### 
### Experimental version of gradient boosting with componentwise least
### smoothing splines
### 

### Fitting function
gamboost_fit <- function(object, family = GaussReg(), 
                         control = boost_control(), weights = NULL) {

    y <- object$yfit
    check_y_family(object$y, family)
    if (is.null(weights)) {
        weights <- object$w
    } else {
        if (NROW(y) == length(weights))
            object$w <- weights
        else 
            stop(sQuote("weights"), " is not of length ", NROW(y))
    }

    ### hyper parameters
    mstop <- control$mstop
    risk <- control$risk
    constraint <- control$constraint
    nu <- control$nu

    ### extract negative gradient and risk functions
    ngradient <- family@ngradient
    riskfct <- family@risk

    ### unweighted problem
    WONE <- (max(abs(weights - 1)) < .Machine$double.eps)
    if (!family@weights && !WONE)
        stop(sQuote("family"), " is not able to deal with weights")

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- rescale_weights(weights)
    oobweights <- as.numeric(weights == 0)

    ### the ensemble
    ens <- matrix(NA, nrow = mstop, ncol = 1)
    colnames(ens) <- "xselect"
    ensss <- vector(mode = "list", length = mstop)

    ### vector of empirical risks for all boosting iterations
    ### (either in-bag or out-of-bag)
    mrisk <- numeric(mstop)
    mrisk[1:mstop] <- NA   

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

    xnames <- ls(envir = object$menv@env)
    xnames <- xnames[!(xnames %in% c("response", "designMatrix", "responseMatrix"))]
    x <- sapply(xnames, function(xn) object$menv@get(xn))

    ### start boosting iteration
    for (m in 1:mstop) {
  
        sums <- 0
        xselect <- 0
        ### fit least squares to residuals _componentwise_
        for (i in 1:length(x)) {
            ss <- try(attr(x[[i]], "fit")(u))
            if (inherits(ss, "try-error")) next
            tsums <- sum((attr(x[[i]], "fitted")(ss) - u)^2)
            if (tsums < sums || sums == 0) {
                sums <- tsums
                xselect <- i
                basess <- ss
            }
        }

        if (xselect == 0) 
            stop("could not fit base learner in boosting iteration ", m)

        ### update step
        fit <- fit + nu * attr(x[[xselect]], "fitted")(basess)

        ### L2 boost with constraints (binary classification)
        if (constraint)
            fit <- sign(fit) * pmin(abs(fit), 1)

        ### negative gradient vector, the new `residuals'
        u <- ngradient(y, fit, weights)

        ### evaluate risk, either for the learning sample (inbag)
        ### or the test sample (oobag)
        if (risk == "inbag") mrisk[m] <- riskfct(y, fit, weights)
        if (risk == "oobag") mrisk[m] <- riskfct(y, fit, oobweights)

        ### save the model, i.e., the selected coefficient and variance
        ens[m,] <- xselect
        ensss[[m]] <- basess

    }

    updatefun <- function(object, control, weights) 
        gamboost_fit(object, family = family,
                     control = control, weights = weights)

    RET <- list(ensemble = ens,         ### selected variables 
                ensembless = ensss,	### list of smooth.spline fits
                fit = fit,              ### vector of fitted values
                offset = offset,        ### offset
                ustart = ustart,        ### first negative gradients
                risk = mrisk,           ### empirical risks for m = 1, ..., mstop
                control = control,      ### control parameters   
                family = family,        ### family object
                response = y,           ### the response variable
                weights = weights,      ### weights used for fitting   
                update = updatefun     ### a function for fitting with new weights
    )
    ### save learning sample
    if (control$savedata) {
        RET$data <- object
        RET$data$x <- x
    }

    ### prediction function (linear predictor only)
    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {

        if (!is.null(newdata))
            newdata <- x

        lp <- offset
        for (m in 1:mstop)
            lp <- lp + nu * attr(x[[ens[m,]]], "predict")(ensss[[m]], newdata = x)
        if (constraint) lp <- sign(lp) * pmin(abs(lp), 1)    
        return(lp)
    }

    ### function for computing hat matrices of individual predictors
    RET$hat <- function(j) attr(ensss[[which(ens[,1] == j)[1]]], "hatmatrix")()

    class(RET) <- c("gamboost", "gb")
    return(RET)
}

### generic method for gradient boosting with componentwise smoothing splines
### for fitting generalized additive models
gamboost <- function(x, ...) UseMethod("gamboost")

### formula interface
gamboost.formula <- function(formula, data = list(), weights = NULL, ...) {

    ### construct design matrix etc.
    object <- boost_dpp(formula, data, weights)

    ### fit the ensemble
    RET <- gamboost_fit(object, ...)

    RET$call <- match.call()

    return(RET)
}

### methods: print
print.gamboost <- function(x, ...) {

    cat("\n")
    cat("\t Generalized Additive Models Fitted via Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("\n")
    invisible(x)

}
