
### 
### Experimental version of gradient boosting with componentwise least
### smoothing splines
### 

### hat matrix for smoothing splines
hatMatTH <- function(x, w = NULL, df = 5) {
    n <- nrow(x)
    indx <- diag(n)
    apply(indx, 2, function(y) 
        predict(smooth.spline(x = x, y = y, w = w, df = df, 
                keep.data = FALSE)$fit, x = x)$y)
}

predict.smooth.spline.fit <- stats:::predict.smooth.spline.fit

### Fitting function
gamboost_fit <- function(object, dfbase = 4, family = GaussReg(), 
                      control = boost_control(), weights = NULL) {

    ### data
    x <- object$x
    y <- object$yfit
    if (is.null(weights)) {
        weights <- object$w
    } else {
        if (length(y) == length(weights))
            object$w <- weights
        else 
            stop(sQuote("weights"), " is not of length ", length(y))
    }

    ### hyper parameters
    mstop <- control$mstop
    AIC <- control$risk
    constraint <- control$constraint
    nu <- control$nu

    ### extract negative gradient and risk functions
    ngradient <- family@ngradient
    risk <- family@risk

    ### inputs with more than four unique values
    pindx <- (1:ncol(x))[apply(x, 2, function(xx) length(unique(xx)) > 4)]

    ### unweighted problem
    WONE <- (max(abs(weights - 1)) < .Machine$double.eps)
    if (!family@weights && !WONE)
        stop(sQuote("family"), " is not able to deal with weights")

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- rescale_weights(weights)

    ### the ensemble
    ens <- matrix(NA, nrow = mstop, ncol = 2)
    colnames(ens) <- c("xselect", "risk")
    ensss <- vector(mode = "list", length = mstop)
    brisk <- NA

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit)

    ### start boosting iteration
    for (m in 1:mstop) {
  
        sums <- 0
        xselect <- 0
        ### fit least squares to residuals _componentwise_
        for (i in pindx) {
            if (WONE) {
                ss <- try(smooth.spline(x = x[,i], y = u, df = dfbase, 
                                        keep.data = FALSE))
            } else {
                ss <- try(smooth.spline(x = x[,i], y = u, w = weights, 
                                        df = dfbase, keep.data = FALSE))
            }
            if (inherits(ss, "try-error")) next()
            tsums <- sum((predict(ss$fit, x[,i])$y - u)^2)
            if (tsums < sums || sums == 0) {
                sums <- tsums
                xselect <- i
                basess <- ss$fit
            }
        }

        if (xselect == 0) 
            stop("could not fit base learner in boosting iteration ", m)

        ### update step
        fit <- fit + nu * predict(basess, x = x[,xselect])$y

        ### L2 boost with constraints (binary classification)
        if (constraint)
            fit <- sign(fit) * pmin(abs(fit), 1)

        ### negative gradient vector, the new `residuals'
        u <- ngradient(y, fit)   

        ### AIC precomputations
        if (AIC) brisk <- risk(y, fit, weights)

        ### save the model, i.e., the selected coefficient and variance
        ens[m,] <- c(xselect, brisk)
        ensss[[m]] <- basess

    }

    RET <- list(ensemble = ens, ensembless = ensss,
                control = control, fit = fit, offset = offset, 
                ustart = ustart, family = family, dfbase = dfbase, weights = weights)
    if (control$savedata) RET$data <- object
    class(RET) <- c("gamboost", "gb")

    ### prediction function
    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {

        if (!is.null(newdata)) {
            mf <- object$menv@get("input", data = newdata)
            x <- model.matrix(attr(mf, "terms"), data = mf)
        }

        lp <- offset
        for (m in 1:mstop)
            lp <- lp + nu * predict(ensss[[m]], 
                                    x = x[,ens[m,"xselect"]])$y
        if (constraint) lp <- sign(lp) * pmin(abs(lp), 1)    
        return(lp)
    }
    RET$hat <- function(j) hatMatTH(x[,j,drop = FALSE],
                                    w = weights, df = dfbase)
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

### matrix interface
gamboost.matrix <- function(x, y, weights = NULL, ...) {

    if (length(y) != nrow(x))
        stop("number of observations in", sQuote("x"), "and",
             sQuote("y"), "differ")
    if (is.null(weights)) weights <- rep(1, length(y))
    if (length(weights) != nrow(x))
        stop("number of observations in", sQuote("x"), "and",
             sQuote("weights"), "differ")

    object <- gb_xyw(x, y, weights)
    gamboost_fit(object, ...)
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
    cat("Number of boosting iterations: mstop =", nrow(x$ensemble), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Degree of freedom: ", x$dfbase, "\n")
    cat("\n")
    invisible(x)

}
