
### 
### Experimental version of gradient boosting with componentwise least
### squares as base learner, i.e., fitting of generalized linear models
### 


### Fitting function
glmboost_fit <- function(object, family = GaussReg(), control = boost_control(),
                      weights = NULL) {

    ### init data and weights
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

    ### unweighted problem
    WONE <- (max(abs(weights - 1)) < .Machine$double.eps)
    if (!family@weights && !WONE)
        stop(sQuote("family"), " is not able to deal with weights") 

    ### rescale weights (because of the AIC criterion)
    ### <FIXME> is this correct with zero weights??? </FIXME>
    weights <- rescale_weights(weights)

    ### the ensemble
    ens <- matrix(NA, nrow = mstop, ncol = 3)
    colnames(ens) <- c("xselect", "coef", "risk")
    brisk <- NA

    ### some calculations independent of mstop and memory allocation
    ### for each _column_ of the design matrix x, compute the corresponding
    ### Moore-Penrose inverse (which is a scalar in this case) for the raw
    ### and standardized input variables
    xw <- x * weights
    xtx <- colSums(x^2 * weights)
    MPinv <- (1 / xtx) * t(xw)
    MPinvS <- (1 / sqrt(xtx)) * t(xw)

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit)

    ### start boosting iteration
    for (m in 1:mstop) {
  
        ### fit least squares to residuals _componentwise_, i.e.,
        ### compute regression coefficients for each _standardized_
        ### input variable and select the best variable
        xselect <- which.max(abs(MPinvS %*% u))

        ### estimate regression coefficient (not standardized)
        coef <- MPinv[xselect,] %*% u

        ### update step
        fit <- fit + (nu * coef) * x[,xselect]

        ### L2 boost with constraints (binary classification)
        if (constraint)
            fit <- sign(fit) * pmin(abs(fit), 1)

        ### negative gradient vector, the new `residuals'
        u <- ngradient(y, fit)   

        ### AIC precomputations
        if (AIC) brisk <- risk(y, fit, weights)

        ### save the model, i.e., the selected coefficient and variance
        ens[m,] <- c(xselect, coef, brisk)
    }

    RET <- list(ensemble = ens, control = control, fit = fit, 
                ustart = ustart,  offset = offset,
                MPinv = MPinv, family = family, weights = weights)
    if (control$savedata) RET$data <- object
    class(RET) <- c("glmboost", "gb")

    ### prediction function (linear predictor)
    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {

        if (!is.null(newdata)) {
            mf <- object$menv@get("input", data = newdata)
            x <- model.matrix(attr(mf, "terms"), data = mf)
        }

        tmp <- RET
        tmp$ensemble <- tmp$ensemble[1:mstop,,drop = TRUE]
        lp <- offset + x %*% coef(tmp)
        if (constraint) lp <- sign(lp) * pmin(abs(lp), 1)    
        return(drop(lp))
    }
    RET$hat <- function(j) x[,j] %*% MPinv[j, ,drop = FALSE]

    return(RET)
}

### generic method for gradient boosting with componentwise linear models
### for fitting generalized linear models
glmboost <- function(x, ...) UseMethod("glmboost")

### formula interface
glmboost.formula <- function(formula, data = list(), weights = NULL, ...) {

    ### construct design matrix etc.
    object <- boost_dpp(formula, data, weights)

    ### fit the ensemble
    RET <- glmboost_fit(object, ...)

    RET$call <- match.call()

    return(RET)
}

### matrix interface
glmboost.matrix <- function(x, y, weights = NULL, ...) {

    if (length(y) != nrow(x))
        stop("number of observations in", sQuote("x"), "and", 
             sQuote("y"), "differ")
    if (is.null(weights)) weights <- rep(1, length(y))
    if (length(weights) != nrow(x))
        stop("number of observations in", sQuote("x"), "and", 
             sQuote("weights"), "differ")

    object <- gb_xyw(x, y, weights)
    glmboost_fit(object, ...)
}

### methods: coefficients
coef.glmboost <- function(object, ...) {

    ret <- sapply(1:ncol(object$data$x), function(j)
        sum(object$ensemble[object$ensemble[,"xselect"] == j, "coef"]))
    names(ret) <- colnames(object$data$x)
    ret * object$control$nu
}

### methods: hatvalues. For L_2 loss ONLY!
hatvalues.glmboost <- function(model, ...) {

    if (!checkL2(model)) return(hatvalues.gb(model, ...))
    xf <- t(model$MPinv) * model$control$nu
    op <- .Call("R_trace_glmboost", model$data$x, xf,
                as.integer(model$ensemble[, "xselect"]),
                PACKAGE = "mboost")
    RET <- diag(op[[1]])
    attr(RET, "hatmatrix") <- op[[1]]  
    attr(RET, "trace") <- op[[2]] 
    RET
}

### methods: print
print.glmboost <- function(x, ...) {

    cat("\n")
    cat("\t Generalized Linear Models Fitted via Gradient Boosting) \n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", nrow(x$ensemble), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("\n")
    cat("Coefficients: \n")
    print(coef(x))
    cat("\n")
    invisible(x)
}
