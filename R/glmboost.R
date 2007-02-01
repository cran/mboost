
### 
### Experimental version of gradient boosting with componentwise least
### squares as base learner, i.e., fitting of generalized linear models
### 


### Fitting function
glmboost_fit <- function(object, family = GaussReg(), control = boost_control(),
                      weights = NULL) {

    ### init data and weights
    x <- object$x
    if (control$center) {
        x <- object$center(x)
        ### object$x <- x
    }
    y <- object$yfit
    check_y_family(object$y, family)
    if (is.null(weights)) {
        weights <- object$w
    } else {
        if (NROW(x) == length(weights)) 
            object$w <- weights
        else 
            stop(sQuote("weights"), " is not of length ", NROW(x))
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
    ens <- matrix(NA, nrow = mstop, ncol = 2)
    colnames(ens) <- c("xselect", "coef")

    ### vector of empirical risks for all boosting iterations 
    ### (either in-bag or out-of-bag)
    mrisk <- numeric(mstop)
    mrisk[1:mstop] <- NA

    ### some calculations independent of mstop and memory allocation
    ### for each _column_ of the design matrix x, compute the corresponding
    ### Moore-Penrose inverse (which is a scalar in this case) for the raw
    ### and standardized input variables
    xw <- x * weights
    xtx <- colSums(x^2 * weights)
    MPinv <- (1 / xtx) * t(xw)
    MPinvS <- (1 / sqrt(xtx)) * t(xw)

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit, weights)

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
        u <- ngradient(y, fit, weights)

        ### check if learning is still possible
        if (all(u < 0) || all(u > 0))
            warning("All elements of the negative gradient vector have the same sign in iteration ", m, ".")

        ### evaluate risk, either for the learning sample (inbag)
        ### or the test sample (oobag)
        if (risk == "inbag") mrisk[m] <- riskfct(y, fit, weights)
        if (risk == "oobag") mrisk[m] <- riskfct(y, fit, oobweights)

        ### save the model, i.e., the selected coefficient and variance
        ens[m,] <- c(xselect, coef)
    }

    updatefun <- function(object, control, weights)
        glmboost_fit(object, family = family,
                     control = control, weights = weights)

    RET <- list(ensemble = ens,		### coefficients for selected variables
                fit = fit,		### vector of fitted values
                offset = offset,	### offset
                ustart = ustart,	### first negative gradients
                risk = mrisk,		### empirical risks for m = 1, ..., mstop
                control = control, 	### control parameters
                family = family,	### family object
                response = y, 		### the response variable
                weights = weights,	### weights used for fitting
                update = updatefun,	### a function for fitting with new weights
                MPinv = MPinv 		### Moore-Penrose inverse of x
    )
    ### save learning sample
    if (control$savedata) RET$data <- object
    class(RET) <- c("glmboost", "gb")

    ### prediction function (linear predictor only)
    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {

        if (!is.null(newdata)) {
            if (is.null(object$menv)) {
                if (!is.matrix(newdata) || ncol(newdata) != ncol(x))
                    stop(sQuote("newdata"), " is not a matrix with ",
                         ncol(x), " columns")
                x <- newdata
            } else {
                mf <- object$menv@get("input", data = newdata)
                x <- model.matrix(attr(mf, "terms"), data = mf)
            }
            if (control$center) x <- object$center(x)
        }

        tmp <- RET
        tmp$ensemble <- tmp$ensemble[1:mstop,,drop = TRUE]
        lp <- offset + x %*% coef(tmp)
        if (constraint) lp <- sign(lp) * pmin(abs(lp), 1)    
        return(drop(lp))
    }
    ### function for computing hat matrices of individual predictors
    RET$hat <- function(j) x[,j] %*% MPinv[j, ,drop = FALSE]

    return(RET)
}

### generic method for gradient boosting with componentwise linear models
### for fitting generalized linear models
glmboost <- function(x, ...) UseMethod("glmboost")

### formula interface
glmboost.formula <- function(formula, data = list(), weights = NULL, 
                             contrasts.arg = NULL, ...) {

    ### construct design matrix etc.
    object <- boost_dpp(formula, data, weights, contrasts.arg = contrasts.arg)

    object$center <- function(xmat) {
        cm <- colMeans(object$x)
        num <- which(sapply(object$menv@get("input"), class) == "numeric")
        cm[!attr(object$x, "assign") %in% num] <- 0       
        scale(xmat, center = cm, scale = FALSE)
    }

    ### fit the ensemble
    RET <- glmboost_fit(object, ...)

    RET$call <- match.call()

    return(RET)
}

### matrix interface
glmboost.matrix <- function(x, y, weights = NULL, ...) {

    if (NROW(x) != NROW(y))
        stop("number of observations in", sQuote("x"), "and", 
             sQuote("y"), "differ")
    if (is.null(weights)) weights <- rep(1, NROW(x))
    if (length(weights) != NROW(x))
        stop("number of observations in", sQuote("x"), "and", 
             sQuote("weights"), "differ")

    object <- gb_xyw(x, y, weights)
    object$center <- function(xmat) 
        scale(xmat, center = colMeans(x), scale = FALSE)
    glmboost_fit(object, ...)
}

### methods: coefficients
coef.glmboost <- function(object, ...) {

    ret <- numeric(NCOL(object$data$x))
    xselect <- object$ensemble[,"xselect"]
    for (j in unique(xselect))
        ret[j] <- sum(object$ensemble[xselect == j, "coef"])
    names(ret) <- colnames(object$data$x)
    RET <- ret * object$control$nu
    attr(RET, "offset") <- object$offset
    RET
}

### methods: hatvalues. 
hatvalues.glmboost <- function(model, ...) {

    if (!checkL2(model)) return(hatglm(model))
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
    cat("\t Generalized Linear Models Fitted via Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", mstop(x), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("Offset: ", x$offset, "\n")
    cat("\n")
    cat("Coefficients: \n")
    cf <- coef(x)
    attr(x, "offset") <- NULL
    print(cf)
    cat("\n")
    invisible(x)
}
