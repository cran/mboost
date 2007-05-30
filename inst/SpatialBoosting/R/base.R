
mybs2 <- function(x, z = NULL, df = 5, knots = 20, degree = 3, differences = 2) {

    xname <- deparse(substitute(x))
    if (differences < 1 || differences > 3) 
        stop(sQuote("differences"), " are not in 1:3")
    if (df < differences)
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if (is.factor(x)) ols(x)
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")

    if (length(knots) == 1) {
        knots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = knots + 2) 
        knots <- knots[2:(length(knots) - 1)]
    }
    X <- bs(x, knots = knots, degree = degree, intercept = TRUE)
    if (!is.null(z)) {
        zname <- deparse(substitute(z))
        X <- X * z
    }

    K <- diff(diag(ncol(X)), differences = differences)
    K <- crossprod(K, K) 

    if (FALSE) {
        Rmat <- qr.R(qr(X))
        D <- forwardsolve(t(Rmat),K)
        D <- forwardsolve(t(Rmat),t(D))
        evals <- eigen(D, only.values = TRUE)$values
        XtX <- crossprod(X, X)

        ### Eigenwerte?
        foo <- function(lambda) 
            (sum(1 / (1 + lambda * evals)) - df)^2
        fooopt <- optimize(foo, interval = c(0, 100000))
        lambda <- fooopt$minimum
        Xsolve <- tcrossprod(solve(XtX + lambda * K), X)
        rm(K)
    }
    XtX <- crossprod(X, X)
    ### Eigenwerte?
    foo <- function(lambda) 
        (sum(diag(solve(XtX + lambda * K) %*% XtX)) - df)^2
    fooopt <- optimize(foo, interval = c(0, 100000))
    lambda <- fooopt$minimum
    Xsolve <- tcrossprod(solve(XtX + lambda * K), X)
    rm(XtX)
    rm(K)

    fit <- function(u) list(coef = Xsolve %*% u)
    fitted <- function(object) X %*% object$coef
    predict <- function(object, newdata = NULL) {

        if (is.null(newdata)) return(fitted(object))

        newX <- predict(X, newx = newdata[[xname]])
        if (!is.null(z)) newX <- newX * newdata[[zname]]
        newX %*% object$coef
    }
    hatmatrix <- function() X %*% Xsolve

    attr(X, "fit") <- fit
    attr(X, "predict") <- predict
    attr(X, "fitted") <- fitted
    attr(X, "hatmatrix") <- hatmatrix
    attr(X, "lambda") <- lambda
    return(X)
}

ols <- function(x) {

    xname <- deparse(substitute(x))
    X <- model.matrix(~ x)

    Xsolve <- qr(X)

    fit <- function(u) list(coef = qr.coef(Xsolve, u))
    fitted <- function(object) X %*% object$coef
    predict <- function(object, newdata = NULL) {

        if (is.null(newdata)) return(fitted(object))

        newX <- model.matrix(as.formula(paste("~ ", xname, sep = "")), 
                             newdata = newdata)
        newX %*% object$coef
    }
    hatmatrix <- function() X %*% tcrossprod(solve(crossprod(X, X)), X)

    attr(X, "fit") <- fit
    attr(X, "predict") <- predict
    attr(X, "fitted") <- fitted
    attr(X, "hatmatrix") <- hatmatrix
    attr(X, "lambda") <- lambda
    return(X)
}

spatial <- function(x, y, z = NULL, df = 5, xknots = 20, yknots = 20, 
                    degree = 3, differences = 2) {

    xname <- deparse(substitute(x))
    yname <- deparse(substitute(y))
    if (differences < 1 || differences > 3) 
        stop(sQuote("differences"), " are not in 1:3")
    if (df < differences)
        stop(sQuote("df"), " is less than ", sQuote("differences"))
    if (length(unique(x)) < 6)
        stop(sQuote(xname), " has less than 6 unique values")
    if (length(unique(y)) < 6)
        stop(sQuote(yname), " has less than 6 unique values")

    if (length(xknots) == 1) {
        xknots <- seq(from = min(x, na.rm = TRUE), 
                     to = max(x, na.rm = TRUE), length = xknots + 2) 
        xknots <- xknots[2:(length(xknots) - 1)]
    }
    if (length(yknots) == 1) {
        yknots <- seq(from = min(y, na.rm = TRUE), 
                     to = max(y, na.rm = TRUE), length = yknots + 2) 
        yknots <- yknots[2:(length(yknots) - 1)]
    }
    Xmat <- function(x, y) {
        Xx <- bs(x, knots = xknots, degree = degree, intercept = TRUE)
        Xy <- bs(y, knots = yknots, degree = degree, intercept = TRUE)
        X <- kronecker(Xx, matrix(1, nc = ncol(Xy))) * kronecker(matrix(1, nc = ncol(Xx)), Xy)
        return(X)
    }
    X <- Xmat(x, y)

    if (!is.null(z)) {
        zname <- deparse(substitute(z))
        X <- X * z
    }

    xd <- length(xknots) + degree + 1
    yd <- length(yknots) + degree + 1

    Kx <- diff(diag(xd), differences = differences)
    Kx <- crossprod(Kx, Kx) 
    Ky <- diff(diag(yd), differences = differences)
    Ky <- crossprod(Ky, Ky) 
    K <- kronecker(Kx, diag(yd)) + kronecker(diag(xd), Ky)

    XtX <- crossprod(X, X)
    ### Eigenwerte?
    foo <- function(lambda) 
        (sum(diag(solve(XtX + lambda * K) %*% XtX)) - df)^2
    fooopt <- optimize(foo, interval = c(0, 100000))
    lambda <- fooopt$minimum
    Xsolve <- tcrossprod(solve(XtX + lambda * K), X)
    rm(XtX)
    rm(K)

    fit <- function(u) list(coef = Xsolve %*% u)
    fitted <- function(object) X %*% object$coef
    predict <- function(object, newdata = NULL) {

        if (is.null(newdata)) return(fitted(object))

        newX <- Xmat(newdata[[xname]], newdata[[yname]])
        if (!is.null(z)) newX <- newX * newdata[[zname]]
        newX %*% object$coef
    }
    hatmatrix <- function() X %*% Xsolve

    attr(X, "fit") <- fit
    attr(X, "predict") <- predict
    attr(X, "fitted") <- fitted
    attr(X, "hatmatrix") <- hatmatrix
    attr(X, "lambda") <- lambda
    return(X)
}
