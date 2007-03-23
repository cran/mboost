
hatMatTH <- function(x, w = NULL, df = 5) {
    n <- NROW(x)
    indx <- diag(n)
    x <- signif(x, 10)
    apply(indx, 2, function(y)
        predict(smoothbase(x = x, ux = unique(sort(x)), y = y, w = w, df = df),
                x = x)$y)
}

ssp <- function(df = 4) {

    dpp <- function(x, weights) {
        xs <- signif(x, 10)
        ux <- unique(sort(xs))

        fit <- function(y) {
            ret <- smoothbase(x = xs, ux = ux, y = y, w = weights, df = df)
            if (df == 1) {
                ret$hatvalues <- function() {
                    if (length(ret$coef) == 1) return(xs %*% t(xs) / sum(xs^2))
                    X <- cbind(1, xs)
                    X %*% solve(t(X) %*% X) %*% t(X)
                }
                return(ret)
            }
            ret$hatvalues <- function()
                hatMatTH(x = x, w = weights, df = df)
            class(ret) <- "ssp"
            ret
        }
        fit
    }
    return(dpp)
}

fitted.ssp <- function(object, ...)
    object$yfit

predict.ssp <- function(object, newdata, ...)
    stats:::predict.smooth.spline.fit(object, x = newdata, ...)$y

hatvalues.ssp <- function(model, ...)
    model$hatvalues()

hatvalues.lmfit <- function(model, ...)
    model$hatvalues()

bsp <- function(df = 4) {

    dpp <- function(x, weights) {
        if (df > 1) {
            bx <- bs(x, df = df)
        } else {
            bx <- matrix(x, nrow = NROW(x))
            if (length(unique(as.vector(x))) > 5) bx <- cbind(1, bx)
        }

        xtx <- solve(t(bx) %*% diag(weights) %*% bx) %*% t(bx) %*% diag(weights)

        fit <- function(y) {
            ret <- list(coef = xtx %*% y)
            ret$fitted.values <- drop(bx %*% ret$coef)
            ret$predict <- function(newx) {
                if (df > 1) return(drop(bs(newx, df = df, Bound = range(x)) %*% ret$coef))
                return(drop(matrix(newx, ncol = ncol(bx)) %*% ret$coef))
            }
            ret$hatvalues <- function() bx %*% xtx
            class(ret) <- "bsp"
            ret
        }
        fit
    }
    return(dpp)
}

fitted.bsp <- function(object, ...)
    object$fitted.values

predict.bsp <- function(object, newdata, ...)
    object$predict(newdata)

hatvalues.bsp <- function(model, ...)
    model$hatvalues()

ols <- function(df = 1) bsp(df = 1)
