
### two possible interpretations of weights:
### 1) case counts: observation i is w_i times in the sample
### 2) relative weights: observation i is given weight w_i
rescale_weights <- function(w) {
    if (max(abs(w - floor(w))) < sqrt(.Machine$double.eps))
        return(w)
    return(w / sum(w) * sum(w > 0))
}

### data preprocessing
boost_dpp <- function(formula, data, weights = NULL, ...) {

    env <- ModelEnvFormula(formula, data)
    y <- env@get("response")
    if (length(y) != 1)
        stop("cannot deal with multivariate response variables")
    y <- y[[1]]
    x <- env@get("designMatrix")

    if (is.null(weights))
        weights <- rep.int(1, NROW(x))

    if (is.factor(y)) {
        if (nlevels(y) != 2)
            stop("not a binary classification problem")
        yfit <- as.numeric(y) - 1
        yfit[yfit == 0] <- -1
    } else {
        yfit <- y
    }

    RET <- gb_xyw(x, y, weights)
    RET$formula <- formula
    RET$menv <- env
    class(RET) <- "boost_data"
    RET
}

gb_xyw <- function(x, y, w) {

    if (is.null(w))
        w <- rep.int(1, NROW(x))

    if (is.factor(y)) {
        if (nlevels(y) != 2)
            stop("not a binary classification problem")
        yfit <- as.numeric(y) - 1
        yfit[yfit == 0] <- -1
    } else {
        yfit <- y
    }

    list(x = x, y = y, yfit = yfit, w = w)
}

### check for negative gradient corresponding to L2 loss
checkL2 <- function(object)
    isTRUE(all.equal(attributes(object$family)[-c(4:8)], 
                     attributes(GaussReg())[-c(4:8)]))

gm <- function(object) UseMethod("gm")

gm.glmboost <- function(object) {

    mstop <- nrow(object$ensemble)
    x <- object$data$x
    RET <- matrix(0, nrow = NROW(x), ncol = mstop)

    jsel <- object$ensemble[,"xselect"]
    cf <- object$ensemble[,"coef"] * object$control$nu 

    for (m in 1:mstop)
        RET[,m] <- cf[m] * x[,jsel[m]]

    RET[,1] <- RET[,1] + object$offset

    return(RET)
}

gm.gamboost <- function(object) {

    mstop <- nrow(object$ensemble)
    x <- object$data$x
    RET <- matrix(0, nrow = NROW(x), ncol = mstop)
    nu <- object$control$nu

    jsel <- object$ensemble[,"xselect"]

    for (m in 1:mstop)
        RET[,m] <- nu * predict(object$ensembless[[m]], x[,jsel[m]])$y

    RET[,1] <- RET[,1] + object$offset

    return(RET) 
}

### partial fits
gamplot <- function(object) {

     x <- object$data$x
     lp <- matrix(0, ncol = NCOL(x), nrow = NROW(x))
     ens <- object$ensemble
     ensss <- object$ensembless
     nu <- object$control$nu
     mstop <- nrow(ens)
     for (m in 1:mstop) {
         xselect <- ens[m,"xselect"]
         lp[,xselect] <- lp[,xselect] + nu * predict(ensss[[m]], 
                                     x = x[,xselect])$y
     }
     colnames(lp) <- colnames(x)
     lp
}

bhatmat <- function(n, H, xselect, fitm, fW) {
    
    B <- matrix(0, nrow = n, ncol = n)
    I <- diag(n)
    tr <- numeric(length(xselect))

    for (m in 1:length(xselect)) {
        B <- B + (H[[xselect[m]]] * fW(fitm[,m])) %*% (I - B)
        tr[m] <- sum(diag(B))
    }
    list(hatmatrix = B, trace = tr)
}

### fractional polynomials transformation
### all powers `p' of `x', all powers `p' of `x' times `log(x)' and `log(x)'
### see Sauerbrei & Royston (1999), JRSS A (162), 71--94
FP <- function(x, p = c(-2, -1, -0.5, 0.5, 1, 2, 3)) {
    xname <- deparse(substitute(x))
    ### <FIXME> do we need this? 
    ### map x into [1, 2]
    ##if (scale) {
    ##    x <- (x - mean(x))
    ##    x <- x / (max(abs(x)) * 2)
    ##    x <- x - min(x) + 1
    ##}
    ### </FIXME>
    if (any(x < sqrt(.Machine$double.eps)))
        stop("negative values in ", sQuote(xname), "\n")
    Xp <- sapply(p, function(p) x^p)
    Xp <- cbind(Xp, log(x))
    Xp <- cbind(Xp, Xp * log(x))   
    colnames(Xp) <- c(paste(xname, "^", p, sep = ""),
                      paste("log(", xname, ")", sep = ""),
                      paste("log(", xname, ")", xname, "^", p, sep = ""),
                      paste("log(", xname, ")^2", sep = ""))
    Xp
}

### inverse probability of censoring weights
### see van der Laan & Robins (2003)
IPCweights <- function(x, maxweight = 5) {

    if (!extends(class(x), "Surv"))
        stop(sQuote("x"), " is not a Surv object")

    event <- x[,"status"]
    x[,"status"] <- 1 - event
    km <- survfit(x)   
    Ghat <- getsurv(km, time = x[,"time"])
    Ghat[event == 0] <- 1
    w <- event / Ghat
    w[w > maxweight] <- maxweight
    w
}

### extract survival probabilities
### taken from ipred:::getsurv
### DO NOT TOUCH HERE
getsurv <- function(obj, times)
{
    # get the survival probability for times from KM curve j'

    if (!inherits(obj, "survfit")) stop("obj is not of class survfit")
    # <FIXME: methods may have problems with that>
    class(obj) <- NULL
    # </FIXME>
    lt <- length(times)
    nsurv <- times

    # if the times are the same, return the km-curve

    if(length(times) == length(obj$time)) {
        if (all(times == obj$time)) return(obj$surv)
    }

    # otherwise get the km-value for every element of times separatly

    inside <- times %in% obj$time
    for (i in (1:lt)) {
        if (inside[i])
            nsurv[i] <- obj$surv[obj$time == times[i]]
        else  {
            less <- obj$time[obj$time < times[i]]
            if (length(less) == 0)
                nsurv[i] <- 1
            else
                nsurv[i] <- obj$surv[obj$time == max(less)]
        }
    }
    nsurv
}
