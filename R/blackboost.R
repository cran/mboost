
### 
### Experimental version of gradient boosting with conditional trees as 
### as base learner
### 


### Fitting function
blackboost_fit <- function(object, tree_controls, fitmem, family = GaussReg(), 
                        control = boost_control(), weights = NULL) {

    ### number of observations in the learning sample
    ### make sure this gets _copied_
    y <- .Call("copymem", party:::get_variables(object@responses)[[1]], 
               package = "mboost")
    if (is.factor(y)) {
        y <- (2 * (as.numeric(y) - 1)) - 1
        object@responses <- party:::initVariableFrame(data.frame(y = y), NULL)
    }
    if (is.null(weights)) weights <- object@weights
    storage.mode(weights) <- "double"

    ### hyper parameters
    mstop <- control$mstop
    constraint <- control$constraint
    nu <- control$nu

    ### the ensemble
    ens <- vector(mode = "list", length = mstop)

    ### extract negative gradient function
    ngradient <- family@ngradient
    if (!family@weights && any(max(abs(weights - 1))))
        stop(sQuote("family"), " is not able to deal with weights")

    fit <- offset <- family@offset(y, weights)
    u <- ustart <- ngradient(y, fit)

    where <- rep(1, object@nobs)
    storage.mode(where) <- "integer"

    ### start boosting iteration
    for (m in 1:mstop) {
  
        ### fit tree to residuals
        .Call("R_modify_response", as.double(u), object@responses, 
              PACKAGE = "party")
        ens[[m]] <- .Call("R_TreeGrow", object, weights, fitmem, tree_controls,
                          where, PACKAGE = "party")

        ### check if first node is terminal, i.e., if at least 
        ### one split was performed
        if (ens[[m]][[4]])
            warning("could not split root node in iteration ", m, 
                    ", decrease ", sQuote("mincriterion"))

        ### update step
        fit <- fit + nu * unlist(.Call("R_getpredictions", ens[[m]], where, 
                                PACKAGE = "party"))

        ### L2 boost with constraints (binary classification)
        if (constraint)
            fit <- sign(fit) * pmin(abs(fit), 1)

        ### negative gradient vector, the new `residuals'
        u <- ngradient(y, fit)   

    }

    RET <- list(ensemble = ens, control = control, ustart = ustart,
                fit = fit, family = family, offset = offset)
    if (control$savedata) RET$data <- object
    class(RET) <- "blackboost"

    ### prediction function
    RET$predict <- function(newdata = NULL, mstop = mstop, ...) {

        if (is.null(newdata)) {
            newinp <- object@inputs
        } else {
            newinp <- object@menv@get("input", data = newdata)
            newinp <- party:::initVariableFrame(newinp, trafo = NULL)
        }

        p <- offset
        for (m in 1:mstop) {
            wh <- .Call("R_get_nodeID", RET$ensemble[[m]], newinp, 0.0, 
                        PACKAGE = "party")
            p <- p + nu * unlist(.Call("R_getpredictions", 
                 RET$ensemble[[m]], wh, PACKAGE = "party"))
        }
        if (constraint) p <- sign(p) * pmin(abs(p), 1)
        return(p)
    }

    return(RET)
}

### methods: subset
"[.blackboost" <- function(x, i, ...) { 
    if (i == length(x$ensemble)) return(x)
    if (length(i) != 1)
        stop("not a positive integer")
    if (i < 1 || i > length(x$ensemble))
        warning("invalid number of boosting iterations")
    indx <- 1:min(max(i, 1), nrow(x$ensemble))
    x$ensemble <- x$ensemble[indx]
    x$fit <- x$predict(mstop = max(indx))
    x
}

### methods: prediction
predict.blackboost <- function(object, newdata = NULL, 
                              type = c("lp", "response"), ...) {
    y <- object$data$y
    type <- match.arg(type)
    lp <- object$predict(newdata = newdata, mstop = length(object$ensemble), ...)
    if (type == "response" && is.factor(y))
        return(factor(levels(y)[(lp > 0) + 1], levels = levels(y)))
    return(lp)
}

blackboost <- function(formula, data = list(), weights = NULL, 
                      tree_controls = ctree_control(teststat = "max",
                          testtype = "Teststatistic",
                          mincriterion = 0,
                          maxdepth = 2), ...) {

    ### construct design matrix etc.
    object <- party:::ctreedpp(formula, data, ...)
    fitmem <- ctree_memory(object, TRUE)

    ### fit the ensemble
    RET <- blackboost_fit(object, tree_controls = tree_controls, 
                         fitmem = fitmem, weights = weights, ...)

    RET$call <- match.call()

    return(RET)
}

print.blackboost <- function(x, ...) {

    cat("\n")
    cat("\t Tree-Based Gradient Boosting\n")
    cat("\n")
    if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
    show(x$family)
    cat("\n")
    cat("Number of boosting iterations: mstop =", length(x$ensemble), "\n")
    cat("Step size: ", x$control$nu, "\n")
    cat("\n")
    invisible(x)
}
