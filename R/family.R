
setClass("boost_family", representation = representation(
    ngradient  = "function",
    loss       = "function",
    risk       = "function",
    offset     = "function",
    fW         = "function",
    weights    = "logical",
    name       = "character",
    charloss   = "character"
))

setMethod("show", "boost_family", function(object) {
    cat("\n\t", object@name, "\n\n")
    cat("Loss function:", object@charloss, "\n")
})

Family <- function(ngradient, loss = NULL, risk = NULL, 
                   offset = function(y, w) 0, 
                   fW = function(f) rep(1, length(f)),
                   weights = TRUE, name = "user-specified") {

    if (is.null(loss))
        loss <- function(y, f) NA
    if (is.null(risk))
        risk <- function(y, f, w = 1) sum(w * loss(y, f))
    RET <- new("boost_family", ngradient = ngradient, loss = loss, 
               risk = risk, offset = offset, fW = fW, weights = weights, 
               name = name, charloss = paste(deparse(body(loss)), "\n"))
    RET
}

### Gaussian (Regression)
GaussReg <- function()
    Family(ngradient = function(y, f) y - f,
           loss = function(y, f) (y - f)^2,
           offset = weighted.mean,
           name = "Squared Error (Regression)")

### Gaussian (-1 / 1 Binary Classification)
GaussClass <- function()
    Family(ngradient = function(y, f) - 2 * y + 2 * y * f,
           loss = function(y, f) 1 - 2 * y * f + (y * f)^2,
           name = "Squared Error (Classification)")

### Laplace
Laplace <- function()
    Family(ngradient = function(y, f) sign(y - f),
           loss = function(y, f) abs(y - f),
           offset = function(y, w) median(y),
           name = "Absolute Error")

### Binomial
Binomial <- function()
    Family(ngradient = function(y, f) {
               exp2yf <- exp(-2 * y * f)
               -(-2 * y * exp2yf) / (log(2) * (1 + exp2yf))
           },
           loss = function(y, f) {
               p <- exp(f) / (exp(f) + exp(-f))
               y <- (y + 1) / 2
               -y * log(p) - (1 - y) * log(1 - p)
           },
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               1/2 * log(p / (1 - p))
           },
           fW = function(f) {
               p <- exp(f) / (exp(f) + exp(-f))
               4 * p * (1 - p)
           },
           name = "Negative Binomial Likelihood")

### Poisson
Poisson <- function()
    Family(ngradient = function(y, f) y - exp(f),
           loss = function(y, f) -y*f + exp(f),
           offset = weighted.mean,
           name = "Poisson")

### L1Huber
Huber <- function(d = NULL) {
    mc <- match.call()
    if (length(mc) > 2)
        dtxt <- deparse(mc[[2]])
    else
        dtxt <- NULL
    fit <- 0
    Family(ngradient = function(y, f) {
               if (is.null(d)) d <- median(abs(y - fit))
               fit <<- f
               ifelse(abs(y - f) < d, y - f, d * sign(y - f))
           },
           loss = function(y, f)
               ifelse((a <- abs(y - f)) < d, a^2/2, d*(a - d/2)),
           offset = function(y, w) median(y),
           name = paste("Huber Absolute Error", 
               ifelse(is.null(d), "(with adaptive d)", 
                                  paste("(with d = ", dtxt, ")", sep = ""))))
}

### Adaboost
AdaExp <- function()
    Family(ngradient = function(y, f) y * exp(-y * f),
           loss = function(y, f) exp(-y * f),
           offset = function(y, w) {
               p <- weighted.mean(y > 0, w)
               1/2 * log(p / (1 - p))
           },
           name = "Adaboost Exponential Error")
