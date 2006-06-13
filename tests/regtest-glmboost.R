
require("mboost")

set.seed(290875)

dgp <- function(n = 100, beta = rep(0, 10), sd = 1) {
    p <- length(beta) - 1
    x <- cbind(1, matrix(runif(n * p), ncol = p))
    lp <- x %*% beta
    y <- lp + rnorm(n, sd = sd)
    ls <- data.frame(y = y, x[,-1])
    attr(ls, "lp") <- lp
    ls
}

### well-defined problem
mydf <- dgp(beta = c(1, 2.5, rep(0, 2)))

### for easy comparison with lm
fm <- GaussReg()
fm@offset <- function(y, w) 0

mydf.gb <- glmboost(y ~ ., data = mydf, family = fm, 
                    control = boost_control(mstop = 10000))
mydf.lm <- lm(y ~ ., data = mydf)

### compare coefficients
stopifnot(max(abs(coef(mydf.gb) - coef(mydf.lm))) < 1e-10)

### a little bit more difficult
mydf <- dgp(beta = c(1, 2.5, rep(0, 38)))

mydf.gb <- glmboost(y ~ ., data = mydf, family = fm, 
                    control = boost_control(mstop = 10000))
aic <- AIC(mydf.gb, method = "corrected")
mstop(aic)
mydf.lm <- lm(y ~ ., data = mydf)

### compare coefficients
which(abs(coef(mydf.lm)) < abs(coef(mydf.gb[mstop(aic)])))

#### check boosting hat matrix and subsetting / predict
stopifnot(isTRUE(all.equal(drop(attr(aic, "hatmat") %*% mydf$y),
                           as.vector(predict(mydf.gb)))))
stopifnot(isTRUE(all.equal(drop(attr(AIC(mydf.gb[255]), "hatmat") %*% mydf$y),
                           as.vector(predict(mydf.gb[255])))))
stopifnot(isTRUE(all.equal(drop(attr(AIC(mydf.gb[255]), "hatmat") %*% mydf$y),
                           as.vector(fitted(mydf.gb[255])))))

### a simple two-dimensional example from `glmboost.Rd'
data("cars")
cars.gb <- glmboost(dist ~ speed, data = cars, family = fm, 
                    control = boost_control(mstop = 5000))
cars.gb

### coefficients should coincide
coef(cars.gb)
coef(lm(dist ~ speed, data = cars))

### logistic regression
mydf <- data.frame(x = runif(100), y = gl(2, 50))
bmod <- glmboost(y ~ x, data = mydf, family = Binomial(), 
                 control = boost_control(mstop = 10000))
gmod <- glm(y ~ x, data = mydf, family = binomial())
llg <- logLik(gmod)
attributes(llg) <- NULL
stopifnot(all.equal(logLik(bmod), llg))
stopifnot(max(abs(predict(gmod, type = "link")/2 - fitted(bmod))) < 
                  sqrt(.Machine$double.eps))
stopifnot(all.equal(coef(bmod) * 2, coef(gmod)))

### weighted least squares problem

x <- runif(100)  
df <- data.frame(y = 2 + 3 * x + rnorm(length(x)),
                 x = x, z = runif(length(x)),
                 w = runif(length(x)) * 10)

### linear model, classical fit
lmmod <- lm(y ~ x + z, data = df, weights = w)

### linear model, boosting fit
lmb <- glmboost(y ~ x + z, data = df, weights = df$w,
                control = boost_control(mstop = 20000))

### compare fitted values
stopifnot(max(abs(fitted(lmmod) -fitted(lmb))) < sqrt(.Machine$double.eps))

### compare hat matrices
stopifnot(max(abs(hatvalues(lmmod) - hatvalues(lmb))) < sqrt(.Machine$double.eps))      

### compare boosting hat matrix with fitted values
stopifnot(max(abs(attr(AIC(lmb), "hatmatrix") %*% (df$y - lmb$offset) + lmb$offset -
        fitted(lmb))) < sqrt(.Machine$double.eps)) 

