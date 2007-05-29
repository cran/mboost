
library("mboost")
attach(asNamespace("mboost"))

source("../R/gamboost2.R")
source("../R/base.R")
x <- sort(runif(100, min = 0, max = 2 * pi))
tmp <- data.frame(x = x, y = sin(x) + rnorm(length(x), sd = 0.5))

a <- gamboost(y ~ mybs(x), data = tmp)
plot(y ~ x, data = tmp)
lines(tmp$x, fitted(a))

fitted(a) - predict(a, newdata = tmp)
