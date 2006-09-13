
require("mboost")

set.seed(290875)

dummy <- data.frame(y = gl(2, 100), x = runif(200))
pr <- predict(blackboost(y ~ x, data = dummy, family = Binomial()), 
              newdata = dummy, type = "response")
stopifnot(is.factor(pr) && all(levels(pr) %in% levels(dummy$y)))
