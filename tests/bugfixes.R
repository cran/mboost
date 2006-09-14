
require("mboost")

set.seed(290875)

### predict did not return factor levels for blackboost models
dummy <- data.frame(y = gl(2, 100), x = runif(200))
pr <- predict(blackboost(y ~ x, data = dummy, family = Binomial()), 
              newdata = dummy, type = "response")
stopifnot(is.factor(pr) && all(levels(pr) %in% levels(dummy$y)))

### predict for g{al}mboost.matrix did not work
ctrl <- boost_control(mstop = 10)
gb <- glmboost(x = cbind(1, dummy$x), y = dummy$y, family = Binomial(), 
               control = ctrl)
stopifnot(all.equal(predict(gb), predict(gb, newdata = cbind(1, dummy$x))))

gb <- gamboost(x = cbind(1, dummy$x), y = dummy$y, family = Binomial(), 
               control = ctrl)
stopifnot(all.equal(predict(gb), predict(gb, newdata = cbind(1, dummy$x))))

