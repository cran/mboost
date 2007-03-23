
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

### blackboost _did_ touch the response, arg!

data("bodyfat", package = "mboost")
ctrl <- boost_control(mstop = 500, nu = 0.01)
bb <- blackboost(DEXfat ~ ., data = bodyfat, control = ctrl)
n <- nrow(bodyfat)
bs <- rmultinom(3, n, rep(1, n) / n)
x <- seq(from = 10, to = 500, by = 10)
cv <- cvrisk(bb, bs, grid = x)
ctrl$risk <- "oobag"
tmp <- blackboost(DEXfat ~ ., data = bodyfat, control = ctrl,
                 weights = bs[,3])

stopifnot(identical(max(abs(tmp$risk[x] / sum(bs[,3] == 0)  - cv[3,])), 0))

### center = TRUE and cvrisk were broken; same issue with masking original data

gb <- glmboost(DEXfat ~ ., data = bodyfat, control = boost_control(center = TRUE))
cv1 <- cvrisk(gb, folds = bs)
tmp <- glmboost(DEXfat ~ ., data = bodyfat, 
                control = boost_control(center = TRUE, risk = "oobag"), 
                weights = bs[,3])
stopifnot(identical(max(tmp$risk[attr(cv1, "mstop")] / sum(bs[,3] == 0) - cv1[3,]), 0))

### same problem, just another check

indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
cbodyfat <- bodyfat
cbodyfat[indep] <- lapply(cbodyfat[indep], function(x) x - mean(x))
bffm <- DEXfat ~ age + waistcirc + hipcirc + elbowbreadth + kneebreadth +
      anthro3a + anthro3b + anthro3c + anthro4

bf_glm_1 <- glmboost(bffm, data = cbodyfat)
cv1 <- cvrisk(bf_glm_1, folds = bs)
bf_glm_2 <- glmboost(bffm, data = bodyfat, control = boost_control(center = TRUE))
cv2 <- cvrisk(bf_glm_2, folds = bs)

stopifnot(mstop(cv1) == mstop(cv2))

### dfbase=1 was not working correctly for ssp
### spotted by Matthias Schmid <Matthias.Schmid@imbe.imed.uni-erlangen.de>
data("bodyfat", package = "mboost")
ctrl <- boost_control(mstop = 100, center = TRUE)
ga <- gamboost(DEXfat ~ ., data = bodyfat, dfbase = 1, control = ctrl)
gl <- glmboost(DEXfat ~ ., data = bodyfat, control = ctrl)
stopifnot(max(abs(predict(ga) - predict(gl))) < 1e-8)
AIC(gl)

