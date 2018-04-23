
library("trtf")
library("survival")
library("partykit")
data("GBSG2", package = "TH.data")
set.seed(290875)
### Make UL and Windooze happy
options(digits = 5)

yvar <- numeric_var("y", support = c(100, 2000), bounds = c(0, Inf))
By <- Bernstein_basis(yvar, order = 5, ui = "incre")
m <- ctm(response = By, todistr = "MinExt")
GBSG2$y <- with(GBSG2, Surv(time, cens))

tf <- traforest(m, formula = y ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2, 
    control = ctree_control(splitstat = "quad", teststat = "quad",
                    testtype = "Teststatistic", mincriterion = 1, minbucket = 50), 
    ntree = 12, trace = FALSE, cores = 1)

predict(tf, newdata = GBSG2[1:3,], type = "weights")
predict(tf, newdata = GBSG2[1:3,], type = "node")
(cf <- predict(tf, newdata = GBSG2[1:3,], type = "coef"))

logLik(tf, newdata = GBSG2[1:3,])

sapply(1:length(tf$nodes), function(i) logLik(gettree(tf, i)))

varimp(tf)

mod <- mlt(m, data = GBSG2)
logLik(mod)

logLik(mod, newdata = GBSG2[1:3,])

layout(matrix(1:2, nr = 1))
coef(mod) <- cf[[1]]
plot(mod, newdata = data.frame(1), type = "survivor")
coef(mod) <- cf[[2]]
plot(mod, newdata = data.frame(1), type = "survivor", add = TRUE)
coef(mod) <- cf[[3]]
plot(mod, newdata = data.frame(1), type = "survivor", add = TRUE)

simulate(tf, newdata = GBSG2[1:3,], nsim = 3)

cmod <- coxph(Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade +
    pnodes + progrec + estrec, data = GBSG2)

plot(survfit(cmod, newdata = GBSG2[1:3,]))
