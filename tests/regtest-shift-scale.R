
library("tram")
library("trtf")
options(digits = 4)
set.seed(29)

N <- 1000
x1 <- round(runif(N), 2)
x2 <- round(runif(N), 2)
y <- rnorm(N, mean = c(0,1)[1 + (x1 > .5)], 
              sd = sqrt(exp(c(0,1)[1 + (x2 > .5)])))
d <- data.frame(y = y, x1 = x1, x2 = x2)

### split wrt shift and scale
m0 <- Lm(y ~ 1, data = d)
tr0 <- trafotree(m0, formula = y ~ 1 | x1 + x2, data = d,
                 parm = NULL, intercept = "shift-scale", maxdepth = 2)
logLik(tr0)
cr0 <- info_node(node_party(tr0))$criterion

### split wrt coef(as.mlt(<Lm>)): the same as tr0
tr1 <- trafotree(m0, formula = y ~ 1 | x1 + x2, data = d,
                 maxdepth = 2)
logLik(tr1)
cr1 <- info_node(node_party(tr1))$criterion

### split wrt Intercept only
tr2 <- trafotree(m0, formula = y ~ 1 | x1 + x2, data = d,
                 parm = 1, maxdepth = 2)
logLik(tr2)
cr2 <- info_node(node_party(tr2))$criterion

### ctree: the same as tr2
ct <- ctree(y ~ x1 + x2, data = d)
cr3 <- info_node(node_party(ct))$criterion

all.equal(cr0, cr1)
all.equal(cr2, cr3)

### tr1 via mob
normlm <- function(formula, data = list()) {
  rval <- lm(formula, data = data)
  class(rval) <- c("normlm", "lm")
  return(rval)
}
estfun.normlm <- function(obj) {
  res <- residuals(obj)
  ef <- NextMethod(obj)
  sigma2 <- mean(res^2)
  rval <- cbind(ef, res^2 - sigma2)
  colnames(rval) <- c(colnames(ef), "(Variance)")
  return(rval)
}

normlm_fit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  normlm(y ~ 0 + x, ...)
}

tr3 <- mob(y ~ 1 | x1 + x2, data = d, fit = normlm_fit, 
           control = mob_control(maxdepth = 3))

### tests are different, but nodes should be the same
all.equal(predict(tr3, type = "node"),
          predict(tr1, type = "node"),
          check.attributes = FALSE)
 
cf <- coef(tr0)

tmp <- as.mlt(m0)
apply(cf, 1, function(x) {
    coef(tmp) <- x
    class(tmp) <- class(m0)
    cf <- coef(tmp, as.lm = TRUE)
    c(cf, attr(cf, "scale"))
})
sqrt(exp(1))

### richer baseline function: remove conditional normality assumption
m0 <- BoxCox(y ~ 1, data = d)
tr4 <- trafotree(m0, formula = y ~ 1 | x1 + x2, data = d,
                 parm = NULL, intercept = "shift-scale", maxdepth = 2)
logLik(tr4) > logLik(tr0)

