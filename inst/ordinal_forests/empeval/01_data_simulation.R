
# Define global parameters ------------------------------------------------

set.seed(12345)

## number of repetitions
nsim <- 100

## size of learning sample
ntrain <- 250

## size of validation sample
ntest <- 500

## scale and shape parameters of Weibull distribution
wscale <- 1
wshape <- 1

### Friedman function
f1 <- function(x1, x2, x3, x4, x5, scale = TRUE) {
  ret <- 10 * sin(pi * x1 * x2) + 20 * (x3 - 0.5)^2 + 10 * x4 + 5 * x5
  if (scale) {
    ret <- ret - min(ret) 
    ret <- ret / max(ret)
    ret <- 3 * ret - 1.5
  }
  ret
}

### model scale function
## if there is no scale in the model, prod_scale = 0
scale_f <- function(d, prod_scale = 1, ...){
  exp(with(d, f1(X1, X2, X3, X4, X5, scale = TRUE)) * prod_scale)}

### model shift function
## if there is no shift in the model, prod_shift = 0
shift_f <- function(d, prod_shift = 1, ...){
  with(d, f1(X6, X7, X8, X9, X10, scale = TRUE)) * prod_shift}

rmodel <- function(n = 100, shift = 0, scale = 1) {
  # P(Y \eq y \mid x) = F_Z(\epsilon(x) h(y) + \beta(x))
  #                   = expit(\epsilon(x) h(y) + \beta(x))
  FZ <- mlt:::.Logistic()
  q0 <- function(u) qweibull(u, scale = wscale, shape = wshape)
  # h(y) = (FZ$q(u) - shift) / scale
  # expit = h^(-1) ???
  # q1 <- qweibull( expit( h(y))) 
  q1 <- function(u, shift, scale) q0(FZ$p((FZ$q(u) - shift) / scale))
  y <- q1(runif(n), shift = shift, scale = scale)
  
  ### to fix infinite generated values
  inf_idx <- is.infinite(y)
  ninf <- sum(inf_idx)
  while (ninf > 0){
    y[inf_idx] <- q1(runif(ninf), shift = shift[inf_idx], 
                     scale = scale[inf_idx]) 
    inf_idx <- is.infinite(y) 
    ninf <- sum(inf_idx)
  }
  y <- pmax(.Machine$double.eps, y)
  cut(y, breaks = c(0, q0(c(.25, .5, .75)), Inf), ordered_result = TRUE)
}

ry <- function(d, prod_shift, prod_scale, ...)
  rmodel(nrow(d), shift_f(d, prod_shift=prod_shift), 
         scale_f(d, prod_scale=prod_scale))

### data generation
dgp <- function(ntrain, ntest,
                p = 5,
                prod_shift = 0,
                prod_scale = 0, ...) {
  n <- ntrain + ntest
  d <- data.frame(matrix(runif(10 * n), ncol = 10),
                  matrix(runif(n * p, min = 0, max = 1), nrow = n))
  d$y <- ry(d, prod_shift = prod_shift, prod_scale = prod_scale)
  list('train' = d[1:ntrain,], 'test' = d[(ntrain+1):n,],
       prod_shift = prod_shift, prod_scale = prod_scale)
}

### log-likelihood of the model
mylogLik <- function(d){
  shift <- shift_f(d$test, prod_shift = d$prod_shift)
  scale <- scale_f(d$test, prod_scale = d$prod_scale)
  FZ <- mlt:::.Logistic()
  h <- c(-Inf, FZ$q(c(.25, .5, .75)), Inf)
  hu <- h[unclass(d$test$y) + 1]
  hl <- h[unclass(d$test$y)]
  sum(log(FZ$p(hu*scale + shift) - FZ$p(hl*scale + shift) ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data simulation process -------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

args <- expand.grid(p = c(5, 50), 
                    prod_shift = c(0, 1), 
                    prod_scale = c(0, 1))

simdat <- vector(mode = "list", length = nrow(args))
loglik <- vector(mode = "list", length = nrow(args))

for (i in 1:nrow(args)) {

    simdat[[i]] <- replicate(nsim, dgp(ntrain = ntrain,
                                       ntest = ntest,
                                       p = args[i, "p"],
                                       prod_shift = args[i, "prod_shift"],
                                       prod_scale = args[i, "prod_scale"]),
                             simplify = FALSE)
    loglik[[i]] <- lapply(simdat[[i]], mylogLik)
}

save(args, simdat, loglik, file='simulated_data.rda')
