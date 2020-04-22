
rm(list=ls())

source("competitors.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions used for the data generating process -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Friedman function
f1 <- function(x1, x2, x3, x4, x5, scale = TRUE) {
  ret <- 10 * sin(pi * x1 * x2) + 20 * (x3 - 0.5)^2 + 10 * x4 + 5 * x5
  if (scale) {
    ret <- ret - min(ret) 
    ret <- ret / max(ret)
    ret <- 5 * ret - 2.5 # OR
  }
  ret
}

### model scale function
## if there is no scale in the model, prod_scale = 0
scale_f <- function(d, prod_scale = 1, ...){
  ### roughly exp(-1.5, 1.5)
  exp(with(d, f1(X1, X2, X3, X4, X5, scale = TRUE) / (2.5 / 1.5)) * prod_scale)}

### model shift function
## if there is no shift in the model, prod_shift = 0
shift_f <- function(d, prod_shift = 1, ...){
  with(d, f1(X6, X7, X8, X9, X10, scale = TRUE)) * prod_shift}

rmodel <- function(n, shift, scale) {
  theta <- matrix(logit(c(.15, .5, .85)), nrow = n, ncol = 3, byrow = TRUE)
  lp <- scale * theta + shift
  cdf <- plogis(lp)
  prb <- cbind(cdf[,1], cdf[,2] - cdf[,1], cdf[,3] - cdf[,2], 1 - cdf[,3])
  y <- apply(prb, 1, function(r) which(rmultinom(1, size = 1, prob = r) > 0))
  y <- factor(y, levels = 1:4, labels = 1:4, ordered = TRUE)
  return(list(y = y, prb = prb))
}

ry <- function(d, prod_shift, prod_scale, ...){
  rmodel(n = nrow(d),
         shift = shift_f(d, prod_shift = prod_shift),
         scale = scale_f(d, prod_scale = prod_scale))
}

### data generation
dgp <- function(ntrain, ntest,
                p = 5,
                prod_shift = 0,
                prod_scale = 0, ...) {
  n <- ntrain + ntest
  dat <- data.frame(matrix(runif(10 * n), ncol = 10),
                  matrix(runif(n * p, min = 0, max = 1), nrow = n))
  ry_tmp <- ry(dat, prod_shift = prod_shift,
               prod_scale = prod_scale)
  dat$y <- ry_tmp$y
  prb <- ry_tmp$prb
  return(list('train' = dat[1:ntrain,], 'test' = dat[(ntrain+1):n,],
              'test_prb' = prb[(ntrain+1):n,],
              prod_shift = prod_shift, prod_scale = prod_scale))
}

### log-likelihood of the model
mylogLik <- function(d){
  shift <- shift_f(d$test, prod_shift = d$prod_shift)
  scale <- scale_f(d$test, prod_scale = d$prod_scale)
  FZ <- mlt:::.Logistic()
  # theta von rmodel(...)
  h <- c(-Inf, FZ$q(c(.15, .5, .85)), Inf)
  hu <- h[unclass(d$test$y) + 1]
  hl <- h[unclass(d$test$y)]
  sum(log(FZ$p(hu*scale + shift) - FZ$p(hl*scale + shift) ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data generating process -------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(12345)

if(!dir.exists("rda")){dir.create("rda")}

for(nsim_i in seq_along(unique(sim_para[, "nsim"]))){
  
  # nsim_i=1
  # nsim_i=2
  
  ## number of repetitions
  nsim <- unique(sim_para[, "nsim"])[nsim_i]
  
  ## size of training sample
  ntrain <- 250
  
  ## size of validation sample
  ntest <- 500
  
  ## scale and shape parameters of Weibull distribution
  wscale <- 1
  wshape <- 1
  
  # data simulation process -------------------------------------------------
  
  args <- expand.grid(p = c(5, 50), 
                      prod_shift = c(0, 1), 
                      prod_scale = c(0, 1))
  
  simdat <- vector(mode = "list", length = nrow(args))
  loglik <- vector(mode = "list", length = nrow(args))
  
  for (d in 1:nrow(args)) {
    
    # d = 8
    
    simdat[[d]] <- replicate(unique(sim_para[, "nsim"])[nsim_i],
                             dgp(ntrain = ntrain,
                                 ntest = ntest,
                                 p = args[d, "p"],
                                 prod_shift = args[d, "prod_shift"],
                                 prod_scale = args[d, "prod_scale"]),
                             simplify = FALSE)
    loglik[[d]] <- lapply(simdat[[d]], mylogLik)
  }
  
  save(args, simdat, loglik,
       file = paste("rda/simulated_data_nsim",
                    unique(sim_para[, "nsim"])[nsim_i],
                    ".rda", sep = ""))
  
}

