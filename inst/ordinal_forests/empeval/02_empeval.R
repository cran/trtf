
start_time <- Sys.time()

library("parallel")
RNGkind("L'Ecuyer-CMRG")

ncores <- detectCores()

set.seed(12345)

source("competitors.R")
load("simulated_data.rda")

# Define global parameters ------------------------------------------------

### joint parameters for all forest procedures

### number of repetitions
nsim <- length(simdat[[1]])

## number of trees
ntree <- 250

## maximal depth of trees
nodedepth <- 10

## number of obs in terminal node
minbucket <- 20

### parameters for party/partykit
ctrl_partykit <- partykit:::ctree_control(teststat = "quad", 
                                          testtype = "Univ", 
                                          mincriterion = 0,
                                          minbucket = minbucket,
                                          maxdepth = nodedepth)

ret <- args[rep(1:nrow(args), each = nsim),]
simdat <- do.call("c", simdat)
ret$tf_B <- ret$tf_B_alpha <- ret$cf <- ret$of <- 0
ret$loglik <- unlist(do.call("c", loglik))

### filename to save estimated likelihoods
ret_filename <- paste("results_empeval_", ntree, ".rda", sep = "")

# Ordinal transfromation forest with Bernstein basis and general score --------
# Bs(theta)

trtf_fun <- function(i) {
    ## number of randomly preseleted variables to perform a split
    ## for low-dimensional data it is set to the number of variables
    ## for high-dimensional data it is set to the square root of the number of variables
    mtry <- ifelse(ret[i, "p"] <= 5, 
                   ncol(simdat[[i]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[i]]$train) - 1)))
    
    ## learning and validation
    tf <- mytraforest(simdat[[i]]$train, simdat[[i]]$test, mtry = mtry)
    unclass(tf$logLik_NN)
}

tmp <- runif(10)

### learn and validate for each repetition of the simulated data
tf <- mclapply(1:length(simdat), trtf_fun, mc.cores = ncores)
ret$tf_B <- unlist(tf)
print('Traforest B is done!')
save(ret, file = ret_filename)

# Ordinal transformation forest with Bernstein basis and log-rank score ----------
# Bs(alpha)

trtf_fun_B <- function(i) {
    mtry <- ifelse(ret[i, "p"] <= 5, 
                   ncol(simdat[[i]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[i]]$train) - 1)))
    ## learning and validation
    tf <- mytraforest_B(simdat[[i]]$train, simdat[[i]]$test, mtry)
    unclass(tf$logLik_NN)
}

tmp <- runif(10)

### learn and validate for each repetition of the simulated data
tf <- mclapply(1:length(simdat), trtf_fun_B, mc.cores = ncores)
ret$tf_B_alpha <- unlist(tf)
print('Traforest B_alpha is done!')
save(ret, file = ret_filename)

# Ordinal Forest (equal) --------------------------------------------------

### perffunction = "equal"

library("ranger")
library("ordinalForest")
rg_fun <- function(i) {
    mtry <- ifelse(ret[i, "p"] <= 5, 
                   ncol(simdat[[i]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[i]]$train) - 1)))
    ## learning and validation
    rg <- myordforest(simdat[[i]]$train, simdat[[i]]$test, mtry)
    unclass(rg$logLik_NN)
}

tmp <- runif(10)

### learn and validate for each repetition of the simulated data
rg <- mclapply(1:length(simdat), rg_fun, mc.cores = ncores)
ret$of <- unlist(rg)
print('ordinalForest (equal) is done!')
save(ret, file = ret_filename)

# Ordinal Forest (proportional) -------------------------------------------

### perffunction = "proportional"

library("ranger")
library("ordinalForest")
rg_fun_prop <- function(i) {
  mtry <- ifelse(ret[i, "p"] <= 5, 
                 ncol(simdat[[i]]$train) - 1, 
                 ceiling(sqrt(ncol(simdat[[i]]$train) - 1)))
  ## learning and validation
  rg <- myordforest_prop(simdat[[i]]$train, simdat[[i]]$test, mtry)
  unclass(rg$logLik_NN)
}

tmp <- runif(10)

### learn and validate for each repetition of the simulated data
rg <- mclapply(1:length(simdat), rg_fun_prop, mc.cores = ncores)
ret$of_prop <- unlist(rg)
print('ordinalForest (proportional) is done!')
save(ret, file = ret_filename)

# Conditional Inference Forests -------------------------------------------

cf_fun <- function(i) {
    mtry <- ifelse(ret[i, "p"] <= 5, 
                   ncol(simdat[[i]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[i]]$train) - 1)))
    cf <- mycforest(simdat[[i]]$train, simdat[[i]]$test, mtry)
    unclass(cf$logLik_NN)
}

tmp <- runif(10)

### learn and validate for each repetition of the simulated data
cf <- mclapply(1:length(simdat), cf_fun, mc.cores = ncores)
ret$cf <- unlist(cf)
print('Cforest is done!')
save(ret, file = ret_filename)

### for each repetition of the simulated data
### the parameters of the corresponding model,
### the true model likelihood and
### the estimated likelihoods of all forest models are saved
save(ret, file = ret_filename)

end_time <- Sys.time()
print(end_time - start_time)

sessionInfo()
