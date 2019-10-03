rm(list=ls())

start_time <- Sys.time()

set.seed(12345)

load("ALS_samples.rda")

idx <- i <- INDEX

# number of trees
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

res_filename <- paste('ALS_results_combined_', i, '.rda', sep = "")

NSim <- length(learn)

source("../empeval/competitors.R")

## number of variables to choose from in a node
mtry <- ceiling(sqrt(ncol(tmp) - 2))

ret <- data.frame(polrconst =  rep(0, NSim))
ret$polr <- ret$of <- ret$of_prop <- ret$cf <- ret$tf_B_alpha <- ret$tf_B <- 0
# head(ret)

### Models without treatment
### erase Riluzole from the data
riluzole <- tmp$Riluzole
tmp$Riluzole <- NULL


### Unconditional Polr model: intercept only
polrc_fun <- function(i) {
  w <- 1:nrow(tmp) %in% learn[[i]] + 0L
  cph <- Polr(y ~ 1, data = tmp, weights = w)
  unclass(logLik(cph, w = 1 - w))
}
polrint <- lapply(idx, polrc_fun)
ret$polrconst[idx] <- unlist(polrint)


### Prognostic Polr with exp(b)
polr_fun <- function(i) {
  w <- 1:nrow(tmp) %in% learn[[i]] + 0L
  cph <- Polr(y ~ ., data = tmp, weights = w)
  unclass(logLik(cph, w = 1 - w))
}
polrall <- lapply(idx, polr_fun)
ret$polr[idx] <- unlist(polrall)


### Transformation forest OTF(alpha)
trtf_fun <- function(i) {
  tf <- tryCatch(mytraforest(tmp[learn[[i]],], 
                             tmp[-learn[[i]],], 
                             mtry = mtry),
                 error = function(e) {NA})
  tryCatch(unclass(tf$logLik_NN), error = function(e) {NA})
}
tf <- lapply(idx, trtf_fun)
ret$tf_B[idx] <- unlist(tf)
print('Traforest Bs is done!')
save(ret, file = res_filename)


### Transformation forest OTF(theta)
trtf_fun_B <- function(i) {
  tf <- tryCatch(mytraforest_B(tmp[learn[[i]],],
                               tmp[-learn[[i]],],
                               mtry = mtry),
                 error = function(e) {NA})
  tryCatch(unclass(tf$logLik_NN), error = function(e) {NA})
}

tf <- lapply(idx, trtf_fun_B)
ret$tf_B_alpha[idx] <- unlist(tf)
print('Traforest B_alpha is done!')
save(ret, file = res_filename)


### Ordinal Forest (equally)
library("ranger")
library("ordinalForest")
rg_fun <- function(i) {
  rg <- tryCatch(myordforest(tmp[learn[[i]],], tmp[-learn[[i]],], mtry = mtry),
                 error = function(e) {NA})
  tryCatch(unclass(rg$logLik_NN), error = function(e) {NA})
}

rg <- lapply(idx, rg_fun)
ret$of[idx] <- unlist(rg)
print('Ranger is done!')
save(ret, file = res_filename)


### Ordinal Forest (proportional)
library("ranger")
library("ordinalForest")
rg_fun_prop <- function(i) {
  rg <- tryCatch(myordforest_prop(tmp[learn[[i]],], tmp[-learn[[i]],], mtry = mtry),
                 error = function(e) {NA})
  tryCatch(unclass(rg$logLik_NN), error = function(e) {NA})
}

rg <- lapply(idx, rg_fun_prop)
ret$of_prop[idx] <- unlist(rg)
print('Ranger is done!')
save(ret, file = res_filename)


### ConditionalForest
cf_fun <- function(i) {
  cf <- tryCatch(mycforest(tmp[learn[[i]],], tmp[-learn[[i]],], mtry = mtry),
                 error = function(e) {NA})
  tryCatch(unclass(cf$logLik_NN), error = function(e) {NA})
}

cf <- lapply(idx, cf_fun)
ret$cf[idx] <- unlist(cf)
print('CForest is done!')
save(ret, file = res_filename)

end_time <- Sys.time()

print(end_time - start_time)

sessionInfo()
