rm(list=ls())

start_time <- Sys.time()

set.seed(12345)

load("ALS_samples.rda")

idx <- i <- INDEX

source("ALS_competitors.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define global parameters ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

if(!dir.exists("rda")){dir.create("rda")}

filename <- paste("rda/ALS_results_combined_", i, ".rda", sep = "")

NSim <- length(learn)

## number of variables to choose from in a node
mtry <- ceiling(sqrt(ncol(tmp) - 2))

res <- data.frame(polrconst =  rep(0, NSim))

res$polrconst <- res$polr <- 0
res$of_eq_ll <- res$of_prop_ll <- 0
res$tf_theta_ll <- res$tf_alpha_ll <- res$cf_ll <- 0

res$tf_of_qk <- res$tf_of_lk <- res$tf_of_ck <- 0

# head(res)

### Models without treatment
### erase Riluzole from the data
riluzole <- tmp$Riluzole
tmp$Riluzole <- NULL

### Unconditional Polr model: intercept only
polrc_fun <- function(i) {
  w <- 1:nrow(tmp) %in% learn[[i]] + 0L
  uncond_polr <- Polr(y ~ 1, data = tmp, weights = w)
  unclass(logLik(uncond_polr, w = 1 - w))
}
polrconst <- lapply(idx, polrc_fun)
res$polrconst[idx] <- unlist(polrconst)

### Prognostic Polr with exp(b)
polr_fun <- function(i) {
  w <- 1:nrow(tmp) %in% learn[[i]] + 0L
  prog_polr <- Polr(y ~ ., data = tmp, weights = w)
  unclass(logLik(prog_polr, w = 1 - w))
}
polrall <- lapply(idx, polr_fun)
res$polr[idx] <- unlist(polrall)

### Ordinal Forest (equally)
of_equal_fun <- function(i) {
  rf <- tryCatch(ordforest_equal(train = tmp[learn[[i]],],
                                 test = tmp[-learn[[i]],],
                                 mtry = mtry),
                 error = function(e) {NA})
  return(as.numeric(unclass(rf$of_eq_ll)))
}
of_eq_all <- lapply(idx, of_equal_fun)
res$of_eq_ll[idx] <- unlist(of_eq_all)
print('of_equal_fun is done!')
save(res, file = filename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Ordinal Forest (proportional)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

of_prop_fun <- function(i) {
  rf <- tryCatch(ordforest_prop(train = tmp[learn[[i]],],
                                test = tmp[-learn[[i]],],
                                mtry = mtry),
                 error = function(e) {NA})
  return(list(of_prop_ll = rf$of_prop_ll,
              of_prop_y = rf$of_prop_y))
}
of_prop_all <- lapply(idx, of_prop_fun)
res$of_prop_ll[idx] <- as.numeric(of_prop_all[[1]]$of_prop_ll)
of_prop_y <- of_prop_all[[1]]$of_prop_y
print('ordforest_prop is done!')
save(res, file = filename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Transformation forest OTF(theta)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trtf_theta_fun <- function(i) {
  rf <- tryCatch(traforest_theta(train = tmp[learn[[i]],],
                                 test = tmp[-learn[[i]],],
                                 mtry = mtry),
                 error = function(e) {NA})
  return(list(tf_theta_ll = rf$tf_theta_ll,
              tf_theta_y = rf$tf_theta_y))
}
tf_theta_all <- lapply(idx, trtf_theta_fun)
res$tf_theta_ll[idx] <- as.numeric(tf_theta_all[[1]]$tf_theta_ll)
tf_theta_y <- tf_theta_all[[1]]$tf_theta_y
print("trtf_theta is done!")
save(res, file = filename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Cohen's Kappa
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

kappa_tf_of <- ckappa(yobs1 = of_prop_y,
                      yobs2 = tf_theta_y)
res[idx, c("tf_of_ck",                                # unweighted kappa
           "tf_of_lk",                                # linear weighted kappa
           "tf_of_qk")] <- do.call("c", kappa_tf_of)  # quadratic weighted kappa

### save results
print("Cohen's Kappa(s) are done!")
save(res, file = filename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Transformation forest OTF(alpha)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

trtf_alpha_fun <- function(i) {
  rf <- tryCatch(traforest_alpha(train = tmp[learn[[i]],], 
                                 test = tmp[-learn[[i]],], 
                                 mtry = mtry),
                 error = function(e) {NA})
  return(as.numeric(unclass(rf$tf_alpha_ll)))
}
tf_alpha_all <- lapply(idx, trtf_alpha_fun)
res$tf_alpha_ll[idx] <- unlist(tf_alpha_all)
print("trtf_alpha is done!")
save(res, file = filename)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Conditional Inference Forests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cf_fun <- function(i) {
  rf <- tryCatch(mycforest(train = tmp[learn[[i]],],
                           test = tmp[-learn[[i]],],
                           mtry = mtry),
                 error = function(e) {NA})
  return(as.numeric(unclass(rf$cf_ll)))
}
cf_all <- lapply(idx, cf_fun)
res$cf_ll[idx] <- unlist(cf_all)
print('cforest is done!')
save(res, file = filename)

end_time <- Sys.time()

print(end_time - start_time)

sessionInfo()
