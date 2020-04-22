
rm(list=ls())

set.seed(12345)

start_time <- Sys.time()

library("parallel")
RNGkind("L'Ecuyer-CMRG")

ncores <- detectCores()

source("competitors.R")

if(!dir.exists("rda")){dir.create("rda")}

for(sp in seq_len(nrow(sim_para))){          # sp = 1

  load(paste("rda/simulated_data_nsim", sim_para[sp, "nsim"], ".rda", sep = ""))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define global parameters ------------------------------------------------
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ### joint parameters for all forest procedures
  
  ### number of repetitions
  nsim <- sim_para[sp, "nsim"]
  
  ## number of trees
  ntree <- sim_para[sp, "ntree"]
  
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
  
  res <- args[rep(1:nrow(args), each = sim_para[sp, "nsim"]),]
  simdat <- do.call("c", simdat)
  
  # true log-likelihood
  res$loglik <- unlist(do.call("c", loglik))
  
  # Ordinal Forests
  res$of_eq_ll <- res$of_eq_KL_I <- res$of_eq_KL_II <- res$of_eq_KL_III <- 0
  res$of_prop_ll <- res$of_prop_KL_I <- res$of_prop_KL_II <- res$of_prop_KL_III <- 0
  # Ordinal Transformation Forests
  res$tf_theta_ll <- res$tf_theta_KL_I <- res$tf_theta_KL_II <- res$tf_theta_KL_III <- 0
  res$tf_alpha_ll <- res$tf_alpha_KL_I <- res$tf_alpha_KL_II <- res$tf_alpha_KL_III <- 0
  res$cf_ll <- res$cf_KL_I <- res$cf_KL_II <- res$cf_KL_III <- 0
  
  # filename to save estimated likelihoods
  filename <- paste("rda/results_empeval_", rownames(sim_para)[sp], ".rda", sep = "")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Ordinal Forest (equal) --------------------------------------------------
  # ### perffunction = "equal"
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  of_equal_fun <- function(d) {
    mtry <- ifelse(res[d, "p"] <= 5, 
                   ncol(simdat[[d]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[d]]$train) - 1)))
    rf <- ordforest_equal(train = simdat[[d]]$train,
                          test = simdat[[d]]$test,
                          test_prb = simdat[[d]]$test_prb,
                          mtry = mtry)
    return(list(of_eq_ll = rf$of_eq_ll,
                of_eq_KL_I = rf$of_eq_KL_I,
                of_eq_KL_II = rf$of_eq_KL_II,
                of_eq_KL_III = rf$of_eq_KL_III))
  }
  
  tmp <- runif(10)
  
  ### train and test for each repetition of the simulated data
  of_eq_all <- mclapply(1:length(simdat), of_equal_fun, mc.cores = ncores)
  ### summarize and results
  of_eq_colnam <- c("of_eq_ll", "of_eq_KL_I", "of_eq_KL_II", "of_eq_KL_III")
  for(i in seq_along(of_eq_colnam)){
    res[, of_eq_colnam[i]] <- as.vector( do.call("rbind", lapply(of_eq_all, function(l) l[[of_eq_colnam[i]]])) )
  }
  print("ordinalForest (equal) is done!")
  save(res, file = filename)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Ordinal Forest (proportional) -------------------------------------------
  # ### perffunction = "proportional"
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  of_prop_fun <- function(d) {
    mtry <- ifelse(res[d, "p"] <= 5, 
                   ncol(simdat[[d]]$train) - 1, 
                   ceiling(sqrt(ncol(simdat[[d]]$train) - 1)))
    rf <- ordforest_prop(train = simdat[[d]]$train,
                         test = simdat[[d]]$test,
                         test_prb = simdat[[d]]$test_prb,
                         mtry = mtry)
    return(list(of_prop_ll = rf$of_prop_ll,
                of_prop_y = rf$of_prop_y,
                of_prop_KL_I = rf$of_prop_KL_I,
                of_prop_KL_II = rf$of_prop_KL_II,
                of_prop_KL_III = rf$of_prop_KL_III))
  }
  
  tmp <- runif(10)
  
  ### train and test for each repetition of the simulated data
  of_prop_all <- mclapply(1:length(simdat), of_prop_fun, mc.cores = ncores)
  ### summarize and results
  of_prop_colnam <- c("of_prop_ll", "of_prop_KL_I", "of_prop_KL_II", "of_prop_KL_III")
  for(i in seq_along(of_prop_colnam)){
    res[, of_prop_colnam[i]] <- as.vector( do.call("rbind", lapply(of_prop_all, function(l) l[[of_prop_colnam[i]]])) )
  }
  print("ordinalForest (proportional) is done!")
  save(res, file = filename)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Ordinal transformation forest with Bernstein basis and general score --------
  # Bs(theta)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  trtf_theta_fun <- function(d) {
    ## number of randomly preselected variables to perform a split
    ## for low-dimensional data it is set to the number of variables
    ## for high-dimensional data it is set to the square root of the number of variables
    mtry <- ifelse(res[d, "p"] <= 5,
                   ncol(simdat[[d]]$train) - 1,
                   ceiling(sqrt(ncol(simdat[[d]]$train) - 1)))
    rf <- traforest_theta(train = simdat[[d]]$train,
                          test = simdat[[d]]$test,
                          test_prb = simdat[[d]]$test_prb,
                          mtry = mtry)
    return(list(tf_theta_ll = rf$tf_theta_ll,
                tf_theta_KL_I = rf$tf_theta_KL_I,
                tf_theta_KL_II = rf$tf_theta_KL_II,
                tf_theta_KL_III = rf$tf_theta_KL_III))
  }
  
  tmp <- runif(10)
  
  ### train and test for each repetition of the simulated data
  tf_theta_all <- mclapply(1:length(simdat), trtf_theta_fun, mc.cores = ncores)
  ### summarize and results
  tf_theta_colnam <- c("tf_theta_ll", "tf_theta_KL_I", "tf_theta_KL_II", "tf_theta_KL_III")
  for(i in seq_along(tf_theta_colnam)){
    res[, tf_theta_colnam[i]] <- as.vector( do.call("rbind", lapply(tf_theta_all, function(l) l[[tf_theta_colnam[i]]])) )
  }
  print("trtf_theta is done!")
  save(res, file = filename)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Ordinal transformation forest with Bernstein basis and log-rank score ----------
  # Bs(alpha)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  trtf_alpha_fun <- function(d) {
    mtry <- ifelse(res[d, "p"] <= 5,
                   ncol(simdat[[d]]$train) - 1,
                   ceiling(sqrt(ncol(simdat[[d]]$train) - 1)))
    rf <- traforest_alpha(train = simdat[[d]]$train,
                          test = simdat[[d]]$test,
                          test_prb = simdat[[d]]$test_prb,
                          mtry = mtry)
    return(list(tf_alpha_ll = rf$tf_alpha_ll,
                tf_alpha_KL_I = rf$tf_alpha_KL_I,
                tf_alpha_KL_II = rf$tf_alpha_KL_II,
                tf_alpha_KL_III = rf$tf_alpha_KL_III))
  }
  
  tmp <- runif(10)
  
  ### train and test for each repetition of the simulated data
  tf_alpha_all <- mclapply(1:length(simdat), trtf_alpha_fun, mc.cores = ncores)
  ### summarize and results
  tf_alpha_colnam <- c("tf_alpha_ll", "tf_alpha_KL_I", "tf_alpha_KL_II", "tf_alpha_KL_III")
  for(i in seq_along(tf_alpha_colnam)){
    res[, tf_alpha_colnam[i]] <- as.vector( do.call("rbind", lapply(tf_alpha_all, function(l) l[[tf_alpha_colnam[i]]])) )
  }
  print("trtf_alpha is done!")
  save(res, file = filename)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Conditional Inference Forests -------------------------------------------
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cf_fun <- function(d){
    mtry <- ifelse(res[d, "p"] <= 5,
                   ncol(simdat[[d]]$train) - 1,
                   ceiling(sqrt(ncol(simdat[[d]]$train) - 1)))
    rf <- mycforest(train = simdat[[d]]$train,
                    test = simdat[[d]]$test,
                    test_prb = simdat[[d]]$test_prb,
                    mtry = mtry)
    return(list(cf_ll = rf$cf_ll,
                cf_KL_I = rf$cf_KL_I,
                cf_KL_II = rf$cf_KL_II,
                cf_KL_III = rf$cf_KL_III))
  }
  
  tmp <- runif(10)
  
  ### train and test for each repetition of the simulated data
  cf_all <- mclapply(1:length(simdat), cf_fun, mc.cores = ncores)
  ### summarize and results
  cf_colnam <- c("cf_ll", "cf_KL_I", "cf_KL_II", "cf_KL_III")
  for(i in seq_along(cf_colnam)){
    res[, cf_colnam[i]] <- as.vector( do.call("rbind", lapply(cf_all, function(l) l[[cf_colnam[i]]])) )
  }
  print("mycforest is done!")
  save(res, file = filename)
  
  ### delete no longer used variables
  rm(list = c("args", "simdat", "loglik",
              "of_eq_all", "of_prop_all",
              "tf_theta_all", "tf_alpha_all", "cf_all"))

}

end_time <- Sys.time()
print(end_time - start_time)

sessionInfo()
