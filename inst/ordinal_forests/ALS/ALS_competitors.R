
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load packages --------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library("trtf")
library("party")
library("partykit")
library("tram")
library("ranger")
library("ordinalForest")

library("psych")
# https://personality-project.org/r/psych/
# https://personality-project.org/r/psych-manual.pdf

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cohen's Kappa --------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ckappa <- function(yobs1, yobs2) {
  
  x <- data.frame(yobs1 = yobs1, yobs2 = yobs2)
  
  # unweighted kappa
  ck <- tryCatch(cohen.kappa(x)$kappa,
                 error = function(e) NA)
  
  # weighted kappa
  J <- length(unique(c(x$yobs1, x$yobs2)))
  
  # linear weighted kappa
  w_lin <- matrix(0, ncol = J, nrow = J)
  w_lin[] <- abs((col(w_lin) - row(w_lin)))
  w_lin <- 1 - w_lin/(J - 1)
  lk <- tryCatch(cohen.kappa(x, w = w_lin)$weighted.kappa,
                 error = function(e) NA)
  
  # quadratic weighted kappa
  w_quad <- matrix(0, ncol = J, nrow = J)
  w_quad[] <- abs((col(w_quad) - row(w_quad)))^2
  w_quad <- 1 - w_quad/(J - 1)^2
  qk <- tryCatch(cohen.kappa(x, w = w_quad)$weighted.kappa,
                 error = function(e) NA)
  
  return(list(ck = ck,
              lk = lk,
              qk = qk))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal Forest - perffunction = "equal" ------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ordforest_equal <- function(train, test, mtry) {
  
  rf <- ordfor(depvar = "y", data = train,
               nsets = 1000,       # number of considered score sets
               ntreeperdiv = 100,  # number of trees considered per tried score
               ntreefinal = ntree, # number of trees for the forest
               perffunction = "equal",
               min.node.size = max(minbucket, nrow(train) / 2^nodedepth), 
               mtry = mtry, 
               keep.inbag = TRUE)
  
  # model-based probability density function and y predictions
  rf_prb <- predict(rf, newdata = test)$classfreqtree
  
  # log-likelihood
  prb <- rf_prb[cbind(1:nrow(test), unclass(test$y))]
  of_eq_ll <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  return(list(of_eq_ll = of_eq_ll))
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal Forest - perffunction = "proportional" -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ordforest_prop <- function(train, test, mtry) {
  
  rf <- ordfor(depvar = "y", data = train,
               nsets = 1000,       # number of considered score sets
               ntreeperdiv = 100,  # number of trees considered per tried score
               ntreefinal = ntree, # number of trees for the forest
               perffunction = "proportional",
               min.node.size = max(minbucket, nrow(train) / 2^nodedepth), 
               mtry = mtry, 
               keep.inbag = TRUE)
  
  # model-based probability density function and y predictions
  rf_prb <- predict(rf, newdata = test)$classfreqtree
  rf_y <- factor(predict(rf, newdata = test)$ypred,
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  
  # log-likelihood
  prb <- rf_prb[cbind(1:nrow(test), unclass(test$y))]
  of_prop_ll <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  return(list(of_prop_ll = of_prop_ll,
              of_prop_y = rf_y))  # later needed for Cohen's Kappa
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal transfromation forest with Bernstein basis and general score -------
# Bs(theta) --> non-proportional odds deviations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

traforest_theta <- function(train, test, mtry) {
  
  m0 <- as.mlt(Polr(y ~ 1, data = train))
  rf <- traforest(m0, 
                  formula = y ~ ., 
                  data = train,
                  ntree = ntree, 
                  trace = FALSE, 
                  mtry = mtry,
                  control = ctrl_partykit)
  
  # model-based probability density function and y predictions
  rf_prb <- do.call("rbind", lapply(predict(rf, newdata = test, type = "density"), c))
  rf_prb <- data.frame(rf_prb)
  colnames(rf_prb) <- levels(test$y)
  rf_y <- factor(predict.cforest(rf, newdata = test, type = "response"),
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  names(rf_y) <- NULL
  # yprob <- do.call("rbind", lapply(predict(rf, newdata = test, type = "density"), t))
  # yprob <- data.frame(yprob)
  # colnames(yprob) <- levels(test$y)
  # yhat <- factor(levels(test$y)[apply(yprob, 1, which.max)],
  #                levels = levels(test$y), order = TRUE)
  # table(yhat)
  
  # log-likelihood with “nearest neighbour” forest weights prediction
  cf <- predict(rf, newdata = test, type = "coef")
  tf_theta_ll <- logLik(rf, newdata = test, coef = cf)
  
  return(list(tf_theta_ll = tf_theta_ll,
              tf_theta_y = rf_y)) # later needed for Cohen's Kappa
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal transformation forest with Bernstein basis and log-rank score ------
# Bs(alpha) --> proportional odds deviations from the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

traforest_alpha <- function(train, test, mtry) {
  
  ### There is no intercept in Bernstein-basis, therefor to perform
  ### log-rank splitting constant variable is added to the data implicitly
  train$alpha <- 1
  test$alpha <- 1
  m0 <- as.mlt(Polr(y ~ alpha, data = train, fixed = c("alpha" = 0)))
  ### split wrt to intercept ("log-rank scores") only
  rf <- traforest(m0, 
                  formula = y | alpha ~ ., 
                  data = train,
                  parm = "alpha",
                  mltargs = list(fixed = c("alpha" = 0)),
                  ntree = ntree, 
                  trace = FALSE, 
                  mtry = mtry,
                  control = ctrl_partykit)
  
  # model-based probability density function and y predictions
  rf_prb <- do.call("rbind", lapply(predict(rf, newdata = test,
                                            mnewdata = data.frame(alpha = 1),
                                            type = "density"), c))
  rf_prb <- data.frame(rf_prb)
  colnames(rf_prb) <- levels(test$y)
  rf_y <- factor(predict.cforest(rf, newdata = test, type = "response")[,"y"],
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  
  # log-likelihood with “nearest neighbour” forest weights prediction
  cf <- predict(rf, newdata = test, type = "coef")
  tf_alpha_ll <- logLik(rf, newdata = test, coef = cf)
  
  return(list(tf_alpha_ll = tf_alpha_ll))
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conditional Inference Trees --> Forests ------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mycforest <- function(train, test, mtry) {
  
  ctrl_party <- party:::cforest_unbiased(ntree = ntree, 
                                         maxdepth = nodedepth, 
                                         mtry = mtry,
                                         minbucket = minbucket)
  
  ### log-rank splitting
  rf <- party::cforest(y ~ ., data = train, 
                       control = ctrl_party)
  
  # model-based probability density function and y predictions
  rf_prb <- predict(rf, newdata = test, type = "prob")
  rf_prb <- do.call("rbind", rf_prb)
  
  # log-likelihood
  prb <- rf_prb[cbind(1:nrow(test), unclass(test$y))]
  cf_ll <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  return(list(cf_ll = cf_ll))
  
}

