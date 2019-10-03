

# Load packages -----------------------------------------------------------

library("trtf")
library("party")
library("partykit")
library("tram")

library("ranger")
library("ordinalForest")

set.seed(290875)

# Ordinal Forest ----------------------------------------------------------

myordforest <- function(learn, test, mtry) {

    rf <- ordfor(depvar = "y", data = learn,
                 nsets = 1000,       # number of considered score sets
                 ntreeperdiv = 100,  # number of trees considered per tried score
                 ntreefinal = ntree, # number of trees for the forest
                 perffunction = "equal",
                 min.node.size = max(minbucket, nrow(learn) / 2^nodedepth), 
                 mtry = mtry, 
                 keep.inbag = TRUE)

    prb <- predict(rf, newdata = test)$class
    prb <- prb[cbind(1:nrow(test), unclass(test$y))]
    llNN <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))

    list(logLik_NN = llNN)
}

myordforest_prop <- function(learn, test, mtry) {
  
  rf <- ordfor(depvar = "y", data = learn,
               nsets = 1000,       # number of considered score sets
               ntreeperdiv = 100,  # number of trees considered per tried score
               ntreefinal = ntree, # number of trees for the forest
               perffunction = "proportional",
               min.node.size = max(minbucket, nrow(learn) / 2^nodedepth), 
               mtry = mtry, 
               keep.inbag = TRUE)
  
  prb <- predict(rf, newdata = test)$class
  prb <- prb[cbind(1:nrow(test), unclass(test$y))]
  llNN <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  list(logLik_NN = llNN)
}


# Ordinal transfromation forest with Bernstein basis and general score --------
# Bs(theta) --> on-proportional odds deviations

mytraforest <- function(learn, test, mtry) {
    
    ### the largest values of time are not used for Bernstein approximation
    m0 <- as.mlt(Polr(y ~ 1, data = learn))

    ### traforest generates the subsampling weights later to be used
    ### by other methods
    rf <- traforest(m0, 
                    formula = y ~ ., 
                    data = learn,
                    ntree = ntree, 
                    trace = FALSE, 
                    mtry = mtry,
                    control = ctrl_partykit)

    ### nearest neighbour predictions
    cf <- predict(rf, newdata = test, type = "coef")
    llNN <- logLik(rf, newdata = test, coef = cf)
	
    list(logLik_NN = llNN, mtry = mtry)

}


# Ordinal transformation forest with Bernstein basis and log-rank score ----------
# Bs(alpha) --> proportional odds deviations from the model

mytraforest_B <- function(learn, test, mtry) {
    
    ### There is no intercept in Bernstein-basis, therefore to perform
    ### log-rank splitting constant variable is added to the data implicitly
    learn$alpha <- 1
    test$alpha <- 1
    
    m0 <- as.mlt(Polr(y ~ alpha, data = learn, fixed = c("alpha" = 0)))

    ### split wrt to intercept ("log-rank scores") only
    rf <- traforest(m0, 
                    formula = y | alpha ~ ., 
                    data = learn,
                    parm = "alpha",
                    mltargs = list(fixed = c("alpha" = 0)),
                    ntree = ntree, 
                    trace = FALSE, 
                    mtry = mtry,
                    control = ctrl_partykit)

    ### nearest neighbor predictions
    cf <- predict(rf, newdata = test, type = "coef")
    llNN <- logLik(rf, newdata = test, coef = cf)

    list(logLik_NN = llNN)
}


# Conditional Inference Trees --> Forests ---------------------------------

mycforest <- function(learn, test, mtry) {

    ctrl_party <- party:::cforest_unbiased(ntree = ntree, 
                                           maxdepth = nodedepth, 
                                           mtry = mtry,
                                           minbucket = minbucket)

    ### log-rank splitting
    rf <- party::cforest(y ~ ., data = learn, 
                         control = ctrl_party)

    prb <- predict(rf, newdata = test, type = "prob")
    prb <- do.call("rbind", prb)
    prb <- prb[cbind(1:nrow(test), unclass(test$y))]
    llNN <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
    
    list(logLik_NN = llNN)

}


# Example -----------------------------------------------------------------

if (FALSE) {
  
data("mammoexp", package = "TH.data")
colnames(mammoexp)[colnames(mammoexp) == "ME"] <- "y"


ntree <- 25
test <- sample(1:NROW(mammoexp), 25)

tf <- mytraforest(mammoexp[-test,], mammoexp[test,], mtry = 3)
tfB <- mytraforest_B(tf, mammoexp[-test,], mammoexp[test,])
cf <- mycforest(tf, mammoexp[-test,], mammoexp[test,])
of <- myordforest(tf, mammoexp[-test,], mammoexp[test,])

tf[1]
tfB[1]
cf[1]
of[1]

}

