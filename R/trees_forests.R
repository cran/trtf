
.ctmfit <- function(object, parm, mltargs, reparm, 
                    min_update = length(coef(object)) * 2) {
    
    ctmobject <- object

    ### note: control is never used but expected by partykit::ctree
    function(data, weights, control, ...) {
        mf <- model.frame(data, yxonly = TRUE)
        iy <- data[["yx", type = "index"]]

        mltargs$data <- mf
        ctmobject <- do.call("mlt", mltargs)
        thetastart <- coef(ctmobject, fixed = FALSE)

        function(subset = NULL, weights = NULL, info = NULL, model = FALSE, 
                 estfun = TRUE, object = FALSE) {
            if (model) return(list(object = ctmobject, iy = iy))

            if (!is.null(iy)) {
                w <- libcoin::ctabs(iy, weights = weights, subset = subset)[-1L]
                subset <- NULL
            } else {
                if (is.null(weights) || length(weights) == 0) {
                    w <- rep(1, nrow(mf))
                } else {
                    w <- weights
                }
                w[-subset] <- 0 ### mlt >= 1.0-3 allows subset but we
                ### still need the weights to be zero for some operations below
            }
            if (!is.null(info$coef)) {
                thetastart <- info$coef
            } else {
                thetastart <- coef(ctmobject, fixed = FALSE)
            }
            ### update parameters when there are more than min_update
            ### observations
            nw <- if (is.null(subset)) sum(w) else sum(w[subset])
            if (nw > min_update) {
                umod <- suppressWarnings(try(update(ctmobject, weights = w, 
                    subset = subset, theta = thetastart), silent = TRUE))
                if (inherits(umod, "try-error") || umod$convergence != 0) {
                    umod <- suppressWarnings(try(update(ctmobject, weights = w, subset = subset), 
                                                  silent = TRUE))
                    if (inherits(umod, "try-error") || umod$convergence != 0) {
                        mltargs$weights <- w
                        ### [1]: no subset allowed here, so used zero weights (see
                        ### above AND below)!!!
                        umod <- try(do.call("mlt", mltargs))
                    }
                }
            } else {
                ### reuse parameters from root node (next if statement)
                umod <- 0
                class(umod) <- "try-error"
            }
            if (inherits(umod, "try-error")) {
                ### we badly need some estimate in each node, even if fitting
                ### fails
                if (!estfun) 
                    return(list(coef = thetastart, objfun = NA,
                                converged = FALSE))
                ### estfun was requested and thus a tree wants to look for
                ### splits. We cannot update the model, so we simply reuse
                ### the start parameters. This amounts to treating the
                ### scores as fix when building the subtree. This was always
                ### the case in party::ctree (the "model" was once estimated
                ### in the root node) and Stefan Wager suggested to apply
                ### this scheme in model-based trees in a discussion 2019-09-06.
                # mltargs$weights <- w
                # mltargs$doFit <- FALSE
                # umod <- try(do.call("mlt", mltargs))
                umod <- ctmobject
                umod$weights <- w
                umod$subset <- subset
                coef(umod)[names(thetastart)] <- thetastart
                umod$convergence <- 0L ### this is FAKE but extree would stop
                                       ### after non-convergence
            }
            ret <- NULL
            if (estfun) {
                ret <- estfun(umod, parm = coef(umod, fixed = TRUE))[, parm, drop = FALSE]
                if (!is.null(subset)) {
                    if (NROW(ret) == length(subset)) {
                        tmp <- matrix(0, nrow = length(w), 
                                      ncol = ncol(ret))
                        tmp[subset,] <- ret
                        ret <- tmp
                    } else {
                        ### see [1]
                        ret[-subset,] <- 0
                    }                  
                }
                if (!is.null(iy)) ret <- rbind(0, ret)
                if (!is.null(reparm)) ret <- ret %*% reparm
            }
            return(list(estfun = ret, 
                        coefficients = coef(umod, fixed = FALSE), 
                        ### we always minimise risk
                        objfun = -logLik(umod), 
                        object = if (object) umod else NULL,
                        converged = isTRUE(all.equal(umod$convergence, 0))))
        }
    }
} 

trafotree <- function(object, parm = 1:length(coef(object)), reparm = NULL, 
                      min_update = length(coef(object)) * 2, 
                      mltargs = list(maxit = 10000), ...) {

    ### we only work with the ctm object
    if (inherits(object, "mlt")) {
        if (is.null(mltargs$scale))
            mltargs$scale <- object$scale
        object <- object$model
    }
    ### this is tricky because parm is only valid
    ### for this ctm object (not not for tram-like objects)
    mltargs$model <- object
    ### note: weights, offset, cluster etc. are evaluated here !!!
    args <- list(...)
    args$ytrafo <- .ctmfit(object = object, parm = parm, 
                           mltargs = mltargs, reparm = reparm, 
                           min_update = min_update)
    args$update <- TRUE
    ret <- do.call("ctree", args)
    ret$model <- object
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE)
    ret$mltargs <- mltargs

    weights <- data_party(ret)[["(weights)"]]
    if (is.null(weights)) weights <- rep(1, nrow(data_party(ret)))

    ### store coefs and logLik _outside_ tree
    ### <FIXME> this will cause problems with nodeprune </FIXME>
    nd <- predict(ret, type = "node")
    ret$models <- tapply(1:length(nd), factor(nd), function(i) 
        ret$trafo(i, weights = weights, estfun = FALSE)) ### note: trafo needs weights
    ret$coef <- do.call("rbind", lapply(ret$models, function(x) x$coef))
    ret$logLik <- sapply(ret$models, function(x) -x$objfun) ### objfun is neg logLik

    class(ret) <- c("trafotree", class(ret))
    ret
}

traforest <- function(object, parm = 1:length(coef(object)), reparm = NULL,
                      update = TRUE, min_update = length(coef(object)) * 2, 
                      mltargs = list(maxit = 10000), ...) {

    if (inherits(object, "mlt")) object <- object$model
    ### this is tricky because parm is only valid
    ### for this ctm object (not not for tram-like objects)
    mltargs$model <- object
    ### note: weights, offset, cluster etc. are evaluated here !!!
    args <- list(...)
    args$ytrafo <- .ctmfit(object = object, parm = parm, 
                           mltargs = mltargs, reparm = reparm,
                           min_update = min_update)
    args$update <- update
    ret <- do.call("cforest", args)
    ret$model <- object
    ret$mltargs <- mltargs
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE, object = TRUE)
    class(ret) <- c("traforest", class(ret))
    ret
}
