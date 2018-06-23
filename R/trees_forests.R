
.ctmfit <- function(object, parm, mltargs) {
    
    ctmobject <- object

    function(data, weights, control, ...) {
        mf <- model.frame(data, yxonly = TRUE)
        iy <- data[["yx", type = "index"]]

        mltargs$data <- mf
        ctmobject <- do.call("mlt", mltargs)
        thetastart <- coef(ctmobject)

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
                thetastart <- coef(ctmobject)
            }
            umod <- suppressWarnings(try(update(ctmobject, weights = w, subset = subset, theta = thetastart), silent = TRUE))
            if (inherits(umod, "try-error") || umod$convergence != 0) {
                umod <- suppressWarnings(try(update(ctmobject, weights = w, subset = subset), silent = TRUE))
                if (inherits(umod, "try-error") || umod$convergence != 0) {
                    mltargs$weights <- w
                    ### no subset allowed here, so used zero weights (see
                    ### above)!!!
                    umod <- try(do.call("mlt", mltargs))
                }
            }
            if (inherits(umod, "try-error")) {
                ### we badly need some estimate in each node, even if fitting
                ### fails
                if (!estfun) 
                    return(list(coef = thetastart, objfun = NA,
                                converged = FALSE))
                return(list(estfun = matrix(0, nrow = nrow(mf), ncol = length(w)),
                            coef = thetastart, objfun = NA,  converged = FALSE))
            }
            ret <- NULL
            if (estfun) {
                ret <- estfun(umod)[, parm, drop = FALSE]
                if (!is.null(subset)) {
                    tmp <- matrix(0, nrow = length(w), 
                                  ncol = ncol(ret))
                    tmp[subset,] <- ret
                    ret <- tmp
                }
                if (!is.null(iy)) ret <- rbind(0, ret)
            }
            return(list(estfun = ret, 
                        coefficients = coef(umod), 
                        ### we always minimise risk
                        objfun = -logLik(umod), 
                        object = if (object) umod else NULL,
                        converged = isTRUE(all.equal(umod$convergence, 0))))
        }
    }
} 

trafotree <- function(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), ...) {

    if (inherits(object, "mlt")) object <- object$model
    mltargs$model <- object
    ### note: weights, offset, cluster etc. are evaluated here !!!
    args <- list(...)
    args$ytrafo <- .ctmfit(object, parm, mltargs)
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

traforest <- function(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), update = TRUE, ...) {

    if (inherits(object, "mlt")) object <- object$model
    mltargs$model <- object
    ### note: weights, offset, cluster etc. are evaluated here !!!
    args <- list(...)
    args$ytrafo <- .ctmfit(object, parm, mltargs)
    args$update <- update
    ret <- do.call("cforest", args)
    ret$model <- object
    ret$mltargs <- mltargs
    ret$mltobj <- ret$trafo(model = TRUE, estfun = FALSE, object = TRUE)
    class(ret) <- c("traforest", class(ret))
    ret
}
