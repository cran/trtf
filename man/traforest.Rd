\name{traforest}
\alias{traforest}
\title{Transformation Forests}
\description{
  Partitioned and aggregated transformation models
}
\usage{
traforest(object, parm = 1:length(coef(object)), mltargs = list(maxit = 10000), 
          update = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{\link[mlt]{ctm}} or \code{\link[mlt]{mlt}} specifying the
                abstract model to be partitioned.}
  \item{parm}{parameters of \code{object} those corresponding score is
              used for finding partitions.}
  \item{mltargs}{arguments to \code{\link[mlt]{mlt}} for fitting the
                 transformation models.}
  \item{update}{logical, if \code{TRUE}, models and thus scores are updated in
                every node. If \code{FALSE}, the model and scores are
                computed once in the root node. The latter option is faster
                but less accurate.}
  \item{\dots}{arguments to \code{\link[partykit]{cforest}}, at least
               \code{formula} and \code{data}.}
}
\details{
  Conditional inference trees are used for partitioning likelihood-based transformation
  models as described in Hothorn and Zeileis (2017). The method can be seen
  in action in Hothorn (2018) and the corresponding code is available as
  \code{demo("BMI")}. 

}
\value{
  An object of class \code{traforest} with corresponding \code{logLik} and
  \code{predict} methods.
}
\references{
  Torsten Hothorn and Achim Zeileis (2017). Transformation Forests.
  \url{https://arxiv.org/abs/1701.02110}.

  Torsten Hothorn (2018). Top-Down Transformation Choice. \emph{Statistical
  Modelling}, \url{https://arxiv.org/abs/1706.08269}.
}
\examples{

### Example: Personalised Medicine Using Partitioned and Aggregated Cox-Models
### A combination of <DOI:10.1177/0962280217693034> and <arXiv:1701.02110>
### based on infrastructure in the mlt R add-on package described in
### https://cran.r-project.org/web/packages/mlt.docreg/vignettes/mlt.pdf

library("trtf")
library("survival")
### German Breast Cancer Study Group 2 data set
data("GBSG2", package = "TH.data")

### set-up Cox model with overall treatment effect in hormonal therapy
yvar <- numeric_var("y", support = c(100, 2000), bounds = c(0, Inf))
By <- Bernstein_basis(yvar, order = 5, ui = "incre")
m <- ctm(response = By, shifting = ~ horTh, todistr = "MinExt", data = GBSG2)
GBSG2$y <- with(GBSG2, Surv(time, cens))

### overall log-hazard ratio
coef(cmod <- mlt(m, data = GBSG2))["horThyes"]
### roughly the same as 
coef(coxph(y ~ horTh, data = GBSG2))

\dontrun{

### estimate age-dependent Cox models (here ignoring all other covariates)
ctrl <- ctree_control(minsplit = 50, minbucket = 20, mincriterion = 0)
set.seed(290875)
tf_cmod <- traforest(m, formula = y ~ horTh | age, control = ctrl, 
                     ntree = 50, mtry = 1, trace = TRUE, data = GBSG2)

### plot age-dependent treatment effects vs. overall treatment effect
nd <- data.frame(age = 30:70)
cf <- predict(tf_cmod, newdata = nd, type = "coef")
nd$logHR <- sapply(cf, function(x) x["horThyes"])
plot(logHR ~ age, data = nd, pch = 19, xlab = "Age", ylab = "log-Hazard Ratio")
abline(h = coef(cmod <- mlt(m, data = GBSG2))["horThyes"])
### treatment most beneficial in very young patients
### NOTE: scale of log-hazard ratios depends on
### corresponding baseline hazard function which  _differs_
### across age; interpretation of positive / negative treatment effect is,
### however, save.
}

}
\keyword{trees}
