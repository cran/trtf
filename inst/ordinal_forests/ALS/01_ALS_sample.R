
## load ALS data
## The data is available to registered users of https://nctu.partners.org/ProACT
## Register, download the data, and run
## R> library("TH.data"); demo("PROACT");
load("data/ALSFRSdata.rda")

set.seed(12345)

tmp <- ALSFRSdata

# <Field Name="10. Respiratory" ID="1214">
tmp$y <- tmp[["v_1214.halfyearafter"]]
tmp <- tmp[, -grep("halfyearafter", colnames(tmp), ignore.case = TRUE)]
tmp <- tmp[complete.cases(tmp),]

# dim(tmp)
# [1] 1013   71

NSim <- 120

learn <- list()
for (i in 1:NSim){
  while(TRUE) {
      learn[[i]] <- sample(1:nrow(tmp), ceiling(nrow(tmp) * .75))
      if (nlevels(tmp$y[learn[[i]], drop = TRUE]) == nlevels(tmp$y))
      break
  }
}

save(tmp, learn, file="ALS_samples.rda")

