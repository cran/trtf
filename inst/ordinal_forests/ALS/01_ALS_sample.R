
load("../data/ALSFRSdata.rda")

set.seed(12345)

tmp <- ALSFRSdata

# <Field Name="10. Respiratory" ID="1214">
tmp$y <- tmp[["v_1214.halfyearafter"]]
tmp <- tmp[, -grep("halfyearafter", colnames(tmp), ignore.case = TRUE)]
tmp <- tmp[complete.cases(tmp),]

NSim <- 100

learn <- list()
for (i in 1:NSim){
  while(TRUE) {
      learn[[i]] <- sample(1:nrow(tmp), ceiling(nrow(tmp) * .75))
      if (nlevels(tmp$y[learn[[i]], drop = TRUE]) == nlevels(tmp$y))
      break
  }
}

save(tmp, learn, file="ALS_samples.rda")
