
x <- readLines("03_ALS_analysis.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Calculate individiual results.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

idx <- 1:100

for (i in idx){
    f <- paste("run_", i, ".R", sep = "")
    writeLines(gsub("INDEX", i, x), con = f)
    system(paste("R CMD BATCH", f, ifelse(i < max(idx), " & ", "")))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Combine individiual results.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out <- 0

for (i in idx) {
  f <- paste("ALS_results_combined_", i, ".rda", sep = "")
  if (file.exists(f)) {
    load(f)
    out <- out + ret
  } else {
    cat(f, "missing \n")
  }
}

ret <- out

ret[ret == 0] <- NA

save(ret, file = "ALS_results_combined.rda")
