
x <- readLines("ALS_analysis.R")

idx <- 1:120

if(!dir.exists("rda")){dir.create("rda")}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Calculate individiual results.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in idx){
  f <- paste("rda/run_", i, ".R", sep = "")
  writeLines(gsub("INDEX", i, x), con = f)
  system(paste("R CMD BATCH", f, ifelse(i < max(idx), " & ", "")))
  print(i)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Combine individiual results.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out <- 0

for (i in idx) {
  
  f <- paste("rda/ALS_results_combined_", i, ".rda", sep = "")
  
  if (file.exists(f)) {
    load(f)
    out <- out + res
  } else {
    cat(f, "missing \n")
  }
  
}

res <- out

res[res == 0] <- NA

save(res, file = paste("ALS_results_combined.rda", sep = ""))
