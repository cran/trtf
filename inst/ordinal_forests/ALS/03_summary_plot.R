
library("lattice")
library("latticeExtra")

load('examples/ALS_results_combined.rda') 

lwd <- 1.5
col <- rgb(.1, .1, .1, .05)
colR <- rgb(.75, 0, 0, .8)
colRlight <- rgb(.75, 0, 0, .1)
colB <- rgb(0, 0, 0.75, .8)
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

mypanel <- function(x, y, groups, subscripts, ...) {
  panel.abline(h = 0, lty = 3)
  panel.bwplot(x = x, y = y, ...)
  tapply(1:length(y), groups[subscripts], function(i) {
    llines(x = 1:nlevels(x), y = y[i][order(x[i])], 
           col = rgb(0.1, 0.1, 0.1, 0.1))
  })
}

lab <- c('polrconst' = "Polr()", 
         'polr' = "Polr(\U0001D6FC)", 
         'of' = 'Ordinal Forest \n(equal)',
         'of_prop' = 'Ordinal Forest \n(proportional)',
         'cf' = 'CForest',
         'tf_B_alpha' = 'OTF(\U0001D6FC)', 
         'tf_B' = 'OTF(\U0001D6DD)')

i <- 1:nrow(ret)
m <- colnames(ret)
r <- expand.grid(i = i, m = m)
r$m <- factor(r$m, levels = names(lab), label = lab)
r$fll <- do.call("c", ret)
r$id <- factor(i)

plt <- bwplot(fll ~ m, data = r, groups = id, 
              scales = list(x = list(rot = 45, #relation = "free",
                                     label = lab)),
              panel = mypanel, ylab = "Log-likelihood")
plt