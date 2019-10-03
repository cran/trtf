
library("lattice")
library("latticeExtra")

load("../simulation/results_empeval.rda")

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

lab <- c('of' = 'Ordinal Forest \n(equal)',
         'of_prop' = 'Ordinal Forest \n(proportional)',
         'cf' = 'CForest',
         'tf_B_alpha' = 'OTF(\U0001D6FC)', 
         'tf_B' = 'OTF(\U0001D6DD)')

i <- 1:nrow(ret)
nm <- c("p", "prod_shift", "prod_scale", "loglik")
m <- colnames(ret)
m <- m[!(m %in% nm)]

x <- expand.grid(i = i, m = m)
r <- ret[x$i,nm]

r$m <- factor(x$m, levels = names(lab), label = lab)
r$fll <- do.call("c", ret[,m])
r$id <- factor(i)
for (n in nm[1:3]) r[[n]] <- factor(r[[n]])

r$setup <- with(r, interaction(prod_shift, prod_scale))
levels(r$setup) <- lev <- c("No", "PO", "Non-PO", "Combined")

r$setup <- factor(as.character(r$setup), levels = lev, labels = lev)
levels(r$p) <- c("Low-dimensional", "High-dimensional")

r$fll[r$fll == 0] <- NA

plt <- bwplot(I(fll - loglik) ~ m | p + setup, data = r, groups = id,
              scales = list(x = list(rot = 45, label = lab)),
              layout = c(2, 4),
              panel = mypanel, as.table = TRUE,
              ylab = "Log-likelihood Difference")
useOuterStrips(plt)
