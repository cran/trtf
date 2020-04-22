
rm(list=ls())

library("lattice")
library("latticeExtra")

load("ALS_results_combined.rda")

lwd <- 1.5
col <- rgb(.1, .1, .1, .05)
colR <- rgb(.75, 0, 0, .8)
colRlight <- rgb(.75, 0, 0, .1)
colB <- rgb(0, 0, 0.75, .8)
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

mypanel <- function(x, y, groups, subscripts, ...) {
    # panel.abline(h = 0, lty = 3)
    panel.bwplot(x = x, y = y, ...)
    tapply(1:length(y), groups[subscripts], function(i) {
        llines(x = 1:nlevels(x), y = y[i][order(x[i])], 
               col = rgb(0.1, 0.1, 0.1, 0.1))
    })
}

# only take 100 from 120 random splits into account
res <- res[complete.cases(res),][1:100,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Log-Likelihood ----------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lab <- c('polrconst' = "Polr()", 
         'polr' = "Polr(\U0001D6FC)",
         'of_eq_ll' = 'Ordinal Forest \n(equal)',
         'of_prop_ll' = 'Ordinal Forest \n(proportional)',
         'cf_ll' = 'CForest',
         'tf_alpha_ll' = 'OTF(\U0001D6FC)',
         'tf_theta_ll' = 'OTF(\U0001D6DD)')

i <- 1:nrow(res)

colnam.ll <- c("polrconst", "polr",
               "of_eq_ll", "of_prop_ll",
               "cf_ll", "tf_alpha_ll", "tf_theta_ll")
res.ll <- res[, colnam.ll]

m <- colnames(res.ll)
r <- expand.grid(i = i, m = m)
r$m <- factor(r$m, levels = names(lab), label = lab)
r$fll <- do.call("c", res.ll)
r$id <- factor(i)

ll.ALS.100.250.plt <- bwplot(fll ~ m, data = r, groups = id, 
                 scales = list(x = list(rot = 45,
                                        label = lab)),
                 main = "ALS SimStudy with ntree = 250",
                 panel = mypanel, ylab = "Log-likelihood")

### Figure 2
print(ll.ALS.100.250.plt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Quadratic weighted Cohen's Kappa ----------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lab <- c('tf_of_qk' = "Ordinal Forest (proportional) vs.\n OTF(\U0001D6DD)")

colnam.qwck <- c("tf_of_qk")
res.qwck <- data.frame(qwck = (res[ , colnam.qwck]))
res.qwck$lab <- factor(lab)
qwck.ALS.100.250.plt <- bwplot(qwck ~ lab, data = res.qwck, 
                               main = " ",
                               ylab = "Quadratic weighted\n Cohen's Kappa")

### Figure 3
print(qwck.ALS.100.250.plt)

