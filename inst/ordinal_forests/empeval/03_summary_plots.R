
rm(list=ls())

source("competitors.R")

library("lattice")
library("latticeExtra")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### GLobal Variables --------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lwd <- 1.5
col <- rgb(.1, .1, .1, .05)
colR <- rgb(.75, 0, 0, .8)
colRlight <- rgb(.75, 0, 0, .1)
colB <- rgb(0, 0, 0.75, .8)
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme

ltheme$layout.heights$strip <- 2.5
ltheme$layout.heights$axis.panel <- 2.5
ltheme$layout.heights$panel <- 2.5
ltheme$par.sub.text$lineheight <- 2
ltheme$par.main.text$cex <- 1
ltheme$par.main.text$lineheight <- 2.5

ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

LL.KL.mypanel <- function(x, y, groups, subscripts, ...){
    panel.abline(h = 0, lty = 3)
    panel.bwplot(x = x, y = y, ...)
    tapply(1:length(y), groups[subscripts], function(i) {
        llines(x = 1:nlevels(x), y = y[i][order(x[i])], 
               col = rgb(0.1, 0.1, 0.1, 0.1))
    })
}

plot.bwplot <- function(data.tmp, lab.tmp, ylab.tmp,
                        layout.tmp, mypanel,
                        title.tmp = rownames(sim_para)[p]){
    
    i <- 1:nrow(data.tmp)
    nm <- colnames(data.tmp)[param.colnam]
    m <- colnames(data.tmp)[-param.colnam]
    # m <- m[!(m %in% nm)]
    
    x <- expand.grid(i = i, m = m)
    r <- data.tmp[x$i, nm]
    
    r$m <- factor(x$m, levels = names(lab.tmp), label = lab.tmp)
    r$y <- do.call("c", data.tmp[,m])
    
    r$id <- factor(i)
    for (n in nm) r[[n]] <- factor(r[[n]])
    
    r$setup <- with(r, interaction(prod_shift, prod_scale))
    levels(r$setup) <- lev <- c("No", "PO", "Non-PO", "Combined")
    r$setup <- factor(as.character(r$setup), levels = lev, labels = lev)
    levels(r$p) <- c("Low-dimensional", "High-dimensional")
    
    plt <- bwplot(y ~ m | p + setup, data = r, groups = id,
                  between = (list(x = 0.25)),
                  par.settings = list(layout.heights = list(strip = 2.5)), 
                  scales = list(x = list(rot = 45, label = lab.tmp)),
                  as.table = TRUE,
                  layout = layout.tmp,
                  main = title.tmp,
                  panel = mypanel,
                  ylab = ylab.tmp)
    
    return(plt)
}

plot.ll.diff.bwplot <- function(data.tmp, lab.tmp, ylab.tmp,
                                layout.tmp, mypanel,
                                title.tmp = rownames(sim_para)[p]){
    
    i <- 1:nrow(data.tmp)
    nm <- c("p", "prod_shift", "prod_scale", "loglik")
    m <- names(lab.tmp)
    
    x <- expand.grid(i = i, m = m)
    r <- data.tmp[x$i, nm]
    
    r$m <- factor(x$m, levels = names(lab.tmp), label = lab.tmp)
    r$y <- do.call("c", data.tmp[,m])
    
    r$id <- factor(i)
    for (n in  c("p", "prod_shift", "prod_scale")) r[[n]] <- factor(r[[n]])
    
    r$setup <- with(r, interaction(prod_shift, prod_scale))
    levels(r$setup) <- lev <- c("No", "PO", "Non-PO", "Combined")
    r$setup <- factor(as.character(r$setup), levels = lev, labels = lev)
    levels(r$p) <- c("Low-dimensional", "High-dimensional")
    
    r$y[r$y == 0] <- NA
    
    plt <- bwplot(I(y - loglik) ~ m | p + setup, data = r, groups = id,
                  between = (list(x = 0.25)), #as.table = TRUE,
                  par.settings = list(layout.heights = list(strip = 2.5)), 
                  scales = list(x = list(rot = 45, label = lab.tmp)),
                  as.table = TRUE,
                  layout = layout.tmp,
                  main = title.tmp,
                  panel = mypanel,
                  ylab = ylab.tmp)
    
    return(plt)
}

load("../simulation/rda/results_empeval_nsim100_ntree250.rda")
res.100.250 <- res
load("../simulation/rda/results_empeval_nsim100_ntree2000.rda")
res.100.2000 <- res

param.colnam <- which(colnames(res) %in% c("p",
                                           "prod_shift",
                                           "prod_scale"))
ll.colnam <- grep("loglik", colnames(res))
ll.colnam <- c(ll.colnam, grep("_ll", colnames(res)))
ck.colnam <- grep("_ck", colnames(res))
lk.colnam <- grep("_lk", colnames(res))
qk.colnam <- grep("_qk", colnames(res))
KLI.colnam <- which(grepl("_KL_I\\b", colnames(res)))
KLII.colnam <- which(grepl("_KL_II\\b", colnames(res)))
KLIII.colnam <- which(grepl("_KL_III\\b", colnames(res)))

### nsim = 100, ntree = 250
res.100.250.ll <- res.100.250[, c(param.colnam, ll.colnam)]
res.100.250.ck <- res.100.250[, c(param.colnam, ck.colnam, lk.colnam, qk.colnam)]
res.100.250.KLI <- res.100.250[, c(param.colnam, KLI.colnam)]
res.100.250.KLII <- res.100.250[, c(param.colnam, KLII.colnam)]
res.100.250.KLIII <- res.100.250[, c(param.colnam, KLIII.colnam)]

### nsim = 100, ntree = 2000
res.100.2000.ll <- res.100.2000[, c(param.colnam, ll.colnam)]
res.100.2000.ck <- res.100.2000[, c(param.colnam, ck.colnam, lk.colnam, qk.colnam)]
res.100.2000.KLI <- res.100.2000[, c(param.colnam, KLI.colnam)]
res.100.2000.KLII <- res.100.2000[, c(param.colnam, KLII.colnam)]
res.100.2000.KLIII <- res.100.2000[, c(param.colnam, KLIII.colnam)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Log-Likelihood Differences ----------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ll.diff.100.250.plt <- plot.ll.diff.bwplot(data.tmp = res.100.250.ll,
                                           lab.tmp = c('of_eq_ll' = 'Ordinal Forest \n(equal)',
                                                       'of_prop_ll' = 'Ordinal Forest \n(proportional)',
                                                       'cf_ll' = 'CForest',
                                                       'tf_alpha_ll' = 'OTF(\U0001D6FC)', 
                                                       'tf_theta_ll' = 'OTF(\U0001D6DD)'),
                                           ylab.tmp = "Log-Likelihood Differences",
                                           layout.tmp = c(2, 4),
                                           mypanel = LL.KL.mypanel,
                                           title.tmp = NULL)
### Figure 1
print(useOuterStrips(ll.diff.100.250.plt))


ll.diff.100.2000.plt <- plot.ll.diff.bwplot(data.tmp = res.100.2000.ll,
                                            lab.tmp = c('of_eq_ll' = 'Ordinal Forest \n(equal)',
                                                        'of_prop_ll' = 'Ordinal Forest \n(proportional)',
                                                        'cf_ll' = 'CForest',
                                                        'tf_alpha_ll' = 'OTF(\U0001D6FC)', 
                                                        'tf_theta_ll' = 'OTF(\U0001D6DD)'),
                                            ylab.tmp = "Log-Likelihood Differences",
                                            layout.tmp = c(2, 4),
                                            mypanel = LL.KL.mypanel,
                                            title.tmp = NULL)
### Figure 4
print(useOuterStrips(ll.diff.100.2000.plt))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Log-Likelihood ----------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ll.100.250.plt <- plot.bwplot(data.tmp = res.100.250.ll,
                              lab.tmp = c('loglik' = 'True Log-Likelihood',
                                          'of_eq_ll' = 'Ordinal Forest \n(equal)',
                                          'of_prop_ll' = 'Ordinal Forest \n(proportional)',
                                          'cf_ll' = 'CForest',
                                          'tf_alpha_ll' = 'OTF(\U0001D6FC)', 
                                          'tf_theta_ll' = 'OTF(\U0001D6DD)'),
                              ylab.tmp = "Log-Likelihood",
                              layout.tmp = c(2, 4),
                              mypanel = LL.KL.mypanel,
                              title.tmp = NULL)
### Figure 5
print(useOuterStrips(ll.100.250.plt))

ll.100.2000.plt <- plot.bwplot(data.tmp = res.100.2000.ll,
                               lab.tmp = c('loglik' = 'True Log-Likelihood',
                                           'of_eq_ll' = 'Ordinal Forest \n(equal)',
                                           'of_prop_ll' = 'Ordinal Forest \n(proportional)',
                                           'cf_ll' = 'CForest',
                                           'tf_alpha_ll' = 'OTF(\U0001D6FC)', 
                                           'tf_theta_ll' = 'OTF(\U0001D6DD)'),
                               ylab.tmp = "Log-Likelihood",
                               layout.tmp = c(2, 4),
                               mypanel = LL.KL.mypanel,
                               title.tmp = NULL)
### Figure 6
print(useOuterStrips(ll.100.2000.plt))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Kullback-Leibler Divergence ---------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kullback-Leibler divergence KL comparing the true and estimated conditional
# probability densities according to Equation (7)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KLIII.100.250.plt <- plot.bwplot(data.tmp = res.100.250.KLIII * 500,
                                 lab.tmp = c('of_eq_KL_III' = 'Ordinal Forest \n(equal)',
                                             'of_prop_KL_III' = 'Ordinal Forest \n(proportional)',
                                             'cf_KL_III' = 'CForest',
                                             'tf_alpha_KL_III' = 'OTF(\U0001D6FC)', 
                                             'tf_theta_KL_III' = 'OTF(\U0001D6DD)'),
                                 ylab.tmp = "Kullback-Leibler Divergence KL",
                                 layout.tmp = c(2, 4),
                                 mypanel = LL.KL.mypanel,
                                 title.tmp = NULL)
### Figure 7
print(useOuterStrips(KLIII.100.250.plt))

KLIII.100.2000.plt <- plot.bwplot(data.tmp = res.100.2000.KLIII * 500,
                                  lab.tmp = c('of_eq_KL_III' = 'Ordinal Forest \n(equal)',
                                              'of_prop_KL_III' = 'Ordinal Forest \n(proportional)',
                                              'cf_KL_III' = 'CForest',
                                              'tf_alpha_KL_III' = 'OTF(\U0001D6FC)', 
                                              'tf_theta_KL_III' = 'OTF(\U0001D6DD)'),
                                  ylab.tmp = "Kullback-Leibler Divergence KL",
                                  layout.tmp = c(2, 4),
                                  mypanel = LL.KL.mypanel,
                                  title.tmp = NULL)
### Figure 8
print(useOuterStrips(KLIII.100.2000.plt))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Results of the Kullback-Leibler divergence KL_1 according to Equation (8).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KLI.100.250.plt <- plot.bwplot(data.tmp = res.100.250.KLI * 500,
                               lab.tmp = c('of_eq_KL_I' = 'Ordinal Forest \n(equal)',
                                           'of_prop_KL_I' = 'Ordinal Forest \n(proportional)',
                                           'cf_KL_I' = 'CForest',
                                           'tf_alpha_KL_I' = 'OTF(\U0001D6FC)', 
                                           'tf_theta_KL_I' = 'OTF(\U0001D6DD)'),
                               ylab.tmp = expression(paste("Kullback-Leibler Divergence ", KL[1])),
                               layout.tmp = c(2, 4),
                               mypanel = LL.KL.mypanel,
                               title.tmp = NULL)
### Figure 9
print(useOuterStrips(KLI.100.250.plt))

KLI.100.2000.plt <- plot.bwplot(data.tmp = res.100.2000.KLI * 500,
                                lab.tmp = c('of_eq_KL_I' = 'Ordinal Forest \n(equal)',
                                            'of_prop_KL_I' = 'Ordinal Forest \n(proportional)',
                                            'cf_KL_I' = 'CForest',
                                            'tf_alpha_KL_I' = 'OTF(\U0001D6FC)', 
                                            'tf_theta_KL_I' = 'OTF(\U0001D6DD)'),
                                ylab.tmp = expression(paste("Kullback-Leibler Divergence ", KL[1])),
                                layout.tmp = c(2, 4),
                                mypanel = LL.KL.mypanel,
                                title.tmp = NULL)
### Figure 10
print(useOuterStrips(KLI.100.2000.plt))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Results of the Kullback-Leibler divergence KL_2 according to Equation (9).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KLII.100.250.plt <- plot.bwplot(data.tmp = res.100.250.KLII * 500,
                                lab.tmp = c('of_eq_KL_II' = 'Ordinal Forest \n(equal)',
                                            'of_prop_KL_II' = 'Ordinal Forest \n(proportional)',
                                            'cf_KL_II' = 'CForest',
                                            'tf_alpha_KL_II' = 'OTF(\U0001D6FC)', 
                                            'tf_theta_KL_II' = 'OTF(\U0001D6DD)'),
                                ylab.tmp = expression(paste("Kullback-Leibler Divergence ", KL[2])),
                                layout.tmp = c(2, 4),
                                mypanel = LL.KL.mypanel,
                                title.tmp = NULL)
### Figure 11
print(useOuterStrips(KLII.100.250.plt)) 

KLII.100.2000.plt <- plot.bwplot(data.tmp = res.100.2000.KLII * 500,
                                 lab.tmp = c('of_eq_KL_II' = 'Ordinal Forest \n(equal)',
                                             'of_prop_KL_II' = 'Ordinal Forest \n(proportional)',
                                             'cf_KL_II' = 'CForest',
                                             'tf_alpha_KL_II' = 'OTF(\U0001D6FC)', 
                                             'tf_theta_KL_II' = 'OTF(\U0001D6DD)'),
                                 ylab.tmp = expression(paste("Kullback-Leibler Divergence ", KL[2])),
                                 layout.tmp = c(2, 4),
                                 mypanel = LL.KL.mypanel,
                                 title.tmp = NULL)
### Figure 12
print(useOuterStrips(KLII.100.2000.plt))

