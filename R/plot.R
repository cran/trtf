
### very similar to partykit:::node_barplot, maybe we can modularize this
### better?
node_mlt_disc <- function(obj, newdata = data.frame(1), col = "black", bg = "white", fill = "transparent",
                     ylines = 2, id = TRUE, mainlab = NULL, gp = gpar(), 
                     yscale = NULL, beside = NULL, ymax = NULL, widths = 1,
                     gap = NULL, reverse = NULL, rot = 0,
                     just = c("center", "top"), 
                     text = c("none", "horizontal", "vertical"),
                     K = 20, type = c("trafo", "distribution", "survivor",  
                              "density", "logdensity", "hazard",    
                              "loghazard", "cumhazard", "quantile"),
                     flip = FALSE, axes = TRUE, ...)
{

    mod <- obj$model
    q <- mkgrid(mod, n = K)[[mod$response]]
    type <- match.arg(type)

    if (type %in% c("distribution", "survivor")) {
        yscale <- c(0, 1)
    } else {
        if (is.null(yscale)) {
            probs <- predict.trafotree(obj, q = q, type = type, ...)
            yscale <- range(probs, finite = TRUE)
        }
    }

    stopifnot(is.factor(q) || isTRUE(all.equal(round(q), q)))
    
    if(is.factor(q)) {
        ylevels <- levels(q)
	if(is.null(beside)) beside <- if(length(ylevels) < 3L) FALSE else TRUE
        if(is.null(ymax)) ymax <- if(beside) 1.1 else 1
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    } else {
        if(is.null(beside)) beside <- TRUE
        if(is.null(ymax)) ymax <- if(beside) max(probs) * 1.1 else max(probs)
        ylevels <- rownames(probs)
        if(length(ylevels) < 2) ylevels <- ""
	if(is.null(gap)) gap <- if(beside) 0.1 else 0
    }
    if(is.null(reverse)) reverse <- !beside
    if(is.null(fill)) fill <- gray.colors(length(ylevels))
    if(is.null(ylines)) ylines <- if(beside) c(3, 2) else c(1.5, 2.5)

    ## text labels?
    if(isTRUE(text)) text <- "horizontal"
    if(!is.character(text)) text <- "none"
    text <- match.arg(text, c("none", "horizontal", "vertical"))

    ### panel function for barplots in nodes
    rval <- function(node) {
    
        ## id
	nid <- id_node(node)
    
        ## parameter setup
        dat <- data_party(obj, nid)
        wn <- dat[["(weights)"]]   

        cf <- obj$coef[as.character(nid),]
        coef(mod) <- cf
        probs <- predict(mod, newdata = newdata, q = q, type = type)
        if (!is.matrix(probs)) probs <- matrix(probs, ncol = 1)
        if (length(col) != ncol(probs)) col <- rep(col, length.out = ncol(probs))

        pred <- probs
	if(reverse) {
	  pred <- rev(pred)
	  ylevels <- rev(ylevels)
	}
        np <- length(pred)
	nc <- if(beside) np else 1

	fill <- rep(fill, length.out = np)	
        widths <- rep(widths, length.out = nc)
	col <- rep(col, length.out = nc)
	ylines <- rep(ylines, length.out = 2)

	gap <- gap * sum(widths)
        yscale <- c(0, ymax)
        xscale <- c(0, sum(widths) + (nc+1)*gap)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines[1], 1, ylines[2]), c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"), 
                           height = unit(1, "npc") - unit(2, "lines"),
			   name = paste0("node_barplot", nid),
			   gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {	
	  mainlab <- if(id) {
	    function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
  	  } else {
	    function(id, nobs) sprintf("n = %s", nobs)
	  }
        }
	if (is.function(mainlab)) {
          mainlab <- mainlab(names(obj)[nid], sum(wn))
	}
        grid.text(mainlab)
        popViewport()
	
        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale=xscale, yscale=yscale,
			 name = paste0("node_barplot", node$nodeID, "plot"),
			 clip = FALSE)

        pushViewport(plot)
	
	if(beside) {
  	  xcenter <- cumsum(widths+gap) - widths/2
          if(length(xcenter) > 1) grid.xaxis(at = xcenter, label = FALSE)
	  grid.text(ylevels, x = xcenter, y = unit(-1, "lines"), 
                    just = just, rot = rot,
	            default.units = "native", check.overlap = TRUE)
          grid.yaxis()
          grid.rect(gp = gpar(fill = "transparent"))
	  grid.clip()
	  for (i in 1:np) {
            grid.rect(x = xcenter[i], y = 0, height = pred[i], 
                      width = widths[i],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
            if(text != "none") {
              grid.text(x = xcenter[i], y = pred[i] + 0.025,
	        label = paste(format(round(100 * pred[i], 1), nsmall = 1), "%", sep = ""),
	        just = if(text == "horizontal") c("center", "bottom") else c("left", "center"),
	        rot = if(text == "horizontal") 0 else 90,
		default.units = "native")
            }
	  }
	} else {
  	  ycenter <- cumsum(pred) - pred

          if(np > 1) {
	    grid.text(ylevels[1], x = unit(-1, "lines"), y = 0,
                      just = c("left", "center"), rot = 90,
	              default.units = "native", check.overlap = TRUE)
	    grid.text(ylevels[np], x = unit(-1, "lines"), y = ymax,
                      just = c("right", "center"), rot = 90,
	              default.units = "native", check.overlap = TRUE)
	  }
          if(np > 2) {
	    grid.text(ylevels[-c(1,np)], x = unit(-1, "lines"), y = ycenter[-c(1,np)],
                      just = "center", rot = 90,
	              default.units = "native", check.overlap = TRUE)
	  }
          grid.yaxis(main = FALSE)	

          grid.clip()
          grid.rect(gp = gpar(fill = "transparent"))
	  for (i in 1:np) {
            grid.rect(x = xscale[2]/2, y = ycenter[i], height = min(pred[i], ymax - ycenter[i]), 
                      width = widths[1],
	              just = c("center", "bottom"), default.units = "native",
	              gp = gpar(col = col[i], fill = fill[i]))
	  }
	}
	grid.rect(gp = gpar(fill = "transparent"))

	
        upViewport(2)
    }
    
    return(rval)
}
class(node_mlt_disc) <- "grapcon_generator"


node_mlt <- function(obj, newdata = data.frame(1), col = "black", bg = "white", fill = "transparent",
                     ylines = 2, id = TRUE, mainlab = NULL, gp = gpar(), K = 20,
                     yscale = NULL,
                     type = c("trafo", "distribution", "survivor",  
                              "density", "logdensity", "hazard",    
                              "loghazard", "cumhazard", "quantile"),
                     flip = FALSE, axes = TRUE, ...)
{
    mod <- obj$model
    q <- mkgrid(mod, n = K)[[mod$response]]

    if (is.integer(q) || is.factor(q))
        return(node_mlt_disc(obj = obj, newdata = newdata, col = col, bg = bg, fill = fill,
                      ylines = ylines, id = id, mainlab = mainlab, gp = gp, K = K,
                      type = type, flip = flip, axes = axes, ...))

    type <- match.arg(type)

    if (type %in% c("distribution", "survivor")) {
        yscale <- c(0, 1)
    } else {
        if (is.null(yscale)) {
            probs <- predict.trafotree(obj, q = q, type = type, ...)
            yscale <- range(probs, finite = TRUE)
        }
    }

    xscale <- range(q)
    axes <- rep_len(axes, 2)

    ### panel function for ecdf in nodes
    rval <- function(node) {

        nid <- id_node(node)
        dat <- data_party(obj, nid)
        wn <- dat[["(weights)"]]   

        cf <- obj$coef[as.character(nid),]
        coef(mod) <- cf
        y <- predict(mod, newdata = newdata, q = q, type = type)
        if (!is.matrix(y)) y <- matrix(y, ncol = 1)
        if (length(col) != ncol(y)) col <- rep(col, length.out = ncol(y))

        ## set up plot
        q <- q - xscale[1]
        q <- q / diff(xscale)
        y <- y - yscale[1]
        y <- y / diff(yscale)

        top_vp <- viewport(layout = grid.layout(nrow = 2, ncol = 3,
                           widths = unit(c(ylines, 1, 1),
                                         c("lines", "null", "lines")),
                           heights = unit(c(1, 1), c("lines", "null"))),
                           width = unit(1, "npc"),
                           height = unit(1, "npc") - unit(2, "lines"),
                           name = paste("node_mlt", nid, sep = ""), gp = gp)

        pushViewport(top_vp)
        grid.rect(gp = gpar(fill = bg, col = 0))

        ## main title
        top <- viewport(layout.pos.col=2, layout.pos.row=1)
        pushViewport(top)
        if (is.null(mainlab)) {
          mainlab <- if(id) {  
            function(id, nobs) sprintf("Node %s (n = %s)", id, nobs)
          } else {
            function(id, nobs) sprintf("n = %s", nobs)
          }
        }  
        if (is.function(mainlab)) {
          mainlab <- mainlab(nid, sum(wn))
        }
        grid.text(mainlab)
        popViewport()

        plot <- viewport(layout.pos.col=2, layout.pos.row=2,
                         xscale = if(flip) yscale else xscale,
			 yscale = if(flip) xscale else yscale,
                         name = paste0("node_mlt", nid, "plot"),
                         clip = FALSE)

        pushViewport(plot)        
        if(axes[1]) grid.xaxis()
        if(axes[2]) grid.yaxis()
        grid.rect(gp = gpar(fill = "transparent"))
        grid.clip()
        draw <- function(i) {
            y <- y[,i]
            col <- col[i]
            if(flip) {
                if(fill != "transparent") {
                    grid.polygon(c(min(yscale), y, min(yscale)), c(q[1], q, q[K]), gp = gpar(col = col, fill = fill))
                } else {
                    grid.lines(y, q, gp = gpar(col = col))
                }
           } else {
               if(fill != "transparent") {
                   grid.polygon(c(q[1], q, q[K]), c(min(yscale), y, min(yscale)), gp = gpar(col = col, fill = fill))
               } else {
                   grid.lines(q, y, gp = gpar(col = col))
               }
           }
        }
        out <- sapply(1:ncol(y), draw)
        upViewport(2)
    }

    return(rval)
}
class(node_mlt) <- "grapcon_generator"

plot.trafotree <- function(x, newdata = data.frame(1), ...) {
    class(x) <- class(x)[-1L]
    tp <- function(...) node_mlt(newdata = newdata, ...)
    class(tp) <- class(node_mlt)
    plot(x, terminal_panel = tp, ...)
}
