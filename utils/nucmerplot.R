# R for plot Nucmer alignments

args <- commandArgs(trailingOnly=T)
nucmer_show <- args[1] # nucmer show-coordi output
outdir <- args[2] # PDF output directory
outpdf <- args[3] # PDF output file

###########################################################
#' module to determine xaxis
###########################################################
smartaxis <- function(maxnum) {
  numdigits <- nchar(maxnum)
  unit <- 10 ^ (numdigits - 1) / (2- round((maxnum / 10 ^ numdigits), 0)) # 1 or 5 e (numdigits - 1)
  subunit <- unit / 5 
  
  numsat <- unit * (0:10)
  numsat <- numsat[numsat < maxnum]
  
  if (numdigits >= 7) {
    numlabels <- numsat / 1000000
    label.scale <- "Mb"
  } else if (numdigits < 7 & numdigits >= 4) {
    numlabels <- numsat / 1000
    label.scale <- "kb"
  } else {
    numlabels <- numsat
    label.scale <- "bp"
  }
  
  subunits <- seq(0, maxnum, by = subunit)
  subunits <- subunits[!subunits %in% c(numsat, 0)] 
  # return
  list(numsat, numlabels, label.scale, subunits)
}

###########################################################
#' bandconnect
###########################################################
bandconnect <- function(topregion, bottomregion, topheight, bottomheight,
                        border=NULL, bandcol="grey") {
  
  topregion <- as.numeric(as.character(topregion))
  bottomregion <- as.numeric(as.character(bottomregion))
  
  ### module
  transform_curve <- function(startp, endp, npoint=1000) {
    ### computer transform_curve coordinates
    midp <- (startp + endp) / 2 
    beizer_value <- sqrt(1:npoint) / sqrt(npoint)  # sqrt as default
    
    curve.x1 <- seq(startp[1], midp[1], by = (midp[1] - startp[1])/(npoint-1))
    curve.y1 <- startp[2] - beizer_value * (startp[2] - midp[2])
    
    curve.x2 <- seq(midp[1], endp[1], by = (endp[1] - midp[1])/(npoint-1))
    curve.y2 <- rev(endp[2] - beizer_value * (endp[2] - midp[2]))
    
    curve.x <- c(curve.x1, curve.x2)
    curve.y <- c(curve.y1, curve.y2)
    list(x=curve.x, y=curve.y)
  }
  
  p1 <- transform_curve(c(topregion[1], topheight), c(bottomregion[1], bottomheight))
  p2 <- transform_curve(c(topregion[2], topheight), c(bottomregion[2], bottomheight))
  px <- c(p1$x, rev(p2$x))
  py <- c(p1$y, rev(p2$y))
  lines(px, py) 
  polygon(px, py, border=border, col=bandcol)
}


###########################################################
# nucmerplot
###########################################################
nucmerplot <- function(datafile, outpath=".", imageoutfile) {
  
  # data of alignments, qry, and subj
  aln <- read.delim(datafile, stringsAsFactor=F)
  qlen <- unique(as.numeric(as.character(aln$qlen)))
  slen <- unique(as.numeric(as.character(aln$slen)))
  qry <- unique(aln$qry)
  subj <- unique(aln$subj)

  ### plot
  outpdffile <- paste0(outpath, "/", imageoutfile)
  pdf(outpdffile, width=5, height=3)
  
  # plot
  par(mar=c(2.2, 0.5, 1, 0.5), mgp=c(1.5,0.1,0))
  
  # x-range
  max.xsize <- max(qlen, slen)
	xrange <- c(-max.xsize / 50, max.xsize + max.xsize / 50)
	
	# connection colors
	high.ld.col <- "deepskyblue4"
	low.ld.col <- "grey96"
	colfunc <- colorRampPalette(c(low.ld.col, high.ld.col))
	color_num <- 20
	allcols <- colfunc(color_num)
	
  # identity
	identity_high <- max(aln$identity)
	identity_low <- min(aln$identity)
	if (identity_high == identity_low |
	    (identity_high - identity_low < 5)) {
	  identity_low <- identity_high - 5
	}
	
	identity20 <- seq(identity_low, identity_high, by=(identity_high - identity_low)/(color_num -1))
	
	# y-range
	yrange <- c(0, 1)

	# plot canvas
  plot(NULL, NULL, type="n", xlim=xrange, ylim=yrange,
       xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
  
  # x-axis
  xaxis <- smartaxis(max.xsize)
  axis(1, at=xaxis[[1]], labels=xaxis[[2]], tick=F, col="gray50", cex.axis=0.8)
  mtext(paste("position (", xaxis[[3]], ")"), side=1, line=1)
  
  # vertical lines
  abline(v=xaxis[[1]], col="gray70", lwd=0.8)
  abline(v=xaxis[[4]], col="gray90", lwd=0.5)
  
  # region lines
  rect(0, 0.845, qlen, 0.865, col="gray70", border="gray50") # 0.85
  rect(0, 0.14, slen, 0.16, col="gray70", border="gray50") # 0.15
  
  # connections
  for (i in 1:nrow(aln)) {
    band_col <- allcols[which.min(abs(identity20 - aln[i, "identity"]))]
    bandconnect(topregion=aln[i, 1:2], topheight=0.835,
                bottomregion=aln[i, 3:4], bottomheight=0.17,
                border=NA, bandcol=band_col)
  }
  
  # text labels
  xlabels <- qry
  ylabels <- subj
  text(x= -max.xsize / 50, y=0.93, labels=xlabels, cex=0.8, xpd=T, pos=4)
  text(x= -max.xsize / 50, y=0.05, labels=ylabels, cex=0.8, xpd=T, pos=4)

  ### identity legends
  barlen <- max.xsize / 5
  barstep <- barlen / color_num
  barheight <- 0.02
  text(max.xsize - barlen/2, 0.975, labels = "identity", cex=0.6, pos=1, xpd=T)  ### plot LD name
  barlabels <- c(identity_low, identity_high)
  barlabels.num <- floor(barlabels * color_num)
  text(max.xsize-barlen, 0.985, labels=identity_low, cex=0.6, pos=1)  ### plot low identity
  text(max.xsize, 0.985, labels=identity_high, cex=0.6, pos=1)  ### plot low identity
  
  # color bars for identity legends
  col_barpos <- max.xsize-barlen
  for (i in 1:color_num) {
    rect(col_barpos, 0.975, col_barpos + barstep, 0.975 + barheight,
         border=NA, col=allcols[i])
    col_barpos <- col_barpos + barstep
  }
  
  # close image output
  dev.off()
}

# execute plotting
nucmerplot(datafile=nucmer_show, outpath=outdir, imageoutfile=outpdf)
