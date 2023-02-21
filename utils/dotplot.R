# R for plot Nucmer alignments
args <- commandArgs(trailingOnly=T)
nucmer_show <- args[1] # nucmer show-coordi output
line_color <- args[2] # color to be used for plotting
outdir <- args[3] # PDF output directory
outpdf <- args[4] # PDF output file

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
#' check if a color is valid
#' https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
###########################################################
isColor <- function(col) {
  sapply(col, function(incol) {
    tryCatch(is.matrix(col2rgb(incol)), 
             error = function(e) FALSE)
  })
}

###########################################################
#' dotplot
###########################################################
nucmer.dotplot <- function(datafile, line_col="deepskyblue4",
                           outpath=".", imageoutfile,
                           lend.turnoff=T, line.width.factor=5) {
  # check if the color is valid
  if (!isColor(line_col)) {
    cat(line_col, "is not valid. Plot using the default\n")
    line_col="deepskyblue4"
  }
  
  ### datafile is the show-coord output
  aln <- read.delim(datafile, stringsAsFactor=F)
  stopifnot(nrow(aln)>=1)
  
  ### basic information
  qlen <- unique(as.numeric(as.character(aln$qlen)))
  slen <- unique(as.numeric(as.character(aln$slen)))
  qry <- unique(aln$qry)
  stopifnot(length(qry)>=1)
  subj <- unique(aln$subj)
  stopifnot(length(subj)>=1)
  
  ### pdf output
  outpdffile <- paste0(outpath, "/", imageoutfile)
  pdf(outpdffile, width=5, height=5)
  
  # plot
  par(mar=c(2.5, 2.5, 1, 1), mgp=c(1.2,0.1,0))
  
  # x-range
  xrange <- c(0, qlen)
  yrange <- c(0, slen)
  max.size <- max(qlen, slen)
  
  ### plot
  plot(NULL, NULL, type="n", xlim=xrange, ylim=yrange,
       xlab="", ylab="", xaxt="n", yaxt="n")
  
  # x-axis
  xaxis <- smartaxis(qlen)
  axis(1, at=xaxis[[1]], labels=xaxis[[2]], tick=F, col="gray50", cex.axis=0.8)
  mtext(paste0(qry, " (", xaxis[[3]], ")"), side=1, line=1.25)
  # vertical lines
  abline(v=xaxis[[1]], col="gray70", lwd=0.8)
  abline(v=xaxis[[4]], col="gray90", lwd=0.5)
  
  # y-axis
  yaxis <- smartaxis(slen)
  axis(2, at=yaxis[[1]], labels=yaxis[[2]], tick=F, col="gray50", cex.axis=0.8)
  mtext(paste0(subj, " (", yaxis[[3]], ")"), side=2, line=1.25)
  # horizontal lines
  abline(h=yaxis[[1]], col="gray70", lwd=0.8)
  abline(h=yaxis[[4]], col="gray90", lwd=0.5)
  
  col1 <- rgb(0, 0.4, 0, 0.5)
  col2 <- rgb(1, 0, 0, 0.5)
  
  for (i in 1:nrow(aln)) {
    plot.col <- col1
    distance <- abs(aln[i, "qend"] - aln[i, "qstart"])
    lend.val <- 2
    
    if (distance > max.size/50) {
      lend.val <- 1
    }
    
    if (lend.turnoff) {
      lend.val <- 1
    }
    
    # dot_lines
    lines(c(aln[i, "qstart"], aln[i, "qend"]),
          c(aln[i, "sstart"], aln[i, "send"]),
          lwd=line.width.factor,
          col= line_col, lend=lend.val)
  }
  
  # close canvas
  dev.off()
}

### run doplot
nucmer.dotplot(datafile=nucmer_show, line_col=line_color, outpath=outdir, imageoutfile=outpdf)
