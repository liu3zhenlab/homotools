# R for plot Nucmer alignments
args <- commandArgs(trailingOnly=T)
blocktype <- args[1] # blocktying input
genebed <- args[2] # gene BED
title <- args[3] # figure title
outpdf <- args[4] # PDF output file

blocktype <- "hgout/hgout.3.blocktyping"
genebed <- "../0_gtfbed/GRMZM2G171650_T01.adjusted.bed"
title <- "figure"
outpdf <- "out.pdf"

#setwd("/homes/liu3zhen/scripts2/homotools/homograph_dev/test")

#############################################################
# module
#############################################################
block_plot <- function(data, yrange, highcol="red", lowcol="grey80", singlecol=NULL) {
  # plot blocks
  # data include xrange, #poly, #diff
  xrange <- data[1:2]
  if (length(data)>2) {
    if (data[3] > 0) {
      numpoly <- data[3]
      numdiff <- data[4]
      colfunc <- colorRampPalette(c(lowcol, highcol))
      color_num <- 11
      allcols <- colfunc(color_num)
      names(allcols) <- round(seq(0, 1, by=1/(color_num -1)),1)
      
      if (!is.null(singlecol)) {
        sel_col <- singlecol
      } else {
        sel_col <- as.character(allcols[as.character(round(numdiff/numpoly, 1))])
      } 
      rect(xleft=xrange[1], ybottom=yrange[1],
           xright=xrange[2], ytop=yrange[2],
           col=sel_col, border=NA)
    }
  } else {
    if (is.null(singlecol)) {
      sel_col <- "grey"
    } else {
      sel_col <- singlecol
    }
    rect(xleft=xrange[1], ybottom=yrange[1],
         xright=xrange[2], ytop=yrange[2],
         col=sel_col, border=NA)
  }
}
#############################################################

#############################################################
# module to determine xaxis
#############################################################
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
#############################################################

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

#############################################################
# plotting
#############################################################
#d <- read.delim("hgout/hgout.3.blocktyping", header=F)
d <- read.delim(blocktype, header=F)
taxa_names_maxlen <- max(nchar(d[, 1]))
left_margin <- round(taxa_names_maxlen / 5) + 1
xmax <- max(d[, 3])
taxa <- unique(d[, 1])
ymax <- length(taxa)
if (genebed != "NULL") {
  ymax <- ymax + 4
}
yheight <- 2 + (ymax - 2) * 0.2

# plot setting
# par(mar=c(2.2, 0.5, 1, 0.5), mgp=c(1.5,0.1,0))
pdf(outpdf, width=5, height=yheight)
par(mar=c(left_margin, 3, 2.5, 0))
plot(NULL, NULL, xlim=c(0, xmax), ylim=c(1, ymax),
     frame.plot=F, xlab="", ylab="",
     xaxt="n", yaxt="n")

# plot gene structure
if (genebed != "NULL") {
  bed <- read.delim(genebed, header=F, comment.char="#")
  bed_min <- min(bed[, 2])
  bed_max <- max(bed[, 3])
  adjust_height <- height*max(1, ymax/20)
  x_adj <- xmax/2 
  ycenter <- ymax
  lines(c(bed_min, bed_max), rep(ycenter, 2), col="gray50", lwd=1)
  if (nrow(bed)>0) {
    for (i in 1:nrow(bed)) {
      color <- bed[i, 7]
      if (!isColor(color)) {
        color <- "grey"
      }
      height <- bed[i, 5]
      rect(bed[i, 2] + 1, ycenter-adjust_height/2,
           bed[i, 3], ycenter+adjust_height/2,
           col=color, border=color)
    }
  }  
}

# x-axis
xaxis <- smartaxis(xmax)
axis(1, at=xaxis[[1]], labels=xaxis[[2]], tick=T, col="gray50", cex.axis=0.8)
mtext(paste0("position (", xaxis[[3]], ")"), side=1, line=2)

# title
mtext(title, side=3, line=1.2)


# backgrounp of alignment regions (alignment segments)
row <- 0
for (etaxon in taxa) {
  alnseg <- paste(d[, 1], d[, 2], d[, 3])
  taxon_data <- d[!duplicated(alnseg) & d[, 1] == etaxon, 1:3] 
  row <- row + 1
  apply(taxon_data[, 2:3], 1, block_plot, yrange=c(row-0.35, row+0.35), singlecol="lightblue")
  text(-xmax/50, row, adj=1, labels=etaxon, xpd=T, cex=0.6)
}

# block and color-coded difference from the consensus
row <- 0
for (etaxon in taxa) {
  block_taxon_data <- d[d[, 1] == etaxon, 4:7] 
  row <- row + 1
  apply(block_taxon_data[, 1:4], 1, block_plot, yrange=c(row-0.35, row+0.35),
        lowcol="cornsilk1", highcol="darkorange")
}

# dev.off()

