# R for plot Nucmer alignments
args <- commandArgs(trailingOnly=T)
blocktype <- args[1] # blocktying input
adjust_genebed <- args[2] # gene BED with coordinates adjusted
title <- args[3] # figure title
outpdf <- args[4] # PDF output file

cat(blocktype, "\n")
cat(adjust_genebed, "\n")
cat(title, "\n")
cat(outpdf, "\n")

#setwd("/homes/liu3zhen/scripts2/homotools/homograph_dev/test")
#blocktype <- "hgout/hgout.3.blocktyping"
#genebed <- "../0_gtfbed/GRMZM2G171650_T01.adjusted.bed"
#adjust_genebed <- "hgout/hgout.4.gene.bed"
#title <- "haplotype view"
#outpdf <- "out.pdf"

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

######################################################
### BED file
######################################################
is_bed_plot <- F
if (file.exists(adjust_genebed)) {
  file_info <- file.info(adjust_genebed)
  if (file_info$size > 0) {
    is_bed_plot <- T
    ymax <- ymax + 2
  }
}

highest_polymorphisms <- max(d[, 6])
if (highest_polymorphisms >= 1) {
  ymax <- ymax + 3
}

yheight <- 2 + (ymax - 2) * 0.15

# plot setting
# par(mar=c(2.2, 0.5, 1, 0.5), mgp=c(1.5,0.1,0))
pdf(outpdf, width=5, height=yheight)

par(mar=c(left_margin, 3.5, 2.2, 0.5), mgp=c(1, 0.65, 0))
plot(NULL, NULL, xlim=c(0, xmax), ylim=c(1, ymax),
     frame.plot=F, xlab="", ylab="",
     xaxt="n", yaxt="n")

### plot original gene structure
if (is_bed_plot) {
  max_adjust_height <- 0
  
  ### plot coordinate-adjusted gene structure
  newbed <- read.delim(adjust_genebed, header=F, comment.char="#", stringsAsFactors=F)
  is_underscore <- sum(newbed[, 4] == "UNDERSCORE")>0
  if (is_underscore) {
    underscore <- newbed[newbed[, 4] == "UNDERSCORE", ]
    newbed <- newbed[newbed[, 4] != "UNDERSCORE", ]
  }
  seqname <- unique(newbed[, 1])[1]
  newbed_min <- min(newbed[, 2])
  newbed_max <- max(newbed[, 3])
  ycenter_2 <- ymax
  lines(c(newbed_min, newbed_max), rep(ycenter_2, 2), col="gray50", lwd=0.6)
  if (nrow(newbed)>0) {
    for (i in 1:nrow(newbed)) {
      color <- newbed[i, 7]
      if (!isColor(color)) {
        color <- "grey"
      }
      height <- as.numeric(newbed[i, 5])
      adjust_height <- height*max(1, ymax/15)
      max_adjust_height <- max(max_adjust_height, adjust_height)
      rect(newbed[i, 2] + 1, ycenter_2 - adjust_height/2,
           newbed[i, 3], ycenter_2 + adjust_height/2,
           col=color, border=NA)
    }
  }
  text(newbed_min, ycenter_2, adj=2, labels=seqname, cex=0.6, xpd=T)
  
  ### underscore plotting
  if (is_underscore) {
    for (i in 1:nrow(underscore)) {
      lines(underscore[i, 2:3], rep(ycenter_2 - max_adjust_height, 2),
            lwd=1, col=underscore[i, 7], lend=1)
    }
  }
}

######################################################
# x-axis
######################################################
xaxis <- smartaxis(xmax)
axis(1, at=xaxis[[1]], labels=xaxis[[2]], tick=T, col="gray50", cex.axis=0.8)
mtext(paste0("position (", xaxis[[3]], ")"), side=1, line=1.8)

######################################################
# title
######################################################
mtext(title, side=3, line=1)

######################################################
### num of polymorphisms
######################################################
if (highest_polymorphisms >= 1) {
  ypos <- length(taxa) + 1.5
  log10max <- round(log10(highest_polymorphisms + 1), 0)
  pscale <- 2 / log10max
  hlines_val <- seq(0, log10max, by=1)
  lines_pos = ypos + hlines_val*pscale
  for (i in 1:length(lines_pos)) {
    lines(c(-xmax/100, xmax), rep(lines_pos[i], 2), lwd=0.5, lty=2, col="gray80")
    text(-xmax/100, lines_pos[i], adj=2, labels=hlines_val[i], las=1, cex=0.6)
  }
  #axis(2, at=c(ypos, ypos+log10max*pscale), labels=c(0, log10max),
  #     tick=T, col="gray50", cex.axis=0.6, las=1)
  text(-xmax/85, mean(lines_pos), labels="log10", adj=2, las=1, cex=0.7, xpd=T)
  #hlines_val <- seq(0, log10max, by=1)
  #abline(h=ypos + hlines_val*pscale, lwd=0.1, lty=2)
  # nums for plotting
  uniqseg <- !duplicated(paste(d[, 4], d[, 5], d[, 6]))
  if (sum(uniqseg)>0) {
    uniqseq_numpoly <- d[uniqseg, c(4,5,6)]
    for (i in 1:nrow(uniqseq_numpoly)) {
      plot_val <- log10(uniqseq_numpoly[i, 3] + 1)
      rect(uniqseq_numpoly[i, 1], ypos,
           uniqseq_numpoly[i, 2], ypos+plot_val*pscale,
           col="palevioletred", border=NA)
    }
  }
}
  
######################################################
# backgrounp of alignment regions (alignment segments)
######################################################
row <- 0
for (etaxon in taxa) {
  alnseg <- paste(d[, 1], d[, 2], d[, 3])
  taxon_data <- d[!duplicated(alnseg) & d[, 1] == etaxon, 1:3] 
  row <- row + 1
  apply(taxon_data[, 2:3], 1, block_plot, yrange=c(row-0.35, row+0.35), singlecol="lightblue")
  text(-xmax/80, row, adj=1, labels=etaxon, xpd=T, cex=0.6)
}

######################################################
# block and color-coded difference from the consensus
######################################################
row <- 0
for (etaxon in taxa) {
  block_taxon_data <- d[d[, 1] == etaxon, 4:7] 
  row <- row + 1
  apply(block_taxon_data[, 1:4], 1, block_plot, yrange=c(row-0.35, row+0.35),
        lowcol="cornsilk1", highcol="darkorange")
}

dev.off()

