# R for plot Nucmer alignments

args <- commandArgs(trailingOnly=T)
datalist <- args[1] # nucmer show-coordi output
band_color <- args[2] # color to be used for plotting
identity_min <- args[3]
identity_max <- args[4]
outpdf <- args[5] # PDF output file

#setwd("/homes/liu3zhen/scripts2/homotools/case/homostack_case/")
#datalist <- "/homes/liu3zhen/scripts2/homotools/case/homostack_case/hsout/hsout.1.data.list"
#band_color <- "orange"

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
    
    if (startp[1] == endp[1]) {
      curve.x <- rep(startp[1], 2*npoint)
    } else {
      curve.x1 <- seq(startp[1], midp[1], by = (midp[1] - startp[1])/(npoint-1))
      curve.x2 <- seq(midp[1], endp[1], by = (endp[1] - midp[1])/(npoint-1))
      curve.x <- c(curve.x1, curve.x2)
    }
    
    if (startp[2] == endp[2]) {
      curve.y <- rep(startp[2], 2*npoint)
    } else {
      curve.y1 <- startp[2] - beizer_value * (startp[2] - midp[2])
      curve.y2 <- rev(endp[2] - beizer_value * (endp[2] - midp[2]))
      curve.y <- c(curve.y1, curve.y2)
    }
    list(x=curve.x, y=curve.y)
  }
  
  p1 <- transform_curve(c(topregion[1], topheight), c(bottomregion[1], bottomheight))
  p2 <- transform_curve(c(topregion[2], topheight), c(bottomregion[2], bottomheight))
  px <- c(p1$x, rev(p2$x))
  py <- c(p1$y, rev(p2$y))
  #lines(px, py)
  polygon(px, py, border=border, col=bandcol)
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
 
# check if the color is valid
if (!isColor(band_color)) {
  cat(band_color, "is not valid. Plot using the default\n")
  band_col="deepskyblue4"
}

# connection colors
high_identity_col <- band_color
low_identity_col <- "grey96"
colfunc <- colorRampPalette(c(low_identity_col, high_identity_col))
color_num <- 20
allcols <- colfunc(color_num)

###########################################################
# delta files to information
###########################################################
delta2info <- function(deltafile) {
  # input a delta filename and output
  # a data frame including: subj qry subj_len qry_len
  aln <- read.delim(deltafile, stringsAsFactor=F)
  stopifnot(nrow(aln)>=1)
  
  qlen <- unique(as.numeric(as.character(aln$qlen)))
  slen <- unique(as.numeric(as.character(aln$slen)))
  qry <- unique(aln$qry)
  subj <- unique(aln$subj)
  identity_high <- max(aln$identity)
  identity_low <- min(aln$identity) 

  data.frame(subj, slen, qry, qlen, identity_high, identity_low)	
}

# data list
dlist <- read.delim(datalist, stringsAsFactors=F)
num_aln <- nrow(dlist)

# xmax
xmax <- 0
for (i in 1:num_aln) { 
	aln_df <- delta2info(dlist[i, 2])
	xmax <- max(aln_df[1, c(2,4)], xmax)
}

# identity_min and _max
if (identity_min == "auto") {
  identity_min <- 100
  for (i in 1:num_aln) {
    aln_df <- delta2info(dlist[i, 2])
    identity_min <- min(c(aln_df[1, 6]), identity_min)
  }
}

if (identity_max == "auto") {
  identity_max <- 0
  for (i in 1:num_aln) {
    aln_df <- delta2info(dlist[i, 2])
  	identity_max <- max(c(aln_df[1, 5]), identity_max)
  }
}

identity_min <- as.numeric(identity_min)
identity_max <- as.numeric(identity_max)

# axis of x and y
xrange <- c(-xmax / 50, xmax + xmax / 50)
if ((identity_max - identity_min) < 5) {
	identity_min <- identity_max - 5
}
identity20 <- seq(identity_min, identity_max, by=(identity_max - identity_min)/(color_num -1))

yrange <- c(-0.15, num_aln + 0.5)

###########################################################
# function to add highlights
###########################################################
highlight_module <- function(bedfile, center, plot=TRUE) {  
  #chr start end label height strand color
  #10	1000	2488	gene	gray60	+	0.015
	highlight <- read.delim(bedfile, comment.char="#", stringsAsFactors=F, header=F)
	if (nrow(highlight)>0) {
	  highlight_upper_half_height <- max(highlight[, 5]) / 2
	  if (plot) {
	    highlight_names <- highlight[!duplicated(highlight[, 4]), 4]
	    highlight_colors <- highlight[!duplicated(highlight[, 4]), 7]
      for (i in 1:nrow(highlight)) {
	      color <- highlight[i, 7]
	      if (!isColor(color)) {
		      color <- "grey"
	  	  }
	      height <- highlight[i, 5]
          rect(highlight[i, 2] + 1, center - height / 2,
	           highlight[i, 3], center + height / 2,
	           col=color, border=color)
    }
	  }
	}
	highlight_upper_half_height # return to adjust band positions
}
###########################################################
# canvas
###########################################################
# pdf
#pdf(outpdf, width=4.5, height=1.25*num_aln+0.5)
#pdf(outpdf, width=4.5, height=2*sqrt(num_aln)+0.5)
pdf(outpdf, width=4.5, height=num_aln/1.5+2)

# plot canvas
par(mar=c(2.2, 0.5, 1, 0.5), mgp=c(1.5,0.1,0))
plot(NULL, NULL, type="n", xlim=xrange, ylim=yrange,
     xlab="", ylab="", bty="n", xaxt="n", yaxt="n")
  
# x-axis
xaxis <- smartaxis(xmax)
axis(1, at=xaxis[[1]], labels=xaxis[[2]], tick=F, col="gray50", cex.axis=0.8)
mtext(paste0("position (", xaxis[[3]], ")"), side=1, line=1)
  
# vertical lines
abline(v=xaxis[[1]], col="gray70", lwd=0.8)
abline(v=xaxis[[4]], col="gray90", lwd=0.5)
  
###########################################################
# plot sequences, highlights and connections
###########################################################
block_bottom <- 0
block_top <- 1
upper_half_height_default <- 0.005
upper_half_height <- upper_half_height_default
lower_half_height_default <- 0.005
lower_half_height <- lower_half_height_default
bar_band_dist <- 0.005

for (i in 1:num_aln) {
# region lines
  upper_ycenter <- block_top
  lower_ycenter <- block_bottom
  
  aln_file <- dlist[i, 2]
  subj_bed <- dlist[i, 4]
  qry_bed <- dlist[i, 6]
  
  aln <- read.delim(aln_file, stringsAsFactor=F)
  aln_df <- delta2info(aln_file)
  aln_info <- aln_df[1, ]
  subj <- aln_info$subj
  slen <- aln_info$slen
  qry <- aln_info$qry
  qlen <- aln_info$qlen
  
  # subj highlights
  if (i == 1) {
    rect(0, lower_ycenter - lower_half_height, slen, lower_ycenter + lower_half_height,
      col="gray70", border="gray70")
    lower_half_height <- highlight_module(bedfile=subj_bed, center=lower_ycenter)
  } else {
    lower_half_height <- highlight_module(bedfile=subj_bed, center=lower_ycenter, plot=F)
  }
  
  # qry and qry highlights
  rect(0, upper_ycenter - upper_half_height, qlen, upper_ycenter + upper_half_height,
       col="gray70", border="gray70")
  upper_half_height <- highlight_module(bedfile=qry_bed, center=upper_ycenter)
  
  # connections
  for (j in 1:nrow(aln)) {
    band_gradient_col <- allcols[which.min(abs(identity20 - aln[j, "identity"]))]
    bandconnect(topregion=aln[j, 3:4], bottomregion=aln[j, 1:2],
                topheight=upper_ycenter - upper_half_height - bar_band_dist,
                bottomheight=lower_ycenter + lower_half_height + bar_band_dist,
                border=NA, bandcol=band_gradient_col)
  }
  # text labels
  xlabels <- qry
  ylabels <- subj
  
  text(x= -xmax / 50, y=upper_ycenter - upper_half_height - bar_band_dist - 0.125,
       labels=xlabels, cex=0.8, xpd=T, pos=4)
  
  if (i == 1) {
    text(x= -xmax / 50, y=lower_ycenter - lower_half_height - bar_band_dist - 0.125,
         labels=ylabels, cex=0.8, xpd=T, pos=4)
  }
  #text(x= xmax / 2, y=lower_ycenter + lower_half_height,
	#   labels=transcript, cex=0.8, xpd=T, pos=3)
  block_bottom <- block_bottom + 1
  block_top <- block_top + 1
  lower_half_height <- lower_half_height_default
  upper_half_height <- upper_half_height_default
}

###########################################################
# legends
###########################################################

### identity legends
bar_ypos <- block_bottom + 0.3
barlen <- xmax / 4
barstep <- barlen / color_num
barheight <-  block_top / 20 / log(num_aln)
text(xmax - barlen/2, bar_ypos, labels = "identity", cex=0.8, pos=1, xpd=T)  ### plot LD name

# lowest identity
barlabels <- c(identity_min, identity_max)
barlabels.num <- floor(barlabels * color_num)
identity_min <- round(identity_min, 1)
text(xmax-barlen, bar_ypos+0.01, labels=identity_min, cex=0.7, pos=1)  ### plot low identity

# highest identity
identity_max <- round(identity_max, 1)
if (identity_max == 100) {
  identity_max <- round(identity_max, 0)
}
text(xmax, bar_ypos+0.01, labels=identity_max, cex=0.7, pos=1)  ### plot high identity

# color bars for identity legends
col_barpos <- xmax-barlen
for (i in 1:color_num) {
  rect(col_barpos, bar_ypos,
     col_barpos + barstep, bar_ypos + barheight,
     border=NA, col=allcols[i])
  col_barpos <- col_barpos + barstep
}

# close
dev.off() 
