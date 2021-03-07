bandconnect <- function(topregion, bottomregion, topheight, bottomheight,
                        border=NULL, bandcol="grey", method=c("sqrt", "log")) {
  
  stopifnot(sum(method %in% c("sqrt", "log")) > 0)
  ### module
  transform_curve <- function(startp, endp, method=c("sqrt", "log"), logbase=exp(1), npoint=1000) {
    ### computer transform_curve coordinates
    midp <- (startp + endp) / 2
    if (method == "log") {
      beizer_value <- log(base=logbase, x=1:(npoint)) / log(base=logbase, x=npoint)
    } else {
      beizer_value <- sqrt(1:npoint) / sqrt(npoint)  # sqrt as default
    }
    
    curve.x1 <- seq(startp[1], midp[1], by = (midp[1] - startp[1])/(npoint-1))
    curve.y1 <- startp[2] - beizer_value * (startp[2] - midp[2])
    
    curve.x2 <- seq(midp[1], endp[1], by = (endp[1] - midp[1])/(npoint-1))
    curve.y2 <- rev(endp[2] - beizer_value * (endp[2] - midp[2]))
    
    curve.x <- c(curve.x1, curve.x2)
    curve.y <- c(curve.y1, curve.y2)
    list(x=curve.x, y=curve.y)
  }
  
  p1 <- transform_curve(c(topregion[1], topheight), c(bottomregion[1], bottomheight), method=method)
  p2 <- transform_curve(c(topregion[2], topheight), c(bottomregion[2], bottomheight), method=method)
  px <- c(p1$x, rev(p2$x))
  py <- c(p1$y, rev(p2$y))
  lines(px, py)
  polygon(px, py, border=border, col=bandcol)
  #polygon(px, py, border=NA, col=bandcol)
}


region2vec <- function(inregion, inchr) {
# convert coordinates (e.g., "1:1000-2000") to position range (e.g., c(1000, 2000))
  .chr <- gsub("\\:.*", "", inregion)
  .range <- gsub(".*\\:", "", inregion)
  .start <- gsub("\\-.*", "", .range);
  .end <- gsub(".*\\-", "", .range)
  stopifnot(.chr==inchr)
  as.numeric(c(.start, .end))
}

matchplot <- function(match.file, top.bed, bottom.bed, top.name=NULL, bottom.name=NULL,
                     top.chr=NULL, top.start=NULL, top.end=NULL,
                     bottom.chr=NULL, bottom.start=NULL, bottom.end=NULL,
                     strand.col=c("cornflowerblue", "tan"), band.col="gray80",
                     hightlight.top.region=NULL, hightlight.top.region.pos=0.7, hightlight.top.region.col="seagreen",
                     hightlight.bottom.region=NULL, hightlight.bottom.region.pos=0.15, hightlight.bottom.region.col="seagreen",
                     highlight.top=NULL, highlight.top.col="brown1", highlight.top.text=F,
                     highlight.bottom=NULL, highlight.bottom.col="brown1", highlight.bottom.text=F,
                     highlight.band=NULL, highlight.band.col="orange",
                     gene.height=0.1, main="", hightlight.region.lwd=10,
                     top.xpos=0.1, bottom.xpos=0.1,
                     top.ypos=0.7, bottom.ypos=0.15,
                     top.scale=1, bottom.scale=1
                     ) {
  options(stringsAsFactors = F, scipen=999)
  par(mar=c(0,0,0,0))
  # read data
  match <- read.delim(match.file, comment.char="#")
  if (is.null(top.name)) { top.name <- colnames(match)[1] }
  if (is.null(bottom.name)) { bottom.name <- colnames(match)[2] }
  colnames(match) <- c(top.name, bottom.name)
  
  ### top information
  top <- read.delim(top.bed, header=F, comment.char="#")
  top <- top[, 1:6]
  colnames(top) <- c("topChr", "topStart", "topEnd", "topGene", "topOther", "topStrand")
  
  # limite genes in the match list and find the range
  top.tmp <- top[top[, 4] %in% match[, 1], ]
  
  if (!is.null(top.chr)) {
    top.tmp <- top.tmp[top.tmp[,1] == top.chr, ]
  } else {
    top.chr <- unique(top.tmp[, 1])
  }
  
  if (!is.null(top.start)) {
    top.tmp <- top.tmp[top.tmp[,2] >= top.start, ]
  } else {
    top.start <- min(top.tmp[,2] + 1)
  }
  
  if (!is.null(top.end)) {
    top.tmp <- top.tmp[top.tmp[,3] <= top.end, ]
  } else {
    top.end <-  max(top.tmp[,3])
  }
  
  ### find top genes in the range
  top <- top[top[,1]==top.chr & top[,2]>=top.start & top[,3]<=top.end, ]
  
  ### bottom information
  bottom <- read.delim(bottom.bed, header=F, comment.char="#")
  bottom <- bottom[, 1:6]
  colnames(bottom) <- c("bottomChr", "bottomStart", "bottomEnd", "bottomGene", "bottomOther", "bottomStrand")
  # limite genes in the match list and find the range
  bottom.tmp <- bottom[bottom[, 4] %in% match[, 2], ]
  
  if (!is.null(bottom.chr)) {
    bottom.tmp <- bottom.tmp[bottom.tmp[,1] == bottom.chr, ]
  } else {
    bottom.chr <- unique(bottom.tmp[, 1])
  }
  
  if (!is.null(bottom.start)) {
    bottom.tmp <- bottom.tmp[bottom.tmp[,2] >= bottom.start, ]
  } else {
    bottom.start <- min(bottom.tmp[,2] + 1)
  }
  
  if (!is.null(bottom.end)) {
    bottom.tmp <- bottom.tmp[bottom.tmp[,3] <= bottom.end, ]
  } else {
    bottom.end <-  max(bottom.tmp[,3])
  }

  ### find top genes in the range
  bottom <- bottom[bottom[,1]==bottom.chr & bottom[,2]>=bottom.start & bottom[,3]<=bottom.end, ]
  
  # qc data
  stopifnot(length(unique(top[,1])) == 1) # require data from a single chr
  stopifnot(length(unique(bottom[,1])) == 1) # # require data from a single chr
  
  # plot xrange
  top.min <- top.start
  top.max <- top.end
  top.range <- top.max - top.min
  
  bottom.min <- bottom.start
  bottom.max <- bottom.end
  bottom.range <- bottom.max - bottom.min

  # plot unit (bp if unscaled)
  top.unit <- (top.range * top.scale) / (1 - top.xpos)
  bottom.unit <- (bottom.range * bottom.scale) / (1 - bottom.xpos)
  plot.unit <- top.unit
  if (bottom.unit > top.unit) {
    plot.unit <- bottom.unit
  }
  
  # coordinate conversion modules
  top.conversion <- function(pos) { (pos - top.min) / plot.unit + top.xpos }
  bottom.conversion <- function(pos) { (pos - bottom.min) / plot.unit + bottom.xpos }
  
  # plot canvas and labels
  half.gene.height <- gene.height / 2
  plot(NULL, NULL, xlim=c(0,1), ylim=c(0,1), axes=F)
  text(0.5, 1, pos=1, labels=main, cex=1.2, col="palevioletred4")
  lines(top.conversion(c(top.start, top.end)), rep(top.ypos, 2), col="grey", lwd=2)
  top.label <- paste0(top.name,"  ", top.chr, ":", top.start, "-", top.end)
  text(0, top.ypos+2*half.gene.height, labels=top.label, pos=4, col="palevioletred4")
  lines(bottom.conversion(c(bottom.start, bottom.end)), rep(bottom.ypos, 2), col="grey", lwd=2)
  bottom.label <- paste0(bottom.name,"  ",bottom.chr, ":", bottom.start, "-", bottom.end)
  text(0, bottom.ypos-2.2*half.gene.height, labels=bottom.label, pos=4, col="palevioletred4")
  
  # top genes
  for (i in 1:nrow(top)) {
    if (top[i, 6] == "+") {
      rect(top.conversion(top[i,2]), top.ypos - half.gene.height, top.conversion(top[i,3]),
           top.ypos + half.gene.height, col=strand.col[1], border = strand.col[1])
    } else {
      rect(top.conversion(top[i,2]), top.ypos - half.gene.height, top.conversion(top[i,3]),
           top.ypos + half.gene.height, col=strand.col[2], border = strand.col[2])
    }
  }
  
  # bottom genes
  for (i in 1:nrow(bottom)) {
    if (bottom[i, 6] == "+") {
      rect(bottom.conversion(bottom[i,2]), bottom.ypos - half.gene.height, bottom.conversion(bottom[i,3]),
           bottom.ypos + half.gene.height, col=strand.col[1], border = strand.col[1])
    } else {
      rect(bottom.conversion(bottom[i,2]), bottom.ypos - half.gene.height, bottom.conversion(bottom[i,3]),
           bottom.ypos + half.gene.height, col=strand.col[2], border = strand.col[2])
    }
  }
  
  # bands
  link <- merge(match, top, by.x=top.name, by.y="topGene")
  link <- merge(link, bottom, by.x=bottom.name, by.y="bottomGene")
  link <- link[!duplicated(paste(link[, top.name], link[, bottom.name])), ]
  for (i in 1:nrow(link)) {
    topregion <- top.conversion(link[i, c("topStart", "topEnd")])
    topregion <- as.numeric(topregion)
    bottomregion <- bottom.conversion(link[i, c("bottomStart", "bottomEnd")])
    bottomregion <- as.numeric(bottomregion)
    bandconnect(topregion=topregion, bottomregion=bottomregion, method="sqrt",
                topheight=top.ypos-half.gene.height*1.1, bottomheight=bottom.ypos+half.gene.height*1.1,
                border=band.col, bandcol=band.col)
  }
  
  # highlight top region
  if (!is.null(hightlight.top.region)) {
    hightlight.top.coords <- region2vec(hightlight.top.region, top.chr)
    lines(top.conversion(hightlight.top.coords), rep(hightlight.top.region.pos, 2),
          col=hightlight.top.region.col, lwd=hightlight.region.lwd, lend=1)
  }
  
  # highlight top region
  if (!is.null(hightlight.bottom.region)) {
    hightlight.bottom.coords <- region2vec(hightlight.bottom.region, bottom.chr)
    lines(bottom.conversion(hightlight.bottom.coords), rep(hightlight.bottom.region.pos, 2),
          col=hightlight.bottom.region.col, lwd=hightlight.region.lwd, lend=1)
  }
  
  # highlight top
  if (!is.null(highlight.top)) {
    highlight.top.data <- read.delim(highlight.top, header = F)
    highlight.top.genes <- highlight.top.data[, 1]
    # top genes
    for (i in 1:nrow(top)) {
      htop <- top[top$topGene %in% highlight.top.genes, ]
      if (nrow(htop) > 0) {
        rect(top.conversion(htop[i,2]), top.ypos - half.gene.height, top.conversion(htop[i,3]),
             top.ypos + half.gene.height, col=highlight.top.col, border = highlight.top.col)
      }
    }
  }
  
  # highlight bottom
  if (!is.null(highlight.bottom)) {
    highlight.bottom.data <- read.delim(highlight.bottom, header = F)
    highlight.bottom.genes <- highlight.bottom.data[, 1]
    # bottom genes
    for (i in 1:nrow(bottom)) {
      hbottom <- bottom[bottom$bottomGene %in% highlight.bottom.genes, ]
      if (nrow(hbottom) > 0) {
        rect(bottom.conversion(hbottom[i,2]), bottom.ypos - half.gene.height, bottom.conversion(hbottom[i,3]),
             bottom.ypos + half.gene.height, col=highlight.bottom.col, border = highlight.bottom.col)
      }
    }
  }
  
  # highlight bands
  if (!is.null(highlight.band)) {
    highlight.band.data <- read.delim(highlight.band)
    colnames(highlight.band.data) <- c(top.name, bottom.name) 
    
    highlight.link <- merge(link, highlight.band.data, by=c(top.name, bottom.name))
    
    if (nrow(highlight.link) > 1) {
      for (i in 1:nrow(highlight.link)) {
        topregion <- top.conversion(highlight.link[i, c("topStart", "topEnd")])
        topregion <- as.numeric(topregion)
        bottomregion <- bottom.conversion(highlight.link[i, c("bottomStart", "bottomEnd")])
        bottomregion <- as.numeric(bottomregion)
        bandconnect(topregion=topregion, bottomregion=bottomregion, method="sqrt",
                    topheight=top.ypos-half.gene.height*1.1, bottomheight=bottom.ypos+half.gene.height*1.1,
                    border=highlight.band.col, bandcol=highlight.band.col)
      }
    }
  }
}
