# R for plot Nucmer alignments

args <- commandArgs(trailingOnly=T)

bedfile <- args[1] # input BED file
stopifnot(file.exists(bedfile))

data_col <- args[2] # column of data
height <- args[3] # height in the 5th column of the output
lowcolor <- args[4] # color for low values
highcolor <- args[5] # color for high values
lowvalue <- args[6] # the low value for data_col for lowcolor
highvalue <- args[7] # the high value for data_col for highcol
gradient_layer <- args[8] # number for gradient layers
outfile <- args[9] # output file

if (length(args) != 9) {
	stop("Usage: Rscript script.R <bedfile> <data_col> <height> <lowcolor> <highcolor> <lowvalue> <highvalue> <gradient_layer> <outfile>")
}

# Print details for debugging
cat("<bedfile>, input BED file:", bedfile, "\n")
cat("<data_col>, data column number in BED:", data_col, "\n")
cat("<height>, height:", height, "\n")
cat("<lowcolor>, R compatible color for the lowest value:", lowcolor, "\n")
cat("<highcolor>,R compatible color for the highest value:", highcolor, "\n")
cat("<lowvalue>, the input lowest value overriding the lowest value in BED:", lowvalue, "\n")
cat("<highvalue>, the input lowest value overriding the highest value in BED:", highvalue, "\n")
cat("<gradient_layer> number for gradient layers:", gradient_layer, "n")
cat("<outfile>, output BED file:", outfile, "\n")

# convert characters to numbers
data_col <- as.numeric(data_col)
height <- as.numeric(height)
gradient_layer <- as.numeric(gradient_layer)
lowvalue <- as.numeric(lowvalue)
highvalue <- as.numeric(highvalue)

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
if (!isColor(lowcolor)) {
  cat(lowcolor, "is not valid. Plot using the default\n")
  lowcolor="gray90"
}

if (!isColor(highcolor)) {
  cat(highcolor, "is not valid. Plot using the default\n")
  highcolor="blue"
}

# gradient colors
colfunc <- colorRampPalette(c(lowcolor, highcolor))
color_num <- gradient_layer
allcols <- colfunc(color_num)

###########################################################
# BED
###########################################################
inbed <- read.delim(bedfile, header=F)
# clean up BED content
inbed[, 2] <- suppressWarnings(as.numeric(as.character(inbed[, 2])))
inbed[, 3] <- suppressWarnings(as.numeric(as.character(inbed[, 3])))
inbed[, data_col] <- suppressWarnings(as.numeric(as.character(inbed[, data_col])))
inbed <- inbed[!is.na(inbed[, 2]) & !is.na(inbed[, 3]) & !is.na(inbed[, data_col]), ]

if (is.null(lowvalue)) {
	lowvalue <- min(inbed[, data_col], na.rm=T)
}

if (is.null(highvalue)) {
	highvalue <- max(inbed[, data_col], na.rm=T)
}

value_points <- seq(lowvalue, highvalue, by=(highvalue - lowvalue)/(color_num -1))

gradient_cols <- array(1:nrow(inbed))
for (i in 1:nrow(inbed)) {
	gradient_cols[i] <- allcols[which.min(abs(value_points - inbed[i, data_col]))]
}
###########################################################
# output
###########################################################
outbed <- inbed[, 1:3]
outbed[, 4] <- "."
outbed[, 5] <- height

if (ncol(inbed)>5) {
	strand_col_values <- unique(inbed[, 6]) # strands
	if (sum( ! strand_col_values %in% c("+", "-")) > 0) { # either + or -
		outbed[, 6] <- inbed[, 6]
	} else {
		outbed[, 6] <- "+"
	}
} else {
	outbed[, 6] <- "+"
}

outbed[, 7] <- gradient_cols

write.table(outbed, outfile, row.names=F, quote=F, sep="\t", col.names=F)

