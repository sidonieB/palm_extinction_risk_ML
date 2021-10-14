#
# command to run in parallel in a unix server
# input_tsv.txt contains two columns, one with the forest loss tile name, and the second with the mask tile name corresponding to the same geographic extent
# names in the input file have to include full paths
# parallel --colsep '\t' env TMPDIR=/tmp_dir/ R -f Bellot_resample_tiles.R --args Hansen/{1} mask/{2} :::: input_tsv.txt


library(raster)

# set different tempdir to ensure enough space
tempfile(tmpdir = "tmp_dir/")

setwd("forest_loss/")

args <- commandArgs(trailingOnly = TRUE)

hansen_tile = args[1]
mask_tile = args[2]

t <- raster(hansen_tile)
m <- raster(mask_tile)
t2 <- mask(t, m, maskvalue="0", updatevalue=NA)
t <- ""
t3 <- mask(t2, m, maskvalue="2", updatevalue=NA)
t2 <- ""
t4 <- raster::reclassify(t3, matrix(c(0,Inf,1), ncol=3))
t3 <- ""
t_lr4a <- aggregate(t4, fact=c(666.66666666667,666.66666666667), fun="mean", na.rm=T)
filename <- paste(gsub(".tif", "", hansen_tile), "_res20.tif", sep="")
writeRaster(t_lr4a, filename)

