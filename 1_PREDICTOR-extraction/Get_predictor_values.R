###################################################
# main script and workflow used to obtain all predictor values for each species, 
# based on species occurrence points
# other inputs are publicly available data layers
###################################################


###############################################################################
##obtain ConR variables
# EOO, AOO, number of populations, number of locations
# Later we will replace conR EOO and AOO by rCAT ones which seem more accurate
# and work with species spanning more than 180 degrees of longitude
# But we will keep the conR EOO for the species with only 2 occurrences,
# because rCAT does not calculate it
# Meanwhile we use a trick to keep species spanning more than 180 degrees in the table
# by calculating their EOO in two times, one for positive longitudes and one for negative longitudes
# then summing both
###############################################################################
library("ConR")

### Find species with occurrences spanning more than 180 degrees longitude
dat <- read.table("palm_occurences.txt", sep="\t", header=T)
str(dat)
dat_name <- as.character(dat$accepted_name)
name_list <- as.character(unique(dat$accepted_name))
str(name_list)
lon_list <- as.numeric(dat$lon)

name_lon_over180 <- c()
for (name in name_list) {
  position <- grep(name, dat_name, fixed = TRUE)
  lon <- lon_list[position]
  if (length(position) > 1) {
    lon_max <- max(lon)
    lon_min <- min(lon)
    lon_range <- lon_max - lon_min
    if (lon_range > 180) {
      name_lon_over180 <- c(name_lon_over180, name)
      print(name)
      print(lon_max)
      print(lon_min)
      print(lon_range)
    }
  }
}
str(name_lon_over180)

write.table(name_lon_over180, "palm_occurences_LONover180.txt", sep="\t", row.names = FALSE)

###rename species with occurrences spanning more than 180 degrees longitude so that EOO gets calculated for all minus and 
# for all plus coordinates separately. Not ideal but better than nothing.

dat2 <- dat
dat2$accepted_name <- as.character(dat2$accepted_name)

for (s in 1:length(dat2$accepted_name)) {
  spe <- dat2$accepted_name[s]
  lon <- as.numeric(dat2$lon[s])
  if (spe %in% name_lon_over180) {
    print(spe)
    if (lon > 0) {
      dat2$accepted_name[s] <- paste(dat2$accepted_name[s], "PLUS", sep="")
    } else {
      dat2$accepted_name[s] <- paste(dat2$accepted_name[s], "MINUS", sep="")
    }
  }
}

write.table(dat2, "palm_occurences_renamedForEOO.txt", sep="\t", row.names = F)

# generate conR variables for all species
# could also have started here if no species spanning >180 degrees longitude
dat <- read.table("palm_occurences_renamedForEOO.txt", sep="\t", header=T)
head(dat)
dat2 <- dat[,1:3]
dat2[,1] <- dat$lat
dat2[,2] <- dat$lon
dat2[,3] <- dat$accepted_name
colnames(dat2) <- c("lat", "lon", "accepted_name")

dat2$lon[47451] <- -179.4 # was rounded to -180 giving an error

dat2_iucn <- IUCN.eval(dat2)

dat3 <- dat2_iucn[1:6] # keep only relevant variables for us

# generate EOO for species with 2 occurrences (use conR EOO function and generate a new EOO column)
dat3_eoo <- EOO.computing(dat2, method.less.than3 = "arbitrary")

head(dat3_eoo)

dat3$EOO2 <- dat3_eoo$EOO

# some EOO remain uncomputed even with more than 1 occurrence when the occurrences are super close (probably rounded as equal)

# fuse values corresponding to + and - for each species with larger span than 180

dat4 <- dat3
dat4$taxa <- as.character(dat4$taxa)

for (x in 1:length(dat4$taxa)) {
  spe <- dat4$taxa[x]
  #print(spe)
  if(length(grep("PLUS",spe)) > 0) {
    print(spe)
    spe2 <- gsub("PLUS", "", spe)
    dat4$taxa[x] <- spe2
  } else if (length(grep("MINUS",spe)) > 0){
    print(spe)
    spe2 <- gsub("MINUS", "", spe)
    dat4$taxa[x] <- spe2
  }
}


dat4$Nbe_loc <- as.numeric(dat4$Nbe_loc)

dat4_sub <- split(dat4, dat4$taxa)

dat5 <- dat4

dat5$Nbe_loc <- as.numeric(dat5$Nbe_loc)

for (s in 1:length(dat4_sub)) {
  if (length(dat4_sub[[s]]$taxa) > 1) {
    # fuse rows by adding values
    fused <- dat4_sub[[s]][1,]
    for (x in 2:length(fused)) {
      fused[x] <- sum(dat4_sub[[s]][,x])
    }
    print(fused)
    # replace all copies by fused row
    copies <- grep(fused[1], dat5$taxa)
    print(copies)
    for (c in 1:length(copies)) {
      dat5[copies[c],] <- fused
    }
  }
}

# remove identical rows

dat6 <- unique(dat5)
length(dat6$taxa)


write.table(dat6, "conR_variables.txt", sep="\t", row.names = F)


#################################################
# Get EOO and AOO from rCAT
#################################################
# generate rCAT EOO and AOO, as well as a column with EOO from rCAT except for species with 2 points 
# where conR is kept because rCAT cannot generate it

library(rCAT)

#Import data (could also go on from above)

dat7 <- read.table("conR_variables.txt", sep = "\t", header=T)
#EOO: from conR with NA when 1 or 2 occurrences
#EOO2: from conR but with value when 2 occurrences  
#AOO: from conR

coords <- read.table("palm_occurences.txt", sep="\t", header=T, stringsAsFactors = F)
head(dat7)
head(coords)

results <- ConBatch(coords$accepted_name,coords$lat,coords$lon,2000,FALSE) 

dat8 <- merge(dat7,results,by="taxa")

head(dat8)

# create a new EOO column where rCAT EOO is replaced by conR EOO when 2 occs (because rCAT gives 0 for these)

dat8$EOO_rCAT_conR2 <- dat8$EOOkm2
for (x in 1:length(dat8$taxa)) {
  if (dat8$Nbe_unique_occ.[x] == 2) {
    dat8$EOO_rCAT_conR2[x] <- dat8$EOO2[x]
  }
}


write.table(dat8, "conR-rCAT_variables.txt", sep="\t", row.names = F)



##################################################################
##obtain other predictor variables and corresponding species values
# script is only for information
# input layers are not provided but can be obtained from sources
# indicated in the article. Output is provided to go on.
##################################################################

#load packgages
library(rgdal)
library(raster)
library(gridsample)
library(maps)
library(gdalUtils)

# details about variable sources are provided in the article or suppl. material
# variable / species value to get from each variable
# 19 bioclim variables (continuous) / mean for each - do it below
# human footprint index (continuous) / mean - do it below
# human population density (continuous) / mean - do it below
# ecoregions (polygons) / number of different ecoregions covered by the species - done with Ondo_extract_ecoregions_from_occ_mod.R
# TDWG count (from powo website) / number of different TDWG3 regions covered by the species - done with Bachman_tdwg_for_palms.R
# forest loss / proportion of forest loss in species range - see below for explanation


# Instead of extracting species values from high resolution layers using a buffer for the uncertainty of each occurrence point, 
# we will extract the values from low resolution (20km=10min) layers without buffering during the extraction, because the low resolution will already buffer this.
# This way is a bith rougher than extracting+buffering from high resolution layers but some layers are really too big at high resolution
# and coordinates uncertainty as well as layer data are rough too anyway.


### Prepare layers for the extraction of variables that can be processed here

##bioclimatic variables (already at 20km ~ 10 min)

# just need bio1 to use as a template below
bio1 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")


##human footprint index

# needs reprojection to WGS84
hfi <- raster("predictors/human_footprint/hfp_global_geo_grid/hf_v2geo/w001001.adf")
crs(hfi) <- "+proj=longlat +ellps=clrk66 +no_defs"
newproj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
pr2 <- projectRaster(hfi, crs=newproj)
writeRaster(pr2, "predictors/human_footprint/hfp_global_geo_grid/hf_v2geo/human_footprint_index_wgs84.tif")

# decrease resolution to 20km (by using bio1 as a template and resampling with bilinear method = same calculations than aggregate according to 
# https://gis.stackexchange.com/questions/255150/using-resample-vs-aggregate-extend-in-r-to-have-rasters-of-matching-resolutio/255155 but perfect overlap with bio: more convenient for later)
hfi <- raster("predictors/human_footprint/hfp_global_geo_grid/hf_v2geo/human_footprint_index_wgs84.tif")
hfi_lr <- resample(hfi, bio1)
writeRaster(hfi_lr, "predictors/human_footprint/hfp_global_geo_grid/hf_v2geo/human_footprint_index_wgs84_res20.tif")


##human population density

hpd <- raster("predictors/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_30_sec.tif")

# decrease resolution to 20km (by using bio1 as a template and resampling)
hpd_lr <- resample(hpd, bio1)
writeRaster(hpd_lr, "predictors/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_res20.tif")


### Preparation of the forest loss dataset

# The data are split into subset (tiles) of 10*10 degrees
# https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.7.html
# use wget -P . url to download each tile
# download the tile of the data (lossyear) and the corresponding tile of no-data zone (datamask)
# each tile is very heavy, so they cannot be readily all merged together
# tiles of lower resolution have to be first created as above, then these can hopefully be merged for species value extraction

## resample the tiles (will have to be performed in another computer, using Bellot_resample_tiles.R)

# for each tile of the forest loss, import it, import corresponding mask tile (*datamask* files, downloaded in the same link)
# mask the NA and water pixels in the forest tile using the mask tile
# reclassify the forest tile as 0 (no loss) and 1 (loss) + NAs, instead of having more values for loss
# resample using the aggregate function and a factor of 666.66666667, and the mean, so that we obtain larger pixels (same resolution as bio1)
# where each large pixel is the mean of the smaller pixels that constitute it in the original forest tile
# That way, because the values are 0 or 1, the mean will be equal to the proportion of the square that has some forest loss

# run Bellot_resample_tiles.R in a bigger computer in parallel mode

# obtain a low resolution tile *_res20.tif for each tile

## merge all low resolution tiles together (need first to correct extent so that perfectly adjacent)

all_tiles = list.files(pattern="*.tif")
all_tiles_c = vector(mode="list", length(all_tiles))

for (r in 1:length(all_tiles)) {
  a <- raster(all_tiles[r])
  a_e <- extent(extent(a)[1],extent(a)[2]-0.005, extent(a)[3]+0.005, extent(a)[4])
  a_c <- setExtent(a, a_e, keepres=FALSE, snap=FALSE)
  all_tiles_c[[r]] <- a_c
}

all_tiles_c$fun <- mean
all_tiles_c$na.rm <- TRUE
r_mos <- do.call(mosaic, all_tiles_c)

plot(r_mos)
writeRaster(r_mos, file="Forest_loss_all_tiles.tif", format="GTiff")


### extract species values from all the low resolution layers (bioclim, human footprit, human density and forest loss)
# import low resolution layers
bio1 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
bio2 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_2.tif")
bio3 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_3.tif")
bio4 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_4.tif")
bio5 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_5.tif")
bio6 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_6.tif")
bio7 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_7.tif")
bio8 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_8.tif")
bio9 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_9.tif")
bio10 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_10.tif")
bio11 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_11.tif")
bio12 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")
bio13 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_13.tif")
bio14 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_14.tif")
bio15 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_15.tif")
bio16 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_16.tif")
bio17 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_17.tif")
bio18 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_18.tif")
bio19 <- raster("predictors/wc2.1_10m_bio/wc2.1_10m_bio_19.tif")
hfi_lr <- raster("predictors/human_footprint/hfp_global_geo_grid/hf_v2geo/human_footprint_index_wgs84_res20.tif")
hpd_lr <- raster("predictors/gpw-v4-population-density-rev11_2015_30_sec_tif/gpw_v4_population_density_rev11_2015_res20.tif")
fol_lr <- raster("predictors/forest_loss/resampled/Forest_loss_all_tiles.tif")

# import species coordinates
spe_occs <- read.table("palm_occurences.txt", sep="\t", header=T)

# subset the species data and prepare a result table
spe_occs_sub <- split(spe_occs, spe_occs$accepted_name) 
results <- data.frame(names(spe_occs_sub))
results$HFI <- rep(NA, length(results[,1]))
results$HPD <- rep(NA, length(results[,1]))
results$FLO <- rep(NA, length(results[,1]))
for (x in 1:19) {
  results[,x+4] <- rep(NA, length(results[,1])) 
  names(results)[x+4] <- paste("bio",x,sep="")
}
for (x in 1:19) {
  results[,x+23] <- rep(NA, length(results[,1])) 
  names(results)[x+23] <- paste("bio",x, "_sd", sep="")
}
results$bio5max <- rep(NA, length(results[,1]))
results$bio6min <- rep(NA, length(results[,1]))
results$bio13max <- rep(NA, length(results[,1]))
results$bio14min <- rep(NA, length(results[,1]))
results$Tflex_ext <- rep(NA, length(results[,1]))
results$Pflex_ext <- rep(NA, length(results[,1]))
results$Tflex_mean <- rep(NA, length(results[,1]))
results$Pflex_mean <- rep(NA, length(results[,1]))
names(results)[1] <- "accepted_name"

# loop to extract the species values from each layer

for (s in 1:length(names(spe_occs_sub))) {
  
 
  spe_occs_spe <- spe_occs_sub[[s]]
  coords <- spe_occs_spe[,2:3]
  results$HFI[s] <- mean(extract(hfi_lr, coords), na.rm=T)
  results$HPD[s] <- mean(extract(hpd_lr, coords), na.rm=T)
  results$FLO[s] <- mean(extract(fol_lr, coords), na.rm=T)
  results$bio1[s] <- mean(extract(bio1, coords), na.rm=T)
  results$bio2[s] <- mean(extract(bio2, coords), na.rm=T)
  results$bio3[s] <- mean(extract(bio3, coords), na.rm=T)
  results$bio4[s] <- mean(extract(bio4, coords), na.rm=T)
  results$bio5[s] <- mean(extract(bio5, coords), na.rm=T)
  results$bio6[s] <- mean(extract(bio6, coords), na.rm=T)
  results$bio7[s] <- mean(extract(bio7, coords), na.rm=T)
  results$bio8[s] <- mean(extract(bio8, coords), na.rm=T)
  results$bio9[s] <- mean(extract(bio9, coords), na.rm=T)
  results$bio10[s] <- mean(extract(bio10, coords), na.rm=T)
  results$bio11[s] <- mean(extract(bio11, coords), na.rm=T)
  results$bio12[s] <- mean(extract(bio12, coords), na.rm=T)
  results$bio13[s] <- mean(extract(bio13, coords), na.rm=T)
  results$bio14[s] <- mean(extract(bio14, coords), na.rm=T)
  results$bio15[s] <- mean(extract(bio15, coords), na.rm=T)
  results$bio16[s] <- mean(extract(bio16, coords), na.rm=T)
  results$bio17[s] <- mean(extract(bio17, coords), na.rm=T)
  results$bio18[s] <- mean(extract(bio18, coords), na.rm=T)
  results$bio19[s] <- mean(extract(bio19, coords), na.rm=T)
  results$bio1_sd[s] <- sd(extract(bio1, coords), na.rm=T)
  results$bio2_sd[s] <- sd(extract(bio2, coords), na.rm=T)
  results$bio3_sd[s] <- sd(extract(bio3, coords), na.rm=T)
  results$bio4_sd[s] <- sd(extract(bio4, coords), na.rm=T)
  results$bio5_sd[s] <- sd(extract(bio5, coords), na.rm=T)
  results$bio6_sd[s] <- sd(extract(bio6, coords), na.rm=T)
  results$bio7_sd[s] <- sd(extract(bio7, coords), na.rm=T)
  results$bio8_sd[s] <- sd(extract(bio8, coords), na.rm=T)
  results$bio9_sd[s] <- sd(extract(bio9, coords), na.rm=T)
  results$bio10_sd[s] <- sd(extract(bio10, coords), na.rm=T)
  results$bio11_sd[s] <- sd(extract(bio11, coords), na.rm=T)
  results$bio12_sd[s] <- sd(extract(bio12, coords), na.rm=T)
  results$bio13_sd[s] <- sd(extract(bio13, coords), na.rm=T)
  results$bio14_sd[s] <- sd(extract(bio14, coords), na.rm=T)
  results$bio15_sd[s] <- sd(extract(bio15, coords), na.rm=T)
  results$bio16_sd[s] <- sd(extract(bio16, coords), na.rm=T)
  results$bio17_sd[s] <- sd(extract(bio17, coords), na.rm=T)
  results$bio18_sd[s] <- sd(extract(bio18, coords), na.rm=T)
  results$bio19_sd[s] <- sd(extract(bio19, coords), na.rm=T)
  results$bio5max[s] <- max(extract(bio5, coords), na.rm=T)
  results$bio6min[s] <- min(extract(bio6, coords), na.rm=T)
  results$bio13max[s] <- max(extract(bio13, coords), na.rm=T)
  results$bio14min[s] <- min(extract(bio14, coords), na.rm=T)
}

head(results)

# calculate Tflex and Pflex
# need first to shift all T values so that all >0, otherwise ratio won't mean the same for ranges spanning 0
# also, when the minimal (numerator) is 0, replace it by 0.00001 (arbitrary, chosen because 5 decimals in the values) 
# so that the ratio remains meaningful
# for Pflex, only concerns 78 species if use min and 5 if use mean (we will use mean in the article)
# for Tflex, only concerns the species that had the minimal temperature

minT <- abs(min(c(results$bio6min, results$bio5max, results$bio6, results$bio5)))

for (s in 1:length(results$accepted_name)) {
  if ((results$bio6min[s] + minT) != 0) {
    results$Tflex_ext[s] <- (results$bio6min[s] + minT) / (results$bio5max[s] + minT)
  } else {
    results$Tflex_ext[s] <- 0.00001 / (results$bio5max[s] + minT)  
  }
  
  if ( results$bio14min[s] != 0) {
    results$Pflex_ext[s] <- results$bio14min[s] / results$bio13max[s]
  } else {
    results$Pflex_ext[s] <- 0.00001 / results$bio13max[s]
  }

  if ( (results$bio6[s] + minT) != 0) {
    results$Tflex_mean[s] <- (results$bio6[s] + minT) / (results$bio5[s] + minT)
  } else {
    results$Tflex_mean[s] <- 0.00001 / (results$bio5[s] + minT)
  }
  
  if (results$bio14[s] != 0) {
    results$Pflex_mean[s] <- results$bio14[s] / results$bio13[s]
  } else {
    results$Pflex_mean[s] <- 0.00001 / results$bio13[s]
  }

}

# check that ratios are superior to 0, i.e. that shift worked
grep("TRUE", (results$Pflex_ext < 0))
grep("TRUE", (results$Pflex_mean < 0))
grep("TRUE", (results$Tflex_ext < 0))
grep("TRUE", (results$Tflex_mean < 0))

# check that ratios are inferior to 1 (should be fine but just in case of weird anomaly)
grep("TRUE", (results$Pflex_ext > 1))
grep("TRUE", (results$Pflex_mean > 1))
grep("TRUE", (results$Tflex_ext > 1))
grep("TRUE", (results$Tflex_mean > 1))

# check that ratios are not = 0
grep("TRUE", (results$Pflex_ext == 0))
grep("TRUE", (results$Pflex_mean == 0))
grep("TRUE", (results$Tflex_ext == 0))
grep("TRUE", (results$Tflex_mean == 0))


write.table(results, "HFI-HPD-FLO-BIOFLEX_PerSpecies.txt", sep="\t", row.names = F)


########################################################
# combine all predictors
# inputs (outputs from above steps) are provided
########################################################
### Import conR/rCAT values, ecoregions and TDWG counts and add them to result table

conR <- read.table("conR-rCAT_variables.txt", sep="\t", header = T)
H <- read.table("HFI-HPD-FLO-BIOFLEX_PerSpecies.txt", sep="\t", header = T)
ECO <- read.table("NumEcoregionsPerSpecies.txt", sep="\t", header = T)
TAD <- read.table("Palms_TDWG3_count.txt", sep="\t", header = T)

names(conR)[1] <- "accepted_name"
names(TAD)[1] <- "accepted_name"
names(TAD)[2] <- "TDWG3"

library(plyr)
all <- merge(conR,H, by="accepted_name")
all2 <- merge(all,ECO, by="accepted_name")
all3 <- join(all2, TAD, by="accepted_name", type="inner")

write.table(all3, "RANGE-HFI-HPD-FLO-BIO-ECO-TDWG_PerSpecies.txt", sep="\t", row.names = F)


### last tweak: replace EOO=0 by EOO=NA in the original column if numOcc <=2

all4 <- read.table("RANGE-HFI-HPD-FLO-BIO-ECO-TDWG_PerSpecies.txt", sep="\t", header = T)

for (x in 1:length(all4$EOO_rCAT_conR2)) {
  if (is.na(all4$EOO_rCAT_conR2[x])) {
    all4$EOO_rCAT_conR2[x] <- "NA"    
  } else if (as.character(all4$EOO_rCAT_conR2[x]) == "0" && all4$Nbe_unique_occ.[x] <= 2) {
    all4$EOO_rCAT_conR2[x] <- "NA"
  }
}

write.table(all4, "RANGE-HFI-HPD-FLO-BIO-ECO-TDWG_PerSpecies_EOONA.txt", sep="\t", row.names = F)


