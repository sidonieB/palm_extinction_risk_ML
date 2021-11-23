### combine extinction risk predictor data with published (IUCN or ThreatSearch) assessment data
cons <- read.table("palm_summary_merge29_simple.txt", sep="\t", header=T)
vars <- read.table("RANGE-HFI-HPD-FLO-BIO-ECO-TDWG_PerSpecies_EOONA.txt", sep="\t", header=T)

# make new conservation column with all detailed categories, keeping iucn new, or, if NE/DD, replacing it by TS new if the latter global exists
cons$IUCN_TS_detail <- as.character(cons$IUCN_category_new)
cons$IUCN_category_new <- as.character(cons$IUCN_category_new)
cons$TS_category_new <- as.character(cons$TS_category_new)

for (x in 1:length(cons$accepted_name)) {
  if (cons$IUCN_category_new[x] == "DD") {
    cons$IUCN_TS_detail[x] <- "NE"
  }
  if (cons$IUCN_TS_detail[x] == "NE") {
    if (cons$TS_global_new[x] > 0) {
      cons$IUCN_TS_detail[x] <- cons$TS_category_new[x]
    }
  }
  
}

cons$IUCN_TS_sum <- as.character(cons$IUCN_TS_detail)

cons$IUCN_TS_sum <- gsub(pattern="Not Threatened",replacement="LC",cons$IUCN_TS_sum)
cons$IUCN_TS_sum <- gsub(pattern="CR|EN|NT|VU|Threatened,Data Deficient|Near Threatened|Possibly Threatened|Threatened",replacement="nonLC",cons$IUCN_TS_sum)
cons$IUCN_TS_sum <- gsub(pattern="Data Deficient",replacement="NE",cons$IUCN_TS_sum)


levels(as.factor(cons$IUCN_TS_detail))
levels(as.factor(cons$IUCN_TS_sum))


# make similar columns but including iucn old if no iucn/TS new available, or TS old global if nothing else

cons$IUCN_TS_detail_InclOld <- as.character(cons$IUCN_TS_detail)
cons$IUCN_category_old <- as.character(cons$IUCN_category_old)
cons$TS_category_old <- as.character(cons$TS_category_old)

for (x in 1:length(cons$accepted_name)) {
  if (cons$IUCN_TS_sum[x] == "NE") {
    
    cons$IUCN_TS_detail_InclOld[x] <- cons$IUCN_category_old[x]
    if (cons$IUCN_TS_detail_InclOld[x] == "DD") {
      cons$IUCN_TS_detail_InclOld[x] <- "NE"
    }
    if (cons$IUCN_TS_detail_InclOld[x] == "NE") {
      if (cons$TS_global_new[x] > 0) {
        cons$IUCN_TS_detail_InclOld[x] <- cons$TS_category_new[x]
      } else if (cons$TS_global_old[x] > 0) {
        cons$IUCN_TS_detail_InclOld[x] <- cons$TS_category_old[x]
      }
    }
  }
}

cons$IUCN_TS_sum_InclOld <- as.character(cons$IUCN_TS_detail_InclOld)

cons$IUCN_TS_sum_InclOld <- gsub(pattern="Threatened,Not Threatened",replacement="nonLC",cons$IUCN_TS_sum_InclOld)
cons$IUCN_TS_sum_InclOld <- gsub(pattern="Not Threatened",replacement="LC",cons$IUCN_TS_sum_InclOld)
cons$IUCN_TS_sum_InclOld <- gsub(pattern="EW|LC,EN|CR,EN|CR|EN|NT|VU|Possibly Threatened,Threatened|Threatened,Possibly Threatened|Threatened,Data Deficient|Near Threatened|Possibly Threatened|Threatened",replacement="nonLC",cons$IUCN_TS_sum_InclOld)
cons$IUCN_TS_sum_InclOld <- gsub(pattern="Data Deficient",replacement="NE",cons$IUCN_TS_sum_InclOld)


levels(as.factor(cons$IUCN_TS_detail_InclOld))
levels(as.factor(cons$IUCN_TS_sum_InclOld))

# combine with predictors, keep species with no predictors (no occs), to have a table with everything, can subset later
# keep an indication of what species have predictors data or not (and therefore will be included in the ML) by creating a new preds column
# which will be NA for species with no predictors after merging with iucn data

library(plyr)

vars$preds <- rep("YES", length(vars$accepted_name))
cons_sub <- cons[,10:13] # keep only the new columns created in cons
cons_sub$accepted_name <- cons$accepted_name

all <- join(vars, cons_sub, by="accepted_name", type="full")
for (x in 1:length(all$preds)) {
  if (is.na(all$preds[x])) {
    all$preds[x] <- "NO"
  }
}

# keep only relevant predictors (remove many columns that were from the conR and rCAT analyses, as well as most bioclim vars)
drops <- c(names(all)[2:3], names(all)[7:10], names(all)[12:13], names(all)[18:20], names(all)[22:31], names(all)[33:61]) # columns to remove
all2 <- all[, !(names(all) %in% drops)]

# export data, separating what will go in the ML or not

all_notML <- split(all2, all2$preds)[[1]] # preds="NO": no occurrences, no predictors: not for ML
all_forML <- split(all2, all2$preds)[[2]] # preds="YES": occurrences and predictor values: good for ML

write.table(all2, "All_species_all_preds_and_cons.txt", sep="\t", row.names = F)
write.table(all_forML, "MLspecies_all_preds_and_cons.txt", sep="\t", row.names = F)
write.table(all_notML, "NoMLspecies_all_preds_and_cons.txt", sep="\t", row.names = F)

