# Steve Bachman, RBG Kew
# Script for Sidonie - part of Palms paper
# from list of all palms supplied by sidonie - use POWO API to extract TDWG distributions 
# and get count of TDWG per species

# run the function at the bottom first

library(httr)
library(dplyr)
library(tidyr)
library(jsonlite)
library(progress)
library(tidyverse)
library(plyr)


# import the raw file containing species names
input <- read.table("TDWG-extract_input.txt", stringsAsFactors=FALSE)

# drop the unwanted columns (if any)
reduced <- input[c(-1,-2)]

# deduplicate names
dedup <- (unique(reduced))

# now write down the results
results = write.table(
  dedup, # file name
  paste0("results.csv"), # file name
  row.names = FALSE, 
  na = "",
  sep = ","
)

# use  http://namematch.science.kew.org/ to get POWO ID - upload results
name_match_clean_results <- read.csv("name_match_clean_results.csv", stringsAsFactors=FALSE)

# test on single species
ID <- name_match_clean_results$id[800]
test_single = TDWG_from_POWO_ID(ID)

# now run on batch
# get IDs first
IDs = name_match_clean_results[, 3] #1220
#test_batch <- purrr::map_df(IDs, TDWG_from_POWO_ID) < not working? why?

# now run with lapply
tdwg_batch = lapply(IDs,TDWG_from_POWO_ID)
tdwg_batch_df = do.call(rbind, tdwg_batch)

# remove the NA species
tdwg_batch_NA = filter(tdwg_batch_df, tdwgCode == "NA")

# remove NAs from tdwg_batch_df and get count
tdwg_batch_native = filter(tdwg_batch_df, tdwgCode != "NA")

# get the count of natives
sp_count = tdwg_batch_native %>% dplyr::group_by(ID) %>% 
  dplyr::count(ID)

# subset on the NA IDs
NA_IDs = subset(tdwg_batch_NA, select = ID)
NA_IDs$n <- "NA"

# merge empty NAs with sp COUnt
sp_count_with_NA = rbind(sp_count_df, NA_IDs)

# join back to original table and save
final = merge(name_match_clean_results, sp_count_with_NA, by.x = "id", by.y = "ID")

# now write down the results
results = write.table(
  final, # file name
  paste0("Palms_TDWG3_count.txt"), # file name
  row.names = FALSE, 
  na = "",
  sep = "\t"
)

##########################################################
# functions - modified from version provided by baz walker
TDWG_from_POWO_ID = function(ID){
  
  # taxon <- "/taxon/urn:lsid:ipni.org:names:594092-1"
  search_url <- "http://powo.science.kew.org/api/2"
  stub <- "/taxon/"
  filters <- "?fields=distribution"
  
  query = paste0(search_url, stub, ID, filters)
  
  response <- GET(query)
  results <- fromJSON(content(response, as="text"))
  
  if (length(results$distribution$doubtful) >=1) {
    nat_dis = data.frame(
      establishment = "NA",
      featureId = "NA",
      tdwgCode = "NA",
      tdwgLevel = "NA",
      name = "NA",
      ID = ID)
  }
  
  else if (length(results$distribution$extinct) >=1) {
    nat_dis = data.frame(
      establishment = "NA",
      featureId = "NA",
      tdwgCode = "NA",
      tdwgLevel = "NA",
      name = "NA",
      ID = ID)
  }
  
  else if (length(results$distribution) == 0) {
    nat_dis = data.frame(
      establishment = "NA",
      featureId = "NA",
      tdwgCode = "NA",
      tdwgLevel = "NA",
      name = "NA",
      ID = ID)
  }
  
  #else if (!exists(results$distribution)) {
  #  nat_dis = data.frame(
  #    establishment = "NA",
  #    featureId = "NA",
  #    tdwgCode = "NA",
  #    tdwgLevel = "NA",
  #    name = "NA",
  #    ID = ID)
  #}
  
  else{
    nat_dis = results$distribution$natives
    nat_dis$ID = ID
  }
  print(ID)
  Sys.sleep(0.1)
  return(nat_dis)
  
}




