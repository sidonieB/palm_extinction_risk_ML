library(plyr)
library(ggplot2)
library(ape)
library(hash)
library(ggtree)
library(ggstance)
library(ggpubr)
library(gridExtra)
library(caper)
library(funrar)
library(dplyr)
library(reshape2)
library(DMwR)
library(scales)
library(sf)
library(rgdal)
library(hrbrthemes)
library(viridis)
library(GGally)
library(ggrepel)
library(phytools)
library(phangorn)
library(stringr)
library(aplot)
library(ggplotify)

####################
# custom functions useful later
####################

# replace <NA> by "NA"
NA_rep <- function(table) {
  for ( x in 1:length(table[,1])) {
    for (y in 2:length(names(table))) {
      if (is.na(table[x,y])) {
        table[x,y] <- "NA"
      }
    }
  }
  return(table)
}

# make a function for the opposite (useful later)
NA_rep_opp <- function(table) {
  for ( x in 1:length(table[,1])) {
    for (y in 2:length(names(table))) {
      if (table[x,y] == "NA") {
        table[x,y] <- NA
      }
    }
  }
  return(table)
}

######################


#####################
# Proportions and numbers of threatened species
# Multiple biodiversity facets
# World and continents
######################
# We use the files compiled separately for higher precision and higher specificity data after having manually appended 
# respectively "_MOST" and "_LEAST" after their variable headers (only for the "Predicted*" and "subset_plot_*" variables
# since other variables, such as Continent or trait etc, do not differ between the two files)

df_P <- read.table("All_species_cons_and_traits_newly_combined_Prec.txt", sep="\t", header=T, stringsAsFactors = F)
df_S <- read.table("All_species_cons_and_traits_newly_combined_Sens.txt", sep="\t", header=T, stringsAsFactors = F)

# make a dataset with most and least conservative (i.e. best and worst case scenario) predictions by appending the *_LEAST variables of df_S to df_P

df_S_sb <- df_S[,c(1,40:42,44:48)]
df_PS <- left_join(df_P, df_S_sb, by="species")


# add FD/PD as classes
FD <- read.table("Functional_distinctiveness.txt", sep=",", stringsAsFactors = F, header=T)
ED <- read.table("Evolutionary_distinctiveness_Average.txt", sep="\t", stringsAsFactors = F, header=T)
ED2 <- as.data.frame(ED[,c(2,4)])
FD2 <- as.data.frame(FD[,c(2,1)])
names(ED2)[1] <- "species"

# a few tree names are synonyms under our taxonomic backbone, so we average their ED
# ideally one would drop tips of all synonyms and recalculate ED
# 93 species from our backbone are not in the tree: no ED for them and potential bias for the remaining species
# ok for our purpose, and anyway ED always subject to species concept
# output of the loop is manually exported to "ED_problems.txt"

FD2$ED <- rep(NA, length(FD2$species))

for (x in 1:length(FD2$species)) {
  if (length(grep(FD2$species[x], ED2$species)) > 0) {
    if (length(grep(FD2$species[x], ED2$species)) > 1) {
      av <- c()
      for (z in 1:length(grep(FD2$species[x], ED2$species))) {
        av[z] <- ED2$AverageED_noBurnin[grep(FD2$species[x], ED2$species)[z]]
      }
      FD2$ED[x] <- mean(av)
      print(FD2$species[x])
      print(av)
      print(FD2$ED[x])
    } else {
      FD2$ED[x] <- ED2$AverageED_noBurnin[which(ED2$species %in% FD2$species[x])[1]] 
    }
  } else {
    print(paste(FD2$species[x], " NO MATCH"))
  }
}

df_PS_EF <- left_join(df_PS, FD2, by="species")

# check what is higher vs lower than the average:
mean(df_PS_EF$FD_G)
mean(df_PS_EF$ED, na.rm = T)

# categorize FD and ED as high ("YES"), low ("NO") or NA, depending if higher or lower than average
df_PS_EF$FDcat <- rep(NA, length(df_PS_EF$species))
df_PS_EF$EDcat <- rep(NA, length(df_PS_EF$species))

# FD > mean(df_PS_EF$FD_G) is high; ED > mean(df_PS_EF$ED, na.rm = T) too
for (x in 1:length(df_PS_EF$species)) {
  if (!is.na(df_PS_EF$FD_G[x])) {
    #print("yes")
    if (df_PS_EF$FD_G[x] > mean(df_PS_EF$FD_G)) {
      df_PS_EF$FDcat[x] <- "YES"
    } else {
      df_PS_EF$FDcat[x] <- "NO"
    }
  } else {
    df_PS_EF$FDcat[x] <- "NA"
  }
  
  if (!is.na(df_PS_EF$ED[x])) {
    if (df_PS_EF$ED[x] > mean(df_PS_EF$ED, na.rm = T)) {
      df_PS_EF$EDcat[x] <- "YES"
    } else {
      df_PS_EF$EDcat[x] <- "NO"
    }
  } else {
    df_PS_EF$EDcat[x] <- "NA"
  }
  
}


# make a copy of the data (also to not have to change the df names all below)
df3_c <- df_PS_EF

# remove un-necessary columns for clarity
keep <- c(names(df3_c)[1], names(df3_c)[37:84])
df3_c_2 <- df3_c[,names(df3_c) %in% keep]

# add a few columns to subset the data more easily and clarify who has a old assessment and who hasn't etc

df3_c_2$REC <- rep("NA", length(df3_c_2$species))
df3_c_2$OLD <- rep("NA", length(df3_c_2$species))
df3_c_2$PRED <- rep("NA", length(df3_c_2$species))
df3_c_2$Subset2 <- rep("NA", length(df3_c_2$species))

# replace <NA> by "NA" for consistency across the tables
df3_c_3 <- NA_rep(df3_c_2)




# populate the new columns
for (x in 1:length(df3_c_3$species)) {
  if (df3_c_3$IUCN_TS_sum[x]=="NE") { df3_c_3$REC[x] <- "NO"} else {df3_c_3$REC[x] <- "YES"}
  if (df3_c_3$IUCN_TS_sum_InclOld[x]=="NE") { df3_c_3$OLD[x] <- "NO"} else {
    if (df3_c_3$IUCN_TS_sum[x]=="NE") {df3_c_3$OLD[x] <- "YES"} else {df3_c_3$OLD[x] <- "NO"}
  }
  if (df3_c_3$Predicted_class_LEAST[x]=="NA") {df3_c_3$PRED[x] <- "NO"} else {df3_c_3$PRED[x] <- "YES"}
  
  if (df3_c_3$REC[x]=="NO" && df3_c_3$OLD[x]=="NO" && df3_c_3$PRED[x]=="YES") {df3_c_3$Subset2[x] <- "NEWPRED_data"}
  if (df3_c_3$REC[x]=="NO" && df3_c_3$OLD[x]=="YES" && df3_c_3$PRED[x]=="YES") {df3_c_3$Subset2[x] <- "REPRED_data"}
  if (df3_c_3$REC[x]=="NO" && df3_c_3$OLD[x]=="YES" && df3_c_3$PRED[x]=="NO") {df3_c_3$Subset2[x] <- "NOPRED_data_old_noOcc"}
  if (df3_c_3$REC[x]=="YES" && df3_c_3$OLD[x]=="NO" && df3_c_3$PRED[x]=="NO") {
    if (df3_c_3$Subset[x]=="NO_ML") {df3_c_3$Subset2[x] <- "NOPRED_data_rec_noOcc"} else {df3_c_3$Subset2[x] <- "NOPRED_data_rec_TestTrainDypsis"}
  } 
  if (df3_c_3$Subset2[x]=="NA") {df3_c_3$Subset2[x] <- "NOPRED_data_noAss_noOcc"}
    }

# Add new columns where we will put the value of a species dep. on if using old assessments or not

df3_c_3$RLa <- rep("NA", length(df3_c_3$species))
df3_c_3$RLro <- rep("NA", length(df3_c_3$species))
df3_c_3$RLr <- rep("NA", length(df3_c_3$species))
df3_c_3$RLa_MLp_BS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLro_MLpr_BS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLr_MLpr_BS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLa_MLp_WS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLro_MLpr_WS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLr_MLpr_WS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLro_MLpr_EG_BS <- rep("NA", length(df3_c_3$species))
df3_c_3$RLro_MLpr_EG_WS <- rep("NA", length(df3_c_3$species))

# populate columns accordingly

for (x in 1:length(df3_c_3$species)) {
  df3_c_3$RLa[x] <- df3_c_3$IUCN_TS_sum_InclOld[x]
  if (df3_c_3$RLa[x]=="NE") {df3_c_3$RLa[x] <- "NA"}
  
  if (df3_c_3$IUCN_TS_sum[x] == "NE") {
    if (df3_c_3$PRED[x]=="NO") {
      df3_c_3$RLro[x] <- df3_c_3$IUCN_TS_sum_InclOld[x]
      if (df3_c_3$RLro[x]=="NE") {df3_c_3$RLro[x] <- "NA"}
    } else {df3_c_3$RLro[x] <- "NA"}
  }  else {df3_c_3$RLro[x] <- df3_c_3$IUCN_TS_sum[x]}
  
  df3_c_3$RLr[x] <- df3_c_3$IUCN_TS_sum[x]
  if (df3_c_3$RLr[x]=="NE") {df3_c_3$RLr[x] <- "NA"}
  
  if (df3_c_3$RLa[x] == "NA") {df3_c_3$RLa_MLp_BS[x] <- df3_c_3$Predicted_class_MOST[x]} else {df3_c_3$RLa_MLp_BS[x] <- df3_c_3$RLa[x]}
  if (df3_c_3$RLro[x] == "NA") {df3_c_3$RLro_MLpr_BS[x] <- df3_c_3$Predicted_class_MOST[x]} else {df3_c_3$RLro_MLpr_BS[x] <- df3_c_3$RLro[x]}
  if (df3_c_3$RLr[x] == "NA") {df3_c_3$RLr_MLpr_BS[x] <- df3_c_3$Predicted_class_MOST[x]} else {df3_c_3$RLr_MLpr_BS[x] <- df3_c_3$RLr[x]}
  if (df3_c_3$RLro_MLpr_BS[x] == "NA") {df3_c_3$RLro_MLpr_EG_BS[x] <- df3_c_3$Expert[x]} else {df3_c_3$RLro_MLpr_EG_BS[x] <- df3_c_3$RLro_MLpr_BS[x]}

  if (df3_c_3$RLa[x] == "NA") {df3_c_3$RLa_MLp_WS[x] <- df3_c_3$Predicted_class_LEAST[x]} else {df3_c_3$RLa_MLp_WS[x] <- df3_c_3$RLa[x]}
  if (df3_c_3$RLro[x] == "NA") {df3_c_3$RLro_MLpr_WS[x] <- df3_c_3$Predicted_class_LEAST[x]} else {df3_c_3$RLro_MLpr_WS[x] <- df3_c_3$RLro[x]}
  if (df3_c_3$RLr[x] == "NA") {df3_c_3$RLr_MLpr_WS[x] <- df3_c_3$Predicted_class_LEAST[x]} else {df3_c_3$RLr_MLpr_WS[x] <- df3_c_3$RLr[x]}
  if (df3_c_3$RLro_MLpr_WS[x] == "NA") {df3_c_3$RLro_MLpr_EG_WS[x] <- df3_c_3$Expert[x]} else {df3_c_3$RLro_MLpr_EG_WS[x] <- df3_c_3$RLro_MLpr_WS[x]}
  
}


# check

count(df3_c_3, RLa)
count(df3_c_3, RLro)
count(df3_c_3, RLr)
count(df3_c_3, RLa_MLp_BS)
count(df3_c_3, RLro_MLpr_BS)
count(df3_c_3, RLr_MLpr_BS)
count(df3_c_3, RLro_MLpr_EG_BS)
count(df3_c_3, RLa_MLp_WS)
count(df3_c_3, RLro_MLpr_WS)
count(df3_c_3, RLr_MLpr_WS)
count(df3_c_3, RLro_MLpr_EG_WS)

# export table in case want to use later

write.table(df3_c_3, "Best_case_Worst_case_table.txt", sep="\t", row.names = F)

# join this table with predictors table, to create a complete table for supplementary material

all_preds <- read.table("All_species_all_preds_and_cons.txt", sep="\t", header=T, stringsAsFactors = F)
names(all_preds)[1] <- "species"
all_preds2 <- merge(df3_c_3,all_preds,by="species")

write.table(all_preds2, "TableSx_Best_case_Worst_case_table_withPredInfo.txt", sep="\t", row.names = F)



# make a smaller table just for the lollipop proportion plot
keep <- c(names(df3_c_3)[1], names(df3_c_3)[14], names(df3_c_3)[33:37],names(df3_c_3)[48:49],names(df3_c_3)[54:64])
df3_c_4 <- df3_c_3[,names(df3_c_3) %in% keep]
head(df3_c_4)

# make subtable for all species regardless if used etc
keep <- c(names(df3_c_4)[1:2], names(df3_c_4)[10:20])
df3_c_4_sp <- df3_c_4[,names(df3_c_4) %in% keep]
head(df3_c_4_sp)

# replace NA by <NA>
df3_c_4_sp_na <- NA_rep_opp(df3_c_4_sp)
head(df3_c_4_sp_na)

# summarize - careful if change input, many things (column names) hard-coded
# make it as a function because need to do it 8 times
# the function gets both numbers and proportions

sum_props <- function(df_in) {

freqs <- table(df_in[,2:3], useNA = "no")
props <- (prop.table(freqs, margin=1))*100


glob <- as.data.frame("A_World")
glob[,2:3] <- (prop.table(table(df_in[,3], useNA = "no")))*100
names(glob) <- c("Continent", "LC", "nonLC")
glob[,4:5] <- table(df_in[,3], useNA = "no")
names(glob) <- c("Continent", paste(names(glob)[2],names(df_in)[3], "PER", sep="__"), paste(names(glob)[3],names(df_in)[3], "PER", sep="__"), paste(names(glob)[2],names(df_in)[3], "NUM", sep="__"), paste(names(glob)[3],names(df_in)[3], "NUM", sep="__"))


props2 <- as.data.frame.matrix(props)
names(props2) <- c(paste(names(props2)[1],names(df_in)[3], "PER", sep="__"), paste(names(props2)[2],names(df_in)[3], "PER", sep="__"))
freqs2 <- as.data.frame.matrix(freqs)
names(freqs2) <- c(paste(names(freqs2)[1],names(df_in)[3], "NUM", sep="__"), paste(names(freqs2)[2],names(df_in)[3], "NUM", sep="__"))
props2 <- cbind(props2, freqs2)
props2$Continent <- rownames(props2)

props3 <- rbind(glob, props2)
rownames(props3)<-c(1:4)

for (v in 4:length(names(df_in))) {
df <- as.data.frame.matrix((prop.table(table(df_in[,c(2,v)], useNA = "no"), margin=1))*100)
names(df) <- c(paste(names(df)[1],names(df_in)[v], "PER", sep="__"), paste(names(df)[2],names(df_in)[v], "PER", sep="__"))
df2 <- as.data.frame.matrix(table(df_in[,c(2,v)], useNA = "no"))
names(df2) <- c(paste(names(df2)[1],names(df_in)[v], "NUM", sep="__"), paste(names(df2)[2],names(df_in)[v], "NUM", sep="__"))
df <- cbind(df, df2)
df$Continent <- rownames(df)

glob <- as.data.frame("A_World")
glob[,2:3] <- (prop.table(table(df_in[,v], useNA = "no")))*100
names(glob) <- c("Continent", "LC", "nonLC")
glob[,4:5] <- table(df_in[,v], useNA = "no")
names(glob) <- c("Continent", paste(names(glob)[2],names(df_in)[v], "PER", sep="__"), paste(names(glob)[3],names(df_in)[v], "PER", sep="__"), paste(names(glob)[2],names(df_in)[v], "NUM", sep="__"), paste(names(glob)[3],names(df_in)[v], "NUM", sep="__"))

df <- rbind(glob, df)
rownames(df)<-c(1:4)

props3 <- merge(props3,df,by="Continent")
}
 return(props3)
}



props_sp <- sum_props(df3_c_4_sp_na)



# subset table for categories, do same prop tables and rbind them all for the plot

# subset
df3_c_4_u <- subset(df3_c_4, df3_c_4$Human_use_all == "YES")
df3_c_4_uF <- subset(df3_c_4, df3_c_4$Human_use_food == "YES")
df3_c_4_uMa <- subset(df3_c_4, df3_c_4$Human_use_materials == "YES")
df3_c_4_uMe <- subset(df3_c_4, df3_c_4$Human_use_medicine == "YES")
df3_c_4_uC <- subset(df3_c_4, df3_c_4$Human_use_cultural == "YES")
df3_c_4_f <- subset(df3_c_4, df3_c_4$FDcat == "YES")
df3_c_4_e <- subset(df3_c_4, df3_c_4$EDcat == "YES")
df3_c_4_i <- subset(df3_c_4, df3_c_4$Human_use_all == "YES" | df3_c_4$FDcat == "YES" | df3_c_4$EDcat == "YES")

# make subtables with less columns
keep <- c(names(df3_c_4_u)[1:2], names(df3_c_4_u)[10:20])
df3_c_4_u_sp <- df3_c_4_u[,names(df3_c_4_u) %in% keep]
head(df3_c_4_u_sp)

keep <- c(names(df3_c_4_uF)[1:2], names(df3_c_4_uF)[10:20])
df3_c_4_uF_sp <- df3_c_4_uF[,names(df3_c_4_uF) %in% keep]
head(df3_c_4_uF_sp)

keep <- c(names(df3_c_4_uMa)[1:2], names(df3_c_4_uMa)[10:20])
df3_c_4_uMa_sp <- df3_c_4_uMa[,names(df3_c_4_uMa) %in% keep]
head(df3_c_4_uMa_sp)

keep <- c(names(df3_c_4_uMe)[1:2], names(df3_c_4_uMe)[10:20])
df3_c_4_uMe_sp <- df3_c_4_uMe[,names(df3_c_4_uMe) %in% keep]
head(df3_c_4_uMe_sp)

keep <- c(names(df3_c_4_uC)[1:2], names(df3_c_4_uC)[10:20])
df3_c_4_uC_sp <- df3_c_4_uC[,names(df3_c_4_uC) %in% keep]
head(df3_c_4_uC_sp)

keep <- c(names(df3_c_4_f)[1:2], names(df3_c_4_f)[10:20])
df3_c_4_f_sp <- df3_c_4_f[,names(df3_c_4_f) %in% keep]
head(df3_c_4_f_sp)

keep <- c(names(df3_c_4_e)[1:2], names(df3_c_4_e)[10:20])
df3_c_4_e_sp <- df3_c_4_e[,names(df3_c_4_e) %in% keep]
head(df3_c_4_e_sp)

keep <- c(names(df3_c_4_i)[1:2], names(df3_c_4_i)[10:20])
df3_c_4_i_sp <- df3_c_4_i[,names(df3_c_4_i) %in% keep]
head(df3_c_4_i_sp)


# replace NA by <NA>
df3_c_4_u_sp_na <- NA_rep_opp(df3_c_4_u_sp)
df3_c_4_uF_sp_na <- NA_rep_opp(df3_c_4_uF_sp)
df3_c_4_uMa_sp_na <- NA_rep_opp(df3_c_4_uMa_sp)
df3_c_4_uMe_sp_na <- NA_rep_opp(df3_c_4_uMe_sp)
df3_c_4_uC_sp_na <- NA_rep_opp(df3_c_4_uC_sp)
df3_c_4_f_sp_na <- NA_rep_opp(df3_c_4_f_sp)
df3_c_4_e_sp_na <- NA_rep_opp(df3_c_4_e_sp)
df3_c_4_i_sp_na <- NA_rep_opp(df3_c_4_i_sp)


# summarize each table separately

props_u <- sum_props(df3_c_4_u_sp_na)
props_uF <- sum_props(df3_c_4_uF_sp_na)
props_uMa <- sum_props(df3_c_4_uMa_sp_na)
props_uMe <- sum_props(df3_c_4_uMe_sp_na)
props_uC <- sum_props(df3_c_4_uC_sp_na)
props_f <- sum_props(df3_c_4_f_sp_na)
props_e <- sum_props(df3_c_4_e_sp_na)
props_i <- sum_props(df3_c_4_i_sp_na)

# Add a Category column

props_sp$Category <- rep("All_species", length(props_sp$Continent))
props_u$Category <- rep("Used_species", length(props_u$Continent))
props_uF$Category <- rep("Used_speciesF", length(props_uF$Continent))
props_uMa$Category <- rep("Used_speciesMa", length(props_uMa$Continent))
props_uMe$Category <- rep("Used_speciesMe", length(props_uMe$Continent))
props_uC$Category <- rep("Used_speciesC", length(props_uC$Continent))
props_f$Category <- rep("FD_species", length(props_f$Continent))
props_e$Category <- rep("ED_species", length(props_e$Continent))
props_i$Category <- rep("Interesting_species", length(props_i$Continent))

# combine all tables
props_all <- rbind(props_sp, props_u, props_uF, props_uMa, props_uMe, props_uC, props_f, props_e, props_i)

# remove LC proportions
props_all_nLC <- props_all[,c(1, grep("nonLC",names(props_all)), length(names(props_all)))]

# export proportion table

write.table(props_all, "TableSx_Numbers_and_Proportions2.txt", sep="\t")



# make the plot (first remove LC props, and also it may be that RLa vs RLa_ML* + RLr vs RLr_MLpr* are enough, ie no RLro)

 props_all_nLC2 <- props_all_nLC
 props_all_nLC2$Category <- as.factor(props_all_nLC2$Category)

p <- ggplot(props_all_nLC2) +
  geom_hline(yintercept = 50, linetype="dotted", color = "grey", size=0.25)+
  geom_segment( aes(x=Category, xend=Category, y=nonLC__RLr_MLpr_BS, yend=nonLC__RLr_MLpr_WS), color="#6baed6") +
  geom_point( aes(x=Category, y=nonLC__RLa), colour="#000000", size=5, shape=17 )+
  geom_point( aes(x=Category, y=nonLC__RLr), colour="#636363", size=5, shape=17 )+
  #geom_point( aes(x=Category, y=nonLC__RLa_MLp_BS), colour="#56B4E9", size=3) +      #clear blue
  #geom_point( aes(x=Category, y=nonLC__RLa_MLp_WS), colour="#F0E442", size=3 ) +    #yellow
  geom_point( aes(x=Category, y=nonLC__RLr_MLpr_BS), colour="#0072B2", size=3 ) +   #dark blue
  geom_point( aes(x=Category, y=nonLC__RLr_MLpr_WS), colour="#D55E00", size=3 ) +   #orange
  coord_flip()+
  scale_x_discrete(limits = rev(levels(props_all_nLC2$Category)))+
  facet_grid(~ Continent) +
  theme_classic() +
  xlab("") +
  ylab("Percentage of species assessed/predicted to be threatened")

pdf("Figure1A_Proportions.pdf", 10, 4) 
p
dev.off()

setEPS()
postscript("Figure1A_Proportions.eps", height =4, width = 10)
p
dev.off()


# #############################


########################
# Plots of species without recent assessment 
# [not used, but run in case need dataset later]
#########################

# keep only stuff not recently assessed
df3_c_3_s <- subset(df3_c_3, df3_c_3$Subset2 == "NEWPRED_data" | df3_c_3$Subset2 == "NOPRED_data_noAss_noOcc" | df3_c_3$Subset2 == "NOPRED_data_old_noOcc" | df3_c_3$Subset2 == "REPRED_data")

# subset table with only relevant info
keep <- c(names(df3_c_3_s)[1], names(df3_c_3_s)[14], names(df3_c_3_s)[5], names(df3_c_3_s)[38])
df3_c_3_sub <- df3_c_3_s[,names(df3_c_3_s) %in% keep]
head(df3_c_3_sub)
df3_c_3_sub$Predicted_class_MOST <- as.factor(df3_c_3_sub$Predicted_class_MOST)
df3_c_3_sub$Predicted_class_LEAST <- as.factor(df3_c_3_sub$Predicted_class_LEAST)

# rename status to fit with paper terminology
df3_c_3_sub$Predicted_class_MOST <- mapvalues(df3_c_3_sub$Predicted_class_MOST, from = levels(df3_c_3_sub$Predicted_class_MOST), to = c("Non threatened", "Unknown", "Threatened"))
df3_c_3_sub$Predicted_class_LEAST <- mapvalues(df3_c_3_sub$Predicted_class_LEAST, from = levels(df3_c_3_sub$Predicted_class_LEAST), to = c("Non threatened", "Unknown", "Threatened"))

# rename also columns for nicer names in plot
names(df3_c_3_sub)[c(2,4)] <- c("Best", "Worst")

# add a column "Both" which will classify species as T/NT only is same class across best and worse scenario

df3_c_3_sub$Both <- rep("Unknown", length(df3_c_3_sub$species))

for (x in 1:length(df3_c_3_sub$species)){
  ifelse(df3_c_3_sub$Best[x] == df3_c_3_sub$Worst[x], df3_c_3_sub$Both[x] <- as.character(df3_c_3_sub$Best[x]), df3_c_3_sub$Both[x] <- "Unknown")
}

df3_c_3_sub$Both <- as.factor(df3_c_3_sub$Both)
df3_c_3_sub$Best <- factor(df3_c_3_sub$Best, levels = c("Non threatened", "Threatened", "Unknown"))
df3_c_3_sub$Worst <- factor(df3_c_3_sub$Worst, levels = c("Non threatened", "Threatened", "Unknown"))

# melt table
df3_c_3_sub_m <- melt(df3_c_3_sub, id.vars = c("species", "Continent"))
head(df3_c_3_sub_m)
tail(df3_c_3_sub_m)

# subset by continent
df3_c_3_sub_m_Af <- subset(df3_c_3_sub_m, df3_c_3_sub_m$Continent == "AfricaWAsia")
df3_c_3_sub_m_Am <- subset(df3_c_3_sub_m, df3_c_3_sub_m$Continent == "Americas")
df3_c_3_sub_m_As <- subset(df3_c_3_sub_m, df3_c_3_sub_m$Continent == "EAsiaPacific")

df3_c_3_sub_m$value <- as.factor(df3_c_3_sub_m$value)
df3_c_3_sub_m_Af$value <- as.factor(df3_c_3_sub_m_Af$value)
df3_c_3_sub_m_Am$value <- as.factor(df3_c_3_sub_m_Am$value)
df3_c_3_sub_m_As$value <- as.factor(df3_c_3_sub_m_As$value)

df3_c_3_sub_m$variable <- factor(df3_c_3_sub_m$variable, levels = c("Best", "Worst", "Both"))
df3_c_3_sub_m_Af$variable <- factor(df3_c_3_sub_m_Af$variable, levels = c("Best", "Worst", "Both"))
df3_c_3_sub_m_Am$variable <- factor(df3_c_3_sub_m_Am$variable, levels = c("Best", "Worst", "Both"))
df3_c_3_sub_m_As$variable <- factor(df3_c_3_sub_m_As$variable, levels = c("Best", "Worst", "Both"))


NumbersW <- df3_c_3_sub_m %>% group_by(variable, value) %>% tally()
NumbersAm <- df3_c_3_sub_m_Am %>% group_by(variable, value) %>% tally()
NumbersAf <- df3_c_3_sub_m_Af %>% group_by(variable, value) %>% tally()
NumbersAs <- df3_c_3_sub_m_As %>% group_by(variable, value) %>% tally()

NumbersW$region <- rep("World", 9)
NumbersAm$region <- rep("America", 9)
NumbersAf$region <- rep("Africa", 9)
NumbersAs$region <- rep("Asia", 9)

Numbers_all <- rbind(NumbersW, NumbersAm, NumbersAf, NumbersAs)

write.table(Numbers_all, "TableSx_Numbers_of_species_to_RedList_inclAcross.txt", sep="\t")

# make plots
p1 <- ggplot(data=df3_c_3_sub_m, aes(x=variable, fill=variable)) + 
  geom_bar() + facet_grid(~ value) + 
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"))+
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25)) +
  labs(y = "Number of species", x = "")+
  theme(legend.position="none")+
  ylim(0,1020)

pAf <- ggplot(data=df3_c_3_sub_m_Af, aes(x=variable, fill=variable)) + 
  geom_bar() + facet_grid(~ value) + 
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"))+
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25)) +
  labs(y = "Number of species", x = "")+
  theme(legend.position="none")+
  ylim(0,650)

pAm <- ggplot(data=df3_c_3_sub_m_Am, aes(x=variable, fill=variable)) + 
  geom_bar() + facet_grid(~ value) + 
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"))+
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25)) +
  labs(y = "Number of species", x = "")+
  theme(legend.position="none")+
  ylim(0,650)

pAs <- ggplot(data=df3_c_3_sub_m_As, aes(x=variable, fill=variable)) + 
  geom_bar() + facet_grid(~ value) + 
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73"))+
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25)) +
  labs(y = "Number of species", x = "")+
  theme(legend.position="none")+
  ylim(0,650)


pdf("FigureX_MLpr_WAmAfAs2.pdf", 5, 12) 
grid.arrange(p1, pAm, pAf, pAs, nrow=4)
dev.off()

setEPS()
postscript("FigureX_MLpr_WAmAfAs.eps", height =12, width = 5)
grid.arrange(p1, pAm, pAf, pAs, nrow=4)
dev.off()

##########################



##########################
# Threatened species without recent assessment
# same as above but new trait class as ED/FD/EDFD/EDU/FDU/EDFDU and use only T species
###########################

# go back from original table with all the info
# create a new trait class
df3_c_3_s$AllTraits <- rep(NA, length(df3_c_3_s$species))

for (x in 1:length(df3_c_3_s$species)){
 
  v <- c(df3_c_3_s$Human_use_all[x], df3_c_3_s$FDcat[x], df3_c_3_s$EDcat[x])
  v2 <- c("NO", "NO", "NO")
  if (v[1] == "YES") {v2[1] <- "U"}
  if (v[2] == "YES") {v2[2] <- "FD"}
  if (v[3] == "YES") {v2[3] <- "ED"}
  res <- v2[1]
  for (z in 2:length(v2)){
    res <- paste(res, v2[z], sep="_")
  }
  df3_c_3_s$AllTraits[x] <- res
    
}

df3_c_3_s$AllTraits <- as.factor(df3_c_3_s$AllTraits)

df3_c_3_s$AllTraits <- mapvalues(df3_c_3_s$AllTraits, from = levels(df3_c_3_s$AllTraits), to = c("FD_ED", "FD", "ED", "LP", "U_FD_ED", "U_FD", "U_ED", "U"))

  
# separate T and NT species (just do plot with species predicted as T in one or the other scenario for main figure)  

# subset table with only relevant info
keep <- c(names(df3_c_3_s)[1], names(df3_c_3_s)[14], names(df3_c_3_s)[5], names(df3_c_3_s)[38], names(df3_c_3_s)[65])
df3_c_3_sub2 <- df3_c_3_s[,names(df3_c_3_s) %in% keep]
head(df3_c_3_sub2)
df3_c_3_sub2$Predicted_class_MOST <- as.factor(df3_c_3_sub2$Predicted_class_MOST)
df3_c_3_sub2$Predicted_class_LEAST <- as.factor(df3_c_3_sub2$Predicted_class_LEAST)

# rename status to fit with paper terminology
df3_c_3_sub2$Predicted_class_MOST <- mapvalues(df3_c_3_sub2$Predicted_class_MOST, from = levels(df3_c_3_sub2$Predicted_class_MOST), to = c("Non threatened", "Unknown", "Threatened"))
df3_c_3_sub2$Predicted_class_LEAST <- mapvalues(df3_c_3_sub2$Predicted_class_LEAST, from = levels(df3_c_3_sub2$Predicted_class_LEAST), to = c("Non threatened", "Unknown", "Threatened"))

# rename also columns for nicer names in plot
names(df3_c_3_sub2)[c(2,4)] <- c("Best", "Worst")

# reorder factor levels

df3_c_3_sub2$Best <- factor(df3_c_3_sub2$Best, levels = c("Non threatened", "Threatened", "Unknown"))
df3_c_3_sub2$Worst <- factor(df3_c_3_sub2$Worst, levels = c("Non threatened", "Threatened", "Unknown"))

# melt table
df3_c_3_sub2_m <- melt(df3_c_3_sub2, id.vars = c("species", "Continent", "AllTraits"))
head(df3_c_3_sub2_m)
tail(df3_c_3_sub2_m)

# reorder trait levels for nicer plot
df3_c_3_sub2_m$AllTraits <- factor(df3_c_3_sub2_m$AllTraits, levels = c("U_FD_ED", "U_FD", "U_ED", "FD_ED", "U", "FD", "ED", "LP"))

# subset by T/NT/U
df3_c_3_sub2_m_T <- subset(df3_c_3_sub2_m, df3_c_3_sub2_m$value == "Threatened") # only using this for now
df3_c_3_sub2_m_NT <- subset(df3_c_3_sub2_m, df3_c_3_sub2_m$value == "Non threatened")
df3_c_3_sub2_m_U <- subset(df3_c_3_sub2_m, df3_c_3_sub2_m$value == "Unknown")

# subset by continent
df3_c_3_sub2_m_TAf <- subset(df3_c_3_sub2_m_T, df3_c_3_sub2_m_T$Continent == "AfricaWAsia")
df3_c_3_sub2_m_TAm <- subset(df3_c_3_sub2_m_T, df3_c_3_sub2_m_T$Continent == "Americas")
df3_c_3_sub2_m_TAs <- subset(df3_c_3_sub2_m_T, df3_c_3_sub2_m_T$Continent == "EAsiaPacific")

df3_c_3_sub2_m_T$value <- as.factor(df3_c_3_sub2_m_T$value)
df3_c_3_sub2_m_TAf$value <- as.factor(df3_c_3_sub2_m_TAf$value)
df3_c_3_sub2_m_TAm$value <- as.factor(df3_c_3_sub2_m_TAm$value)
df3_c_3_sub2_m_TAs$value <- as.factor(df3_c_3_sub2_m_TAs$value)



# create tables without continent info

df3_c_3_sub2_m_Tw <- as.data.frame(table(df3_c_3_sub2_m_T[3:5]))
df3_c_3_sub2_m_TAf2 <- as.data.frame(table(df3_c_3_sub2_m_TAf[3:5]))
df3_c_3_sub2_m_TAm2 <- as.data.frame(table(df3_c_3_sub2_m_TAm[3:5]))
df3_c_3_sub2_m_TAs2 <- as.data.frame(table(df3_c_3_sub2_m_TAs[3:5]))


# rearrange data to be able to plot numbers
df3_c_3_sub2_m_Tw_sorted <- arrange(df3_c_3_sub2_m_Tw, "variable", "AllTraits")
df3_c_3_sub2_m_TAf2_sorted <- arrange(df3_c_3_sub2_m_TAf2, "variable", "AllTraits")
df3_c_3_sub2_m_TAm2_sorted <- arrange(df3_c_3_sub2_m_TAm2, "variable", "AllTraits")
df3_c_3_sub2_m_TAs2_sorted <- arrange(df3_c_3_sub2_m_TAs2, "variable", "AllTraits")

df3_c_3_sub2_m_Tw_sorted_cumsum <- ddply(df3_c_3_sub2_m_Tw_sorted, "variable", transform, label_ypos=cumsum(Freq))
df3_c_3_sub2_m_TAf2_sorted_cumsum <- ddply(df3_c_3_sub2_m_TAf2_sorted, "variable", transform, label_ypos=cumsum(Freq))
df3_c_3_sub2_m_TAm2_sorted_cumsum <- ddply(df3_c_3_sub2_m_TAm2_sorted, "variable", transform, label_ypos=cumsum(Freq))
df3_c_3_sub2_m_TAs2_sorted_cumsum <- ddply(df3_c_3_sub2_m_TAs2_sorted, "variable", transform, label_ypos=cumsum(Freq))

df3_c_3_sub2_m_Tw_sorted_cumsum$AllTraits <- factor(df3_c_3_sub2_m_Tw_sorted_cumsum$AllTraits, levels = c("LP", "ED", "FD", "U", "FD_ED", "U_ED", "U_FD", "U_FD_ED"))
df3_c_3_sub2_m_TAf2_sorted_cumsum$AllTraits <- factor(df3_c_3_sub2_m_TAf2_sorted_cumsum$AllTraits, levels = c("LP", "ED", "FD", "U", "FD_ED", "U_ED", "U_FD", "U_FD_ED"))
df3_c_3_sub2_m_TAm2_sorted_cumsum$AllTraits <- factor(df3_c_3_sub2_m_TAm2_sorted_cumsum$AllTraits, levels = c("LP", "ED", "FD", "U", "FD_ED", "U_ED", "U_FD", "U_FD_ED"))
df3_c_3_sub2_m_TAs2_sorted_cumsum$AllTraits <- factor(df3_c_3_sub2_m_TAs2_sorted_cumsum$AllTraits, levels = c("LP", "ED", "FD", "U", "FD_ED", "U_ED", "U_FD", "U_FD_ED"))

df3_c_3_sub2_m_Tw_sorted_cumsum$variable <- factor(df3_c_3_sub2_m_Tw_sorted_cumsum$variable, levels = c("Worst", "Best"))
df3_c_3_sub2_m_TAf2_sorted_cumsum$variable <- factor(df3_c_3_sub2_m_TAf2_sorted_cumsum$variable, levels = c("Worst", "Best"))
df3_c_3_sub2_m_TAm2_sorted_cumsum$variable <- factor(df3_c_3_sub2_m_TAm2_sorted_cumsum$variable, levels = c("Worst", "Best"))
df3_c_3_sub2_m_TAs2_sorted_cumsum$variable <- factor(df3_c_3_sub2_m_TAs2_sorted_cumsum$variable, levels = c("Worst", "Best"))


# plot
p1 <- ggplot(df3_c_3_sub2_m_Tw_sorted_cumsum, aes(fill=AllTraits, y=Freq, x=variable)) + 
  geom_bar(stat="identity")+
  geom_text(aes(y =label_ypos, label=paste(Freq, label_ypos, sep="\n")), vjust=0.5, color="grey", size=2.5)+
  scale_fill_viridis(discrete = T, direction = -1) +
  labs(fill="AllTraits")+
  #ggtitle("World") +
  xlab("Scenario") +
  ylab("Number of species")+
  theme_minimal()+
  theme(legend.position="none", legend.box = "horizontal")+ 
  theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))+
  ylim(0,1020)+ coord_flip()



pAf <- ggplot(df3_c_3_sub2_m_TAf2_sorted_cumsum, aes(fill=AllTraits, y=Freq, x=variable)) + 
  geom_bar(stat="identity")+
  geom_text(aes(y =label_ypos, label=paste(Freq, label_ypos, sep="\n")), vjust=0.5, color="grey", size=2.5)+
  scale_fill_viridis(discrete = T, direction = -1) +
  labs(fill="AllTraits")+
  #ggtitle("World") +
  xlab("Scenario") +
  ylab("Number of species")+
  theme_minimal()+
  theme(legend.position="none", legend.box = "horizontal")+ 
  theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))+
  ylim(0,40)+ coord_flip()

pAm <- ggplot(df3_c_3_sub2_m_TAm2_sorted_cumsum, aes(fill=AllTraits, y=Freq, x=variable)) + 
  geom_bar(stat="identity")+
  geom_text(aes(y =label_ypos, label=paste(Freq, label_ypos, sep="\n")), vjust=0.5, color="grey", size=2.5)+
  scale_fill_viridis(discrete = T, direction = -1) +
  labs(fill="AllTraits")+
  #ggtitle("World") +
  xlab("Scenario") +
  ylab("Number of species")+
  theme_minimal()+
  theme(legend.position="none", legend.box = "horizontal")+ 
  theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))+
  ylim(0,375)+ coord_flip()

pAs <- ggplot(df3_c_3_sub2_m_TAs2_sorted_cumsum, aes(fill=AllTraits, y=Freq, x=variable)) + 
  geom_bar(stat="identity")+
  geom_text(aes(y =label_ypos, label=paste(Freq, label_ypos, sep="\n")), vjust=0.5, color="grey", size=2.5)+
  scale_fill_viridis(discrete = T, direction = -1) +
  labs(fill="AllTraits")+
  #ggtitle("World") +
  xlab("Scenario") +
  ylab("Number of species")+
  theme_minimal()+
  theme(legend.position="none", legend.box = "horizontal")+ 
  theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))+
  ylim(0,625)+ coord_flip()

pL <- ggplot(df3_c_3_sub2_m_TAs2_sorted_cumsum, aes(fill=AllTraits, y=Freq, x=variable)) + 
  geom_bar(stat="identity")+
  geom_text(aes(y =label_ypos, label=paste(Freq, label_ypos, sep="\n")), vjust=0.5, color="grey", size=2.5)+
  scale_fill_viridis(discrete = T, direction = -1) +
  labs(fill="AllTraits")+
  #ggtitle("World") +
  xlab("Scenario") +
  ylab("Number of species")+
  theme_minimal()+
  theme(legend.position="bottom", legend.box = "horizontal")+ 
  theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10))+
  ylim(0,625)+ coord_flip()


pdf("FigureS2B_MLpr_WAmAfAs_Tonly.pdf", 6, 10) 
grid.arrange(p1, pAm, pAf, pAs, pL, nrow=5)
dev.off()

setEPS()
postscript("FigureS2B_MLpr_WAmAfAs_Tonly.eps", height =10, width = 6)
grid.arrange(p1, pAm, pAf, pAs, pL, nrow=5)
dev.off()

################################


###############################
# Priority regions maps, tables and bump charts
# Also the gaps maps (not used)
# This is where df12 is generated
###############################

# Map proportion of threatened species for interesting species in each TDWG3 region
# do it RLr only, RLr_MLpr_B and RLr_MLpr_W
# do the same but with detail of uses, i.e score region = %T of species used for food, etc
# do the same just for ED and just for FD and just for all uses
# --> 7 maps*3scenarios
# Plot all 7 proportions for 15? most populated regions, one plot per scenario

# Recode uses as 1 if yes, 0 if NO or NA (so that can later add 1s per region to calculate region total number of used species for that use)
df <- read.table("Best_case_Worst_case_table.txt", header = T, sep = "\t")
df2 <- NA_rep(df)
df2$HUF <- mapvalues(as.factor(df2$Human_use_food), from = levels(as.factor(df2$Human_use_food)), to = c(0,0,1))
df2$HUMa <- mapvalues(as.factor(df2$Human_use_materials), from = levels(as.factor(df2$Human_use_materials)), to = c(0,0,1))
df2$HUMe <- mapvalues(as.factor(df2$Human_use_medicine), from = levels(as.factor(df2$Human_use_medicine)), to = c(0,0,1))
df2$HUC <- mapvalues(as.factor(df2$Human_use_cultural), from = levels(as.factor(df2$Human_use_cultural)), to = c(0,0,1))
df2$HUA <- mapvalues(as.factor(df2$Human_use_all), from = levels(as.factor(df2$Human_use_all)), to = c(0,0,1))
df2$FD2 <- mapvalues(as.factor(df2$FDcat), from = levels(as.factor(df2$FDcat)), to = c(0,1))
df2$ED2 <- mapvalues(as.factor(df2$EDcat), from = levels(as.factor(df2$EDcat)), to = c(0,0,1))
df2$COM <- rep(NA, length(df2$species))

# change factors to nums for HUF etc or below does not work
df2$HUF <- as.numeric(as.character(df2$HUF))
df2$HUMa <- as.numeric(as.character(df2$HUMa))
df2$HUMe <- as.numeric(as.character(df2$HUMe))
df2$HUC <- as.numeric(as.character(df2$HUC))
df2$HUA <- as.numeric(as.character(df2$HUA))
df2$FD2 <- as.numeric(as.character(df2$FD2))
df2$ED2 <- as.numeric(as.character(df2$ED2))
df2$COM <- as.numeric(as.character(df2$COM))

# populate the COM (combined) column that will have 1 if sum(HUA,FD2,ED2) > 0, thereby flagging "interesting" species
for (x in 1:length(df2$species)) {
  ifelse(sum(df2$FD2[x], df2$ED2[x], df2$HUA[x], na.rm = T)>0, df2$COM[x] <- 1, df2$COM[x] <- 0)
}

 
# import data to get tdwg3 region codes

ipni <- read.table("../TDWG3/species_ipni-id.txt", sep="\t", stringsAsFactors = F, header=T)
powo <- read.table("../TDWG3/checklist_species_palms.txt", sep="\t", stringsAsFactors = F, header=T)
dist <- read.table("../TDWG3/dist.txt", sep="\t", stringsAsFactors = F, header=T)

# had to manually replace a few apostrophs by nothing in the checklist otherwise imported file was wrong 
# (because R considered apostrophes as quotes)

head(powo)
length(powo$plant_name_id) #should be around 7416

df6 <- left_join(df2, ipni, by="species")
df6$db_id <- rep(NA, length(df6$species))
for (x in 1:length(df6$species)) {
  #df6$db_id[x] <- powo$db_id[which(powo$ipni_id == df6$ipni_id[x])]
  if (length(which(paste("_",powo$ipni_id,"_",sep="") %in% paste("_",df6$ipni_id[x],"_",sep = "")))>0) {  # just to be sure that the match will be strict
    df6$db_id[x] <- powo$db_id[which(paste("_",powo$ipni_id,"_",sep="") %in% paste("_",df6$ipni_id[x],"_",sep = ""))]
  } else {
    df6$db_id[x] <- NA
    print(df6$species[x])
  }
}

# one species (Dypsis declivium) with no db_id: will have NA in its distribution below
# could manually correct to Madagascar (checked location using ipni; record on ipni was accessed using ipni-id but is flagged as suppressed)
# but better remove it as we also remove lines with location doubtful etc

# get distribution data from dist using join

df7 <- left_join(df6,dist, by="db_id")

summary(as.factor(df7$introduced))
summary(as.factor(df7$extinct))
summary(as.factor(df7$location_doubtful))
summary(as.factor(df7$location_to_be_added))
summary(as.factor(df7$plant_not_present))

# remove Dypsis declivium

df8 <- df7[-grep("Dypsis declivium", df7$species),]

# filter all these out except "introduced" (because what matters for our purpose is that the species is in the region
# not that it is native to it)

df9 <- df8[which(df8$extinct == 0),]
df10 <- df9[which(df9$location_doubtful == 0),]
df11 <- df10[which(df10$location_to_be_added == 0),]
df12 <- df11[which(df11$plant_not_present == 0),]

summary(as.factor(df12$introduced))
summary(as.factor(df12$extinct))
summary(as.factor(df12$location_doubtful))
summary(as.factor(df12$location_to_be_added))
summary(as.factor(df12$plant_not_present))

# get the tdwg3 shape file
tdwg3 <- st_read("../TDWG3/level3/level3.shp")


# build dictionaries with threatened and all species used - do each once for each use/interest, RLr, best case and worst case

# NA and 0 allow later to distinguish regions that do not have palms from those that do not have concerned palms

# dict with the sum of species with uses
regions_h_a_HUF <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
# dict with the sum of threatened species with uses (best and worst case scenario and RL recent only)
regions_h_tB_HUF <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_HUF <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_HUF <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_HUMa <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_HUMa <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_HUMa <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_HUMa <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_HUMe <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_HUMe <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_HUMe <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_HUMe <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_HUC <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_HUC <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_HUC <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_HUC <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_HUA <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_HUA <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_HUA <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_HUA <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_FD2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_FD2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_FD2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_FD2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_ED2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_ED2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_ED2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_ED2 <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))

regions_h_a_COM <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_tB_COM <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_tW_COM <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))
regions_h_r_COM <- hash(keys=tdwg3$LEVEL3_COD, values=rep(0, length(tdwg3$LEVEL3_COD)))


# fill the dictionaries

for (x in 1:length(df12$species)) {
  
  if (is.na(regions_h_a_HUF[[df12$area_code_l3[x]]])) {
    regions_h_a_HUF[[df12$area_code_l3[x]]] <- 0 + df12$HUF[x]
  } else {
    regions_h_a_HUF[[df12$area_code_l3[x]]] <- regions_h_a_HUF[[df12$area_code_l3[x]]] + df12$HUF[x]  
  }
  
  if (is.na(regions_h_a_HUMa[[df12$area_code_l3[x]]])) {
    regions_h_a_HUMa[[df12$area_code_l3[x]]] <- 0 + df12$HUMa[x]
  } else {
    regions_h_a_HUMa[[df12$area_code_l3[x]]] <- regions_h_a_HUMa[[df12$area_code_l3[x]]] + df12$HUMa[x]  
  }
  
  if (is.na(regions_h_a_HUMe[[df12$area_code_l3[x]]])) {
    regions_h_a_HUMe[[df12$area_code_l3[x]]] <- 0 + df12$HUMe[x]
  } else {
    regions_h_a_HUMe[[df12$area_code_l3[x]]] <- regions_h_a_HUMe[[df12$area_code_l3[x]]] + df12$HUMe[x]  
  }
  
  if (is.na(regions_h_a_HUC[[df12$area_code_l3[x]]])) {
    regions_h_a_HUC[[df12$area_code_l3[x]]] <- 0 + df12$HUC[x]
  } else {
    regions_h_a_HUC[[df12$area_code_l3[x]]] <- regions_h_a_HUC[[df12$area_code_l3[x]]] + df12$HUC[x]  
  }

  if (is.na(regions_h_a_HUA[[df12$area_code_l3[x]]])) {
    regions_h_a_HUA[[df12$area_code_l3[x]]] <- 0 + df12$HUA[x]
  } else {
    regions_h_a_HUA[[df12$area_code_l3[x]]] <- regions_h_a_HUA[[df12$area_code_l3[x]]] + df12$HUA[x]  
  }

  if (is.na(regions_h_a_FD2[[df12$area_code_l3[x]]])) {
    regions_h_a_FD2[[df12$area_code_l3[x]]] <- 0 + df12$FD2[x]
  } else {
    regions_h_a_FD2[[df12$area_code_l3[x]]] <- regions_h_a_FD2[[df12$area_code_l3[x]]] + df12$FD2[x]  
  }

  if (is.na(regions_h_a_ED2[[df12$area_code_l3[x]]])) {
    regions_h_a_ED2[[df12$area_code_l3[x]]] <- 0 + df12$ED2[x]
  } else {
    regions_h_a_ED2[[df12$area_code_l3[x]]] <- regions_h_a_ED2[[df12$area_code_l3[x]]] + df12$ED2[x]  
  }

  if (is.na(regions_h_a_COM[[df12$area_code_l3[x]]])) {
    regions_h_a_COM[[df12$area_code_l3[x]]] <- 0 + df12$COM[x]
  } else {
    regions_h_a_COM[[df12$area_code_l3[x]]] <- regions_h_a_COM[[df12$area_code_l3[x]]] + df12$COM[x]  
  }
  
  
  
  if (df12$RLr_MLpr_BS[x] == "nonLC") {
    
    if (is.na(regions_h_tB_HUF[[df12$area_code_l3[x]]])) {
      regions_h_tB_HUF[[df12$area_code_l3[x]]] <- 0 + df12$HUF[x]
    } else {
      regions_h_tB_HUF[[df12$area_code_l3[x]]] <- regions_h_tB_HUF[[df12$area_code_l3[x]]] + df12$HUF[x]  
    } 
    
    if (is.na(regions_h_tB_HUMa[[df12$area_code_l3[x]]])) {
      regions_h_tB_HUMa[[df12$area_code_l3[x]]] <- 0 + df12$HUMa[x]
    } else {
      regions_h_tB_HUMa[[df12$area_code_l3[x]]] <- regions_h_tB_HUMa[[df12$area_code_l3[x]]] + df12$HUMa[x]  
    }
    
    if (is.na(regions_h_tB_HUMe[[df12$area_code_l3[x]]])) {
      regions_h_tB_HUMe[[df12$area_code_l3[x]]] <- 0 + df12$HUMe[x]
    } else {
      regions_h_tB_HUMe[[df12$area_code_l3[x]]] <- regions_h_tB_HUMe[[df12$area_code_l3[x]]] + df12$HUMe[x]  
    }
    
    if (is.na(regions_h_tB_HUC[[df12$area_code_l3[x]]])) {
      regions_h_tB_HUC[[df12$area_code_l3[x]]] <- 0 + df12$HUC[x]
    } else {
      regions_h_tB_HUC[[df12$area_code_l3[x]]] <- regions_h_tB_HUC[[df12$area_code_l3[x]]] + df12$HUC[x]  
    }

    if (is.na(regions_h_tB_HUA[[df12$area_code_l3[x]]])) {
      regions_h_tB_HUA[[df12$area_code_l3[x]]] <- 0 + df12$HUA[x]
    } else {
      regions_h_tB_HUA[[df12$area_code_l3[x]]] <- regions_h_tB_HUA[[df12$area_code_l3[x]]] + df12$HUA[x]  
    }    
    
    if (is.na(regions_h_tB_FD2[[df12$area_code_l3[x]]])) {
      regions_h_tB_FD2[[df12$area_code_l3[x]]] <- 0 + df12$FD2[x]
    } else {
      regions_h_tB_FD2[[df12$area_code_l3[x]]] <- regions_h_tB_FD2[[df12$area_code_l3[x]]] + df12$FD2[x]  
    }    
    
    if (is.na(regions_h_tB_ED2[[df12$area_code_l3[x]]])) {
      regions_h_tB_ED2[[df12$area_code_l3[x]]] <- 0 + df12$ED2[x]
    } else {
      regions_h_tB_ED2[[df12$area_code_l3[x]]] <- regions_h_tB_ED2[[df12$area_code_l3[x]]] + df12$ED2[x]  
    }
    
    if (is.na(regions_h_tB_COM[[df12$area_code_l3[x]]])) {
      regions_h_tB_COM[[df12$area_code_l3[x]]] <- 0 + df12$COM[x]
    } else {
      regions_h_tB_COM[[df12$area_code_l3[x]]] <- regions_h_tB_COM[[df12$area_code_l3[x]]] + df12$COM[x]  
    }
    
  }
  
  
  if (df12$RLr_MLpr_WS[x] == "nonLC") {
    
    if (is.na(regions_h_tW_HUF[[df12$area_code_l3[x]]])) {
      regions_h_tW_HUF[[df12$area_code_l3[x]]] <- 0 + df12$HUF[x]
    } else {
      regions_h_tW_HUF[[df12$area_code_l3[x]]] <- regions_h_tW_HUF[[df12$area_code_l3[x]]] + df12$HUF[x]  
    }

    if (is.na(regions_h_tW_HUMa[[df12$area_code_l3[x]]])) {
      regions_h_tW_HUMa[[df12$area_code_l3[x]]] <- 0 + df12$HUMa[x]
    } else {
      regions_h_tW_HUMa[[df12$area_code_l3[x]]] <- regions_h_tW_HUMa[[df12$area_code_l3[x]]] + df12$HUMa[x]  
    }   

    if (is.na(regions_h_tW_HUMe[[df12$area_code_l3[x]]])) {
      regions_h_tW_HUMe[[df12$area_code_l3[x]]] <- 0 + df12$HUMe[x]
    } else {
      regions_h_tW_HUMe[[df12$area_code_l3[x]]] <- regions_h_tW_HUMe[[df12$area_code_l3[x]]] + df12$HUMe[x]  
    }    
    
    if (is.na(regions_h_tW_HUC[[df12$area_code_l3[x]]])) {
      regions_h_tW_HUC[[df12$area_code_l3[x]]] <- 0 + df12$HUC[x]
    } else {
      regions_h_tW_HUC[[df12$area_code_l3[x]]] <- regions_h_tW_HUC[[df12$area_code_l3[x]]] + df12$HUC[x]  
    }
    
    if (is.na(regions_h_tW_HUA[[df12$area_code_l3[x]]])) {
      regions_h_tW_HUA[[df12$area_code_l3[x]]] <- 0 + df12$HUA[x]
    } else {
      regions_h_tW_HUA[[df12$area_code_l3[x]]] <- regions_h_tW_HUA[[df12$area_code_l3[x]]] + df12$HUA[x]  
    }    
    
    if (is.na(regions_h_tW_FD2[[df12$area_code_l3[x]]])) {
      regions_h_tW_FD2[[df12$area_code_l3[x]]] <- 0 + df12$FD2[x]
    } else {
      regions_h_tW_FD2[[df12$area_code_l3[x]]] <- regions_h_tW_FD2[[df12$area_code_l3[x]]] + df12$FD2[x]  
    }    

    if (is.na(regions_h_tW_ED2[[df12$area_code_l3[x]]])) {
      regions_h_tW_ED2[[df12$area_code_l3[x]]] <- 0 + df12$ED2[x]
    } else {
      regions_h_tW_ED2[[df12$area_code_l3[x]]] <- regions_h_tW_ED2[[df12$area_code_l3[x]]] + df12$ED2[x]  
    }
    
    if (is.na(regions_h_tW_COM[[df12$area_code_l3[x]]])) {
      regions_h_tW_COM[[df12$area_code_l3[x]]] <- 0 + df12$COM[x]
    } else {
      regions_h_tW_COM[[df12$area_code_l3[x]]] <- regions_h_tW_COM[[df12$area_code_l3[x]]] + df12$COM[x]  
    }    
    
   }

  if (df12$RLr[x] == "nonLC") {
    
    if (is.na(regions_h_r_HUF[[df12$area_code_l3[x]]])) {
      regions_h_r_HUF[[df12$area_code_l3[x]]] <- 0 + df12$HUF[x]
    } else {
      regions_h_r_HUF[[df12$area_code_l3[x]]] <- regions_h_r_HUF[[df12$area_code_l3[x]]] + df12$HUF[x]  
    }
    
    if (is.na(regions_h_r_HUMa[[df12$area_code_l3[x]]])) {
      regions_h_r_HUMa[[df12$area_code_l3[x]]] <- 0 + df12$HUMa[x]
    } else {
      regions_h_r_HUMa[[df12$area_code_l3[x]]] <- regions_h_r_HUMa[[df12$area_code_l3[x]]] + df12$HUMa[x]  
    }   
    
    if (is.na(regions_h_r_HUMe[[df12$area_code_l3[x]]])) {
      regions_h_r_HUMe[[df12$area_code_l3[x]]] <- 0 + df12$HUMe[x]
    } else {
      regions_h_r_HUMe[[df12$area_code_l3[x]]] <- regions_h_r_HUMe[[df12$area_code_l3[x]]] + df12$HUMe[x]  
    }    
    
    if (is.na(regions_h_r_HUC[[df12$area_code_l3[x]]])) {
      regions_h_r_HUC[[df12$area_code_l3[x]]] <- 0 + df12$HUC[x]
    } else {
      regions_h_r_HUC[[df12$area_code_l3[x]]] <- regions_h_r_HUC[[df12$area_code_l3[x]]] + df12$HUC[x]  
    }
    
    if (is.na(regions_h_r_HUA[[df12$area_code_l3[x]]])) {
      regions_h_r_HUA[[df12$area_code_l3[x]]] <- 0 + df12$HUA[x]
    } else {
      regions_h_r_HUA[[df12$area_code_l3[x]]] <- regions_h_r_HUA[[df12$area_code_l3[x]]] + df12$HUA[x]  
    }    
    
    if (is.na(regions_h_r_FD2[[df12$area_code_l3[x]]])) {
      regions_h_r_FD2[[df12$area_code_l3[x]]] <- 0 + df12$FD2[x]
    } else {
      regions_h_r_FD2[[df12$area_code_l3[x]]] <- regions_h_r_FD2[[df12$area_code_l3[x]]] + df12$FD2[x]  
    }    
    
    if (is.na(regions_h_r_ED2[[df12$area_code_l3[x]]])) {
      regions_h_r_ED2[[df12$area_code_l3[x]]] <- 0 + df12$ED2[x]
    } else {
      regions_h_r_ED2[[df12$area_code_l3[x]]] <- regions_h_r_ED2[[df12$area_code_l3[x]]] + df12$ED2[x]  
    }
    
    if (is.na(regions_h_r_COM[[df12$area_code_l3[x]]])) {
      regions_h_r_COM[[df12$area_code_l3[x]]] <- 0 + df12$COM[x]
    } else {
      regions_h_r_COM[[df12$area_code_l3[x]]] <- regions_h_r_COM[[df12$area_code_l3[x]]] + df12$COM[x]  
    }    
    
  }
}



# combine our data to the other shapefile data and populate with scores

tdwg3$LEVEL3_SPE_HUF_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUA_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_FD2_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_ED2_B <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_COM_B <- rep(NA, length(tdwg3$LEVEL3_NAM))

tdwg3$LEVEL3_SPE_HUF_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUA_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_FD2_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_ED2_W <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_COM_W <- rep(NA, length(tdwg3$LEVEL3_NAM))

tdwg3$LEVEL3_SPE_HUF_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUA_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_FD2_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_ED2_R <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_COM_R <- rep(NA, length(tdwg3$LEVEL3_NAM))

# add columns to discretise it in case nicer for maps
tdwg3$LEVEL3_SPE_HUF_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_B2<- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUA_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_FD2_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_ED2_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_COM_B2 <- rep(NA, length(tdwg3$LEVEL3_NAM))

tdwg3$LEVEL3_SPE_HUF_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUA_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_FD2_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_ED2_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_COM_W2 <- rep(NA, length(tdwg3$LEVEL3_NAM))

tdwg3$LEVEL3_SPE_HUF_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUA_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_FD2_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_ED2_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_COM_R2 <- rep(NA, length(tdwg3$LEVEL3_NAM))



for (x in 1:length(tdwg3$LEVEL3_COD)) {

  if (is.na(regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUF_B[x] <- regions_h_tB_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUF_B2[x] <- NA
  } else {
  if (regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) { # needed to distinguish places without palms (NA) and places without concerned, eg HUF, palms (NaN/NOC)
    tdwg3$LEVEL3_SPE_HUF_B[x] <- NaN
    tdwg3$LEVEL3_SPE_HUF_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUF_B[x] <- regions_h_tB_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]
  if (tdwg3$LEVEL3_SPE_HUF_B[x] < 20) {tdwg3$LEVEL3_SPE_HUF_B2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUF_B[x] < 40) {tdwg3$LEVEL3_SPE_HUF_B2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUF_B[x] < 60) {tdwg3$LEVEL3_SPE_HUF_B2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUF_B[x] < 80) {tdwg3$LEVEL3_SPE_HUF_B2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUF_B2[x] <- "80_100"
      }
    } 
    }
}
  }
  }
  
  if (is.na(regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUMa_B[x] <- regions_h_tB_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUMa_B2[x] <- NA
  } else {
  if (regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUMa_B[x] <- NaN
    tdwg3$LEVEL3_SPE_HUMa_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUMa_B[x] <- regions_h_tB_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUMa_B[x] < 20) {tdwg3$LEVEL3_SPE_HUMa_B2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUMa_B[x] < 40) {tdwg3$LEVEL3_SPE_HUMa_B2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUMa_B[x] < 60) {tdwg3$LEVEL3_SPE_HUMa_B2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUMa_B[x] < 80) {tdwg3$LEVEL3_SPE_HUMa_B2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUMa_B2[x] <- "80_100"
        }
      }
    }
    }
  }
  }
  
  if (is.na(regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUMe_B[x] <- regions_h_tB_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUMe_B2[x] <- NA
  } else {
  if (regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUMe_B[x] <- NaN
    tdwg3$LEVEL3_SPE_HUMe_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUMe_B[x] <- regions_h_tB_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]
   if (tdwg3$LEVEL3_SPE_HUMe_B[x] < 20) {tdwg3$LEVEL3_SPE_HUMe_B2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUMe_B[x] < 40) {tdwg3$LEVEL3_SPE_HUMe_B2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUMe_B[x] < 60) {tdwg3$LEVEL3_SPE_HUMe_B2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUMe_B[x] < 80) {tdwg3$LEVEL3_SPE_HUMe_B2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUMe_B2[x] <- "80_100"
        }
      }
    }
   }
  }
  }
  
  if (is.na(regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUC_B[x] <- regions_h_tB_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUC_B2[x] <- NA
  } else {
  if (regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUC_B[x] <- NaN
    tdwg3$LEVEL3_SPE_HUC_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUC_B[x] <- regions_h_tB_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUC_B[x] < 20) {tdwg3$LEVEL3_SPE_HUC_B2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUC_B[x] < 40) {tdwg3$LEVEL3_SPE_HUC_B2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUC_B[x] < 60) {tdwg3$LEVEL3_SPE_HUC_B2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUC_B[x] < 80) {tdwg3$LEVEL3_SPE_HUC_B2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUC_B2[x] <- "80_100"
        }
      }
    }
    }
  }
  }
  
  if (is.na(regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUA_B[x] <- regions_h_tB_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUA_B2[x] <- NA
  } else {
  if (regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUA_B[x] <- NaN
    tdwg3$LEVEL3_SPE_HUA_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUA_B[x] <- regions_h_tB_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUA_B[x] < 20) {tdwg3$LEVEL3_SPE_HUA_B2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUA_B[x] < 40) {tdwg3$LEVEL3_SPE_HUA_B2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUA_B[x] < 60) {tdwg3$LEVEL3_SPE_HUA_B2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUA_B[x] < 80) {tdwg3$LEVEL3_SPE_HUA_B2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUA_B2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_FD2_B[x] <- regions_h_tB_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_FD2_B2[x] <- NA
  } else {
  if (regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_FD2_B[x] <- NaN
    tdwg3$LEVEL3_SPE_FD2_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_FD2_B[x] <- regions_h_tB_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_FD2_B[x] < 20) {tdwg3$LEVEL3_SPE_FD2_B2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_FD2_B[x] < 40) {tdwg3$LEVEL3_SPE_FD2_B2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_FD2_B[x] < 60) {tdwg3$LEVEL3_SPE_FD2_B2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_FD2_B[x] < 80) {tdwg3$LEVEL3_SPE_FD2_B2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_FD2_B2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_ED2_B[x] <- regions_h_tB_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_ED2_B2[x] <- NA
  } else {
  if (regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_ED2_B[x] <- NaN
    tdwg3$LEVEL3_SPE_ED2_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_ED2_B[x] <- regions_h_tB_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_ED2_B[x] < 20) {tdwg3$LEVEL3_SPE_ED2_B2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_ED2_B[x] < 40) {tdwg3$LEVEL3_SPE_ED2_B2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_ED2_B[x] < 60) {tdwg3$LEVEL3_SPE_ED2_B2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_ED2_B[x] < 80) {tdwg3$LEVEL3_SPE_ED2_B2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_ED2_B2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_COM_B[x] <- regions_h_tB_COM[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_COM_B2[x] <- NA
  } else {
  if (regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_COM_B[x] <- NaN
    tdwg3$LEVEL3_SPE_COM_B2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_COM_B[x] <- regions_h_tB_COM[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_COM_B[x] < 20) {tdwg3$LEVEL3_SPE_COM_B2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_COM_B[x] < 40) {tdwg3$LEVEL3_SPE_COM_B2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_COM_B[x] < 60) {tdwg3$LEVEL3_SPE_COM_B2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_COM_B[x] < 80) {tdwg3$LEVEL3_SPE_COM_B2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_COM_B2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }

  if (is.na(regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUF_W[x] <- regions_h_tW_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUF_W2[x] <- NA
  } else {
  if (regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUF_W[x] <- NaN
    tdwg3$LEVEL3_SPE_HUF_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUF_W[x] <- regions_h_tW_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUF_W[x] < 20) {tdwg3$LEVEL3_SPE_HUF_W2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUF_W[x] < 40) {tdwg3$LEVEL3_SPE_HUF_W2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUF_W[x] < 60) {tdwg3$LEVEL3_SPE_HUF_W2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUF_W[x] < 80) {tdwg3$LEVEL3_SPE_HUF_W2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUF_W2[x] <- "80_100"
        }
      }
    }
    }
  }
  }
  
  if (is.na(regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUMa_W[x] <- regions_h_tW_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUMa_W2[x] <- NA
  } else {
  if (regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUMa_W[x] <- NaN
    tdwg3$LEVEL3_SPE_HUMa_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUMa_W[x] <- regions_h_tW_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUMa_W[x] < 20) {tdwg3$LEVEL3_SPE_HUMa_W2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUMa_W[x] < 40) {tdwg3$LEVEL3_SPE_HUMa_W2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUMa_W[x] < 60) {tdwg3$LEVEL3_SPE_HUMa_W2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUMa_W[x] < 80) {tdwg3$LEVEL3_SPE_HUMa_W2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUMa_W2[x] <- "80_100"
        }
      }
    }
    }
  } 
  }
  if (is.na(regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUMe_W[x] <- regions_h_tW_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUMe_W2[x] <- NA
  } else {
  if (regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUMe_W[x] <- NaN
    tdwg3$LEVEL3_SPE_HUMe_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUMe_W[x] <- regions_h_tW_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUMe_W[x] < 20) {tdwg3$LEVEL3_SPE_HUMe_W2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUMe_W[x] < 40) {tdwg3$LEVEL3_SPE_HUMe_W2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUMe_W[x] < 60) {tdwg3$LEVEL3_SPE_HUMe_W2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUMe_W[x] < 80) {tdwg3$LEVEL3_SPE_HUMe_W2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUMe_W2[x] <- "80_100"
        }
      }
    }
    }
  }  
  }
  
  if (is.na(regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUC_W[x] <- regions_h_tW_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUC_W2[x] <- NA
  } else {
  if (regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUC_W[x] <- NaN
    tdwg3$LEVEL3_SPE_HUC_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUC_W[x] <- regions_h_tW_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUC_W[x] < 20) {tdwg3$LEVEL3_SPE_HUC_W2[x] <- "0_20"} else {
    if (tdwg3$LEVEL3_SPE_HUC_W[x] < 40) {tdwg3$LEVEL3_SPE_HUC_W2[x] <- "20_40"} else {
      if (tdwg3$LEVEL3_SPE_HUC_W[x] < 60) {tdwg3$LEVEL3_SPE_HUC_W2[x] <- "40_60"} else {
        if (tdwg3$LEVEL3_SPE_HUC_W[x] < 80) {tdwg3$LEVEL3_SPE_HUC_W2[x] <- "60_80"} else {
          tdwg3$LEVEL3_SPE_HUC_W2[x] <- "80_100"
        }
      }
    }
  }
  }
  }
  if (is.na(regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUA_W[x] <- regions_h_tW_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUA_W2[x] <- NA
  } else {
  if (regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUA_W[x] <- NaN
    tdwg3$LEVEL3_SPE_HUA_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUA_W[x] <- regions_h_tW_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUA_W[x] < 20) {tdwg3$LEVEL3_SPE_HUA_W2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUA_W[x] < 40) {tdwg3$LEVEL3_SPE_HUA_W2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUA_W[x] < 60) {tdwg3$LEVEL3_SPE_HUA_W2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUA_W[x] < 80) {tdwg3$LEVEL3_SPE_HUA_W2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUA_W2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_FD2_W[x] <- regions_h_tW_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_FD2_W2[x] <- NA
  } else {
  if (regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_FD2_W[x] <- NaN
    tdwg3$LEVEL3_SPE_FD2_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_FD2_W[x] <- regions_h_tW_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_FD2_W[x] < 20) {tdwg3$LEVEL3_SPE_FD2_W2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_FD2_W[x] < 40) {tdwg3$LEVEL3_SPE_FD2_W2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_FD2_W[x] < 60) {tdwg3$LEVEL3_SPE_FD2_W2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_FD2_W[x] < 80) {tdwg3$LEVEL3_SPE_FD2_W2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_FD2_W2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_ED2_W[x] <- regions_h_tW_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_ED2_W2[x] <- NA
  } else {
  if (regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_ED2_W[x] <- NaN
    tdwg3$LEVEL3_SPE_ED2_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_ED2_W[x] <- regions_h_tW_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_ED2_W[x] < 20) {tdwg3$LEVEL3_SPE_ED2_W2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_ED2_W[x] < 40) {tdwg3$LEVEL3_SPE_ED2_W2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_ED2_W[x] < 60) {tdwg3$LEVEL3_SPE_ED2_W2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_ED2_W[x] < 80) {tdwg3$LEVEL3_SPE_ED2_W2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_ED2_W2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_COM_W[x] <- regions_h_tW_COM[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_COM_W2[x] <- NA
  } else {
  if (regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_COM_W[x] <- NaN
    tdwg3$LEVEL3_SPE_COM_W2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_COM_W[x] <- regions_h_tW_COM[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_COM_W[x] < 20) {tdwg3$LEVEL3_SPE_COM_W2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_COM_W[x] < 40) {tdwg3$LEVEL3_SPE_COM_W2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_COM_W[x] < 60) {tdwg3$LEVEL3_SPE_COM_W2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_COM_W[x] < 80) {tdwg3$LEVEL3_SPE_COM_W2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_COM_W2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  if (is.na(regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUF_R[x] <- regions_h_r_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUF_R2[x] <- NA
  } else {  
  if (regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUF_R[x] <- NaN
    tdwg3$LEVEL3_SPE_HUF_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUF_R[x] <- regions_h_r_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUF[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUF_R[x] < 20) {tdwg3$LEVEL3_SPE_HUF_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUF_R[x] < 40) {tdwg3$LEVEL3_SPE_HUF_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUF_R[x] < 60) {tdwg3$LEVEL3_SPE_HUF_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUF_R[x] < 80) {tdwg3$LEVEL3_SPE_HUF_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUF_R2[x] <- "80_100"
          }
        }
      }
    }
  }
  }
  if (is.na(regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUMa_R[x] <- regions_h_r_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUMa_R2[x] <- NA
  } else {  
  if (regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUMa_R[x] <- NaN
    tdwg3$LEVEL3_SPE_HUMa_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUMa_R[x] <- regions_h_r_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMa[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUMa_R[x] < 20) {tdwg3$LEVEL3_SPE_HUMa_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUMa_R[x] < 40) {tdwg3$LEVEL3_SPE_HUMa_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUMa_R[x] < 60) {tdwg3$LEVEL3_SPE_HUMa_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUMa_R[x] < 80) {tdwg3$LEVEL3_SPE_HUMa_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUMa_R2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  if (is.na(regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUMe_R[x] <- regions_h_r_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUMe_R2[x] <- NA
  } else { 
  if (regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUMe_R[x] <- NaN
    tdwg3$LEVEL3_SPE_HUMe_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUMe_R[x] <- regions_h_r_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUMe[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUMe_R[x] < 20) {tdwg3$LEVEL3_SPE_HUMe_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUMe_R[x] < 40) {tdwg3$LEVEL3_SPE_HUMe_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUMe_R[x] < 60) {tdwg3$LEVEL3_SPE_HUMe_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUMe_R[x] < 80) {tdwg3$LEVEL3_SPE_HUMe_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUMe_R2[x] <- "80_100"
          }
        }
      }
    }
  } 
  }
  if (is.na(regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUC_R[x] <- regions_h_r_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUC_R2[x] <- NA
  } else { 
  if (regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUC_R[x] <- NaN
    tdwg3$LEVEL3_SPE_HUC_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUC_R[x] <- regions_h_r_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUC[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUC_R[x] < 20) {tdwg3$LEVEL3_SPE_HUC_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUC_R[x] < 40) {tdwg3$LEVEL3_SPE_HUC_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUC_R[x] < 60) {tdwg3$LEVEL3_SPE_HUC_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUC_R[x] < 80) {tdwg3$LEVEL3_SPE_HUC_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUC_R2[x] <- "80_100"
          }
        }
      }
    }
  }
  }
  if (is.na(regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_HUA_R[x] <- regions_h_r_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_HUA_R2[x] <- NA
  } else { 
  if (regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_HUA_R[x] <- NaN
    tdwg3$LEVEL3_SPE_HUA_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_HUA_R[x] <- regions_h_r_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_HUA[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_HUA_R[x] < 20) {tdwg3$LEVEL3_SPE_HUA_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_HUA_R[x] < 40) {tdwg3$LEVEL3_SPE_HUA_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_HUA_R[x] < 60) {tdwg3$LEVEL3_SPE_HUA_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_HUA_R[x] < 80) {tdwg3$LEVEL3_SPE_HUA_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_HUA_R2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  if (is.na(regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_FD2_R[x] <- regions_h_r_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_FD2_R2[x] <- NA
  } else { 
  if (regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_FD2_R[x] <- NaN
    tdwg3$LEVEL3_SPE_FD2_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_FD2_R[x] <- regions_h_r_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_FD2[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_FD2_R[x] < 20) {tdwg3$LEVEL3_SPE_FD2_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_FD2_R[x] < 40) {tdwg3$LEVEL3_SPE_FD2_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_FD2_R[x] < 60) {tdwg3$LEVEL3_SPE_FD2_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_FD2_R[x] < 80) {tdwg3$LEVEL3_SPE_FD2_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_FD2_R2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  if (is.na(regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_ED2_R[x] <- regions_h_r_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_ED2_R2[x] <- NA
  } else { 
  if (regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_ED2_R[x] <- NaN
    tdwg3$LEVEL3_SPE_ED2_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_ED2_R[x] <- regions_h_r_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_ED2[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_ED2_R[x] < 20) {tdwg3$LEVEL3_SPE_ED2_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_ED2_R[x] < 40) {tdwg3$LEVEL3_SPE_ED2_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_ED2_R[x] < 60) {tdwg3$LEVEL3_SPE_ED2_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_ED2_R[x] < 80) {tdwg3$LEVEL3_SPE_ED2_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_ED2_R2[x] <- "80_100"
          }
        }
      }
    }
  } 
  }
  if (is.na(regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]) == T) {
    tdwg3$LEVEL3_SPE_COM_R[x] <- regions_h_r_COM[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]
    tdwg3$LEVEL3_SPE_COM_R2[x] <- NA
  } else { 
  if (regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]] == 0) {
    tdwg3$LEVEL3_SPE_COM_R[x] <- NaN
    tdwg3$LEVEL3_SPE_COM_R2[x] <- "NOC"
  } else {
  tdwg3$LEVEL3_SPE_COM_R[x] <- regions_h_r_COM[[as.character(tdwg3$LEVEL3_COD[x])]]*100/regions_h_a_COM[[as.character(tdwg3$LEVEL3_COD[x])]]
    if (tdwg3$LEVEL3_SPE_COM_R[x] < 20) {tdwg3$LEVEL3_SPE_COM_R2[x] <- "0_20"} else {
      if (tdwg3$LEVEL3_SPE_COM_R[x] < 40) {tdwg3$LEVEL3_SPE_COM_R2[x] <- "20_40"} else {
        if (tdwg3$LEVEL3_SPE_COM_R[x] < 60) {tdwg3$LEVEL3_SPE_COM_R2[x] <- "40_60"} else {
          if (tdwg3$LEVEL3_SPE_COM_R[x] < 80) {tdwg3$LEVEL3_SPE_COM_R2[x] <- "60_80"} else {
            tdwg3$LEVEL3_SPE_COM_R2[x] <- "80_100"
          }
        }
      }
    }
  }  
  }
  
  }



# plot

# viridis colors for 5 categories: "#0D0887FF" "#7E03A8FF" "#CC4678FF" "#F89441FF" "#F0F921FF"
# need to change them manually below depending on what categories are in the data
# check the categories:

summary(as.factor(tdwg3$LEVEL3_SPE_COM_B2))
summary(as.factor(tdwg3$LEVEL3_SPE_COM_W2))

# plot the map for best case scenario
map_COM_B2 <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_COM_B2), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F0F921FF", "grey40"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of \"interesting\" species known or predicted to be threatened")+ 
  ggtitle("")+
  theme(legend.position="bottom")

pdf("Fig2A_Map_interesting_species_B.pdf", 16, 8) 
map_COM_B2
dev.off()

setEPS()
postscript("Fig2A_Map_interesting_species_B.eps", height = 8, width = 16)
map_COM_B2
dev.off()


# same for worst case scenario
map_COM_W2 <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_COM_W2), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF", "grey40"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of \"interesting\" species known or predicted to be threatened")+ 
  ggtitle("")+
  theme(legend.position="bottom")

pdf("Fig2B_Map_interesting_species_W.pdf", 16, 8) 
map_COM_W2
dev.off()

setEPS()
postscript("Fig2B_Map_interesting_species_W.eps", height = 8, width = 16)
map_COM_W2
dev.off()



# Map research gaps [not used, but run part to get the data just in case, but not the plot]

# new columns with 0 if predicted/known status versus 1 if NA
df12$GAP_RL <- mapvalues(as.factor(df12$RLr), from = levels(as.factor(df12$RLr)), to = c(0,1,0))
df12$GAP_RLML <- mapvalues(as.factor(df12$RLr_MLpr_WS), from = levels(as.factor(df12$RLr_MLpr_WS)), to = c(0,1,0))

df12$GAP_RL <- as.numeric(as.character(df12$GAP_RL))
df12$GAP_RLML <- as.numeric(as.character(df12$GAP_RLML))

# dictionaries per regions

regions_h_g_RL <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
regions_h_g_RLML <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))

for (x in 1:length(df12$species)) {
  
  if (is.na(regions_h_g_RL[[df12$area_code_l3[x]]])) {
    regions_h_g_RL[[df12$area_code_l3[x]]] <- 0 + df12$GAP_RL[x]
  } else {
    regions_h_g_RL[[df12$area_code_l3[x]]] <- regions_h_g_RL[[df12$area_code_l3[x]]] + df12$GAP_RL[x]  
  }
  
  if (is.na(regions_h_g_RLML[[df12$area_code_l3[x]]])) {
    regions_h_g_RLML[[df12$area_code_l3[x]]] <- 0 + df12$GAP_RLML[x]
  } else {
    regions_h_g_RLML[[df12$area_code_l3[x]]] <- regions_h_g_RLML[[df12$area_code_l3[x]]] + df12$GAP_RLML[x]  
  }
  
}

# fill raster data table

tdwg3$LEVEL3_SPE_GAP_RL <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_GAP_RLML <- rep(NA, length(tdwg3$LEVEL3_NAM))

for (x in 1:length(tdwg3$LEVEL3_COD)) {
  tdwg3$LEVEL3_SPE_GAP_RL[x] <- regions_h_g_RL[[as.character(tdwg3$LEVEL3_COD[x])]]
  tdwg3$LEVEL3_SPE_GAP_RLML[x] <- regions_h_g_RLML[[as.character(tdwg3$LEVEL3_COD[x])]]
}

# # maps (continuous scale) [not used, not run]
# 
# min(tdwg3$LEVEL3_SPE_GAP_RL, na.rm = T)
# max(tdwg3$LEVEL3_SPE_GAP_RL, na.rm = T)
# min(tdwg3$LEVEL3_SPE_GAP_RLML, na.rm = T)
# max(tdwg3$LEVEL3_SPE_GAP_RLML, na.rm = T)
# 
# values <- seq(0, 265, length = 265) # 265 is the max score
# 
# map_GAP_RL <- ggplot(tdwg3) +
#   geom_sf(aes(fill = LEVEL3_SPE_GAP_RL), color=NA) + # remove color=NA if want country borders
#   scale_fill_viridis_c(option= "C", values = values, na.value = "grey", rescaler = function(x, ...) x) +
#   theme_classic() + 
#   labs(fill = "Number of species")+ 
#   ggtitle("Species without a recent published assessment")+
#   theme(legend.position="bottom")
# 
# map_GAP_RLML <- ggplot(tdwg3) +
#   geom_sf(aes(fill = LEVEL3_SPE_GAP_RLML), color=NA) + # remove color=NA if want country borders
#   scale_fill_viridis_c(option= "C", values = values, na.value = "grey", rescaler = function(x, ...) x) +
#   theme_classic() + 
#   labs(fill = "Number of species")+ 
#   ggtitle("Species without a recent published assessment or prediction")+
#   theme(legend.position="bottom")
# 
# 
# 
# pdf("FigX_Map_gaps.pdf", 8, 8) 
# grid.arrange(map_GAP_RL, map_GAP_RLML, nrow=2)
# dev.off()
# 
# setEPS()
# postscript("FigX_Map_gaps.eps", height = 8, width = 8)
# grid.arrange(map_GAP_RL, map_GAP_RLML, nrow=2)
# dev.off()

# same with easier distinction along scale when considering RLML only

# values <- seq(0, 102, length = 102) # 265 is the max score
# 
# map_GAP_RLML2 <- ggplot(tdwg3) +
#   geom_sf(aes(fill = LEVEL3_SPE_GAP_RLML), color=NA) + # remove color=NA if want country borders
#   scale_fill_viridis_c(option= "C", values = values, na.value = "grey", rescaler = function(x, ...) x) +
#   theme_classic() + 
#   labs(fill = "Number of species")+ 
#   ggtitle("Species without a recent published assessment or prediction")+
#   theme(legend.position="bottom")
# 
# 
# pdf("FigX_Map_gaps_RLMLonly.pdf", 8, 4) 
# map_GAP_RLML2
# dev.off()
# 
# setEPS()
# postscript("FigX_Map_gaps_RLMLonly.eps", height = 4, width = 8)
# map_GAP_RLML2
# dev.off()



# Bump charts
#"#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"

# add a species number column to the data (do it on the raster and then generate separate table)

regions_h_spe <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))

for (x in 1:length(df12$species)) {
  
  if (is.na(regions_h_spe[[df12$area_code_l3[x]]])) {
    regions_h_spe[[df12$area_code_l3[x]]] <- 1
  } else {
    regions_h_spe[[df12$area_code_l3[x]]] <- regions_h_spe[[df12$area_code_l3[x]]] + 1  
  }

}

tdwg3$LEVEL3_SPE_speNum <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_speNumCat <- rep(NA, length(tdwg3$LEVEL3_NAM))

for (x in 1:length(tdwg3$LEVEL3_COD)) {
  tdwg3$LEVEL3_SPE_speNum[x] <- regions_h_spe[[as.character(tdwg3$LEVEL3_COD[x])]]
  if (is.na(tdwg3$LEVEL3_SPE_speNum[x]) == T) {
    tdwg3$LEVEL3_SPE_speNumCat[x] <- NA
  } else {
  if (tdwg3$LEVEL3_SPE_speNum[x] > 100) {
    tdwg3$LEVEL3_SPE_speNumCat[x] <- ">100"
  } else {
    if (tdwg3$LEVEL3_SPE_speNum[x] > 50) {
      tdwg3$LEVEL3_SPE_speNumCat[x] <- "51-100"
    } else {
      if (tdwg3$LEVEL3_SPE_speNum[x] > 10) {
        tdwg3$LEVEL3_SPE_speNumCat[x] <- "11-50"
      } else {
          tdwg3$LEVEL3_SPE_speNumCat[x] <- "1-10"
      }     
    }
  }
  }
}

# map just to see if results make sense
ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_speNumCat), color=NA) + # remove color=NA if want country borders
 # scale_fill_viridis_c(option= "C", values = values, na.value = "grey", rescaler = function(x, ...) x) +
  theme_classic()



# export data

data_labs <- as.data.frame(tdwg3$LEVEL3_NAM)
data_labs$LEVEL3_COD <- tdwg3$LEVEL3_COD
data_labs$LEVEL2_COD <- tdwg3$LEVEL2_COD
data_labs$LEVEL1_COD <- tdwg3$LEVEL1_COD
for (z in 6:length(names(tdwg3))) {
  y <- z - 1
  data_labs[,y] <- tdwg3[,z][[1]]
}
names(data_labs)[5:length(names(data_labs))] <- names(tdwg3)[6:length(names(tdwg3))]
names(data_labs)[1] <- "LEVEL3_NAM"
write.table(data_labs, "Uses_by_region_Figs1Band3_Realtable.txt", sep="\t", row.names = F)

# in this table, NOC indicates regions with palms but no concerned palms while NA indicates regions without palms (142)
# Table is analyzed in Uses_by_region_Figs1Band3_Realtable_digest.xlsx to produce Table S4(?)
# [did not regenerate for version 6 but instead renamed original table and digested tables from version 4, since original of version 6 exactly identical to original of version 4]
# Could also be used to produce a sup table with info subtending figures


# Bump charts for species of interest

# based on the table, the top45 regions of COM_W will cover the top 10 regions according to all 3 scenarios (R, B, W)
# reorder the data based on COM_W and keep only the top 45 or the whole data (do two plots, one with just these and one with all)
# only using the plot with all for the paper (Fig. S3)

data_labs_COM_W <- data_labs[ order(-data_labs$LEVEL3_SPE_COM_W), ]
data_labs_COM_W_T45 <- data_labs_COM_W[1:45,]

# order also COM_R and COM_10 to get labels for regions in the top 10, i.e. top 32 for W and top 15 for B (due to ties - I know that by looking at the table in excel)
data_labs_COM_W_T32 <- data_labs_COM_W[1:32,]

data_labs_COM_R <- data_labs[ order(-data_labs$LEVEL3_SPE_COM_R), ]
data_labs_COM_R_T15 <- data_labs_COM_R[1:15,]

summary(as.factor(data_labs_COM_W_T45$LEVEL3_SPE_speNumCat))  # change to color by species richness

p <- ggparcoord(data_labs_COM_W, 
                  columns = c(28, 12, 20), 
                  groupColumn = 56,
                  scale = "globalminmax",
                  showPoints = T,
                  title = "Region ranking based on proportion of \"valuable\" species that are threatened",
                  alphaLines = 1) + 
  scale_color_manual(values = c("#D55E00", "#009E73" , "#0072B2", "#E69F00"))+
  theme_classic()+
  #geom_text_repel(data = data_labs_COM_W, aes(x = 1, y = LEVEL3_SPE_COM_R, label = LEVEL3_NAM),
                 # inherit.aes = FALSE, xlim = c(NA, 1), max.overlaps = 29, min.segment.length = 0)+
  #geom_text_repel(data = data_labs_COM_W, aes(x = 3, y = LEVEL3_SPE_COM_W, label = LEVEL3_NAM),
                  #inherit.aes = FALSE, xlim = c(3,NA), max.overlaps = 29, min.segment.length = 0)+
  theme(legend.position = "none") +
  labs(x= "", y = "Percentage of threatened species")


p45 <- ggparcoord(data_labs_COM_W_T45, 
                columns = c(28, 12, 20), 
                groupColumn = 56,   
                scale = "globalminmax",
                showPoints = T,
                title = "Region ranking based on proportion of \"valuable\" species that are threatened",
                alphaLines = 1) + 
  scale_color_manual(values = c("#D55E00", "#009E73" , "#0072B2"))+
  theme_classic()+
  geom_text_repel(data = data_labs_COM_R_T15, aes(x = 1, y = LEVEL3_SPE_COM_R, label = LEVEL3_NAM),
                  inherit.aes = FALSE, xlim = c(NA, 1), max.overlaps = 15, min.segment.length = 0)+
  geom_text_repel(data = data_labs_COM_W_T32, aes(x = 3, y = LEVEL3_SPE_COM_W, label = LEVEL3_NAM),
                  inherit.aes = FALSE, xlim = c(3,NA), max.overlaps = 32, min.segment.length = 0)+
  theme(legend.position = "none") +
  labs(x= "", y = "Percentage of threatened species")


#pdf("Bump_valuable_top10labs.pdf", 8, 8) 
#p45
#dev.off()
#setEPS()
#postscript("Bump_valuable_top10labs.eps", height = 8, width = 8)
#p45
#dev.off()  

pdf("FigS3_Bump_valuable_full.pdf", 8, 16) 
p
dev.off()

setEPS()
postscript("FigS3_Bump_valuable_full.eps", height = 16, width = 8)
p
dev.off()  


##########################


###############################
# Proportions and numbers of threatened species
# Multiple biodiversity facets
# Regions bump charts
###############################

### plot proportions of threatened species for each use category, color lines by continent for few species rich regions
# Then do same as OSM with all regions, and coloring by species richness 
# do it separately for best and worse


# Create new column manually in the data table with continent of the region 
# (based on what a continent is for us, done manually in Continent_coding sheet of Uses_by_region_Figs1Band3_Realtable_digest.xlsx
# and re-exported to Uses_by_region_Figs1Band3_Realtable_wContinent.txt)
# Manual coding was done based on level 1 of tdwg coding except where:
# 1 and 2 = AfricaWAsia (fits our delimitation of species continent based on their coordinates between -25 and 68, except for Azores which according to us would go to Americas but not important because no native palms)
# 3 is split between AfricaWAsia and EAsiaPacific at 68d longitude (few regions split in 2: code following where the majority of the region is. Not important because they do not have palms)
# 4 and 5 = EAsiaPacific (fits our delimitation based on coordinates between 68 and 180, except for Pakistan --> coded as AfricaWAsia)
# 6 split between Americas and EastAsiaPacific at 180d longitude
# 7 and 8 = Americas
# 9 split between the 3 regions based on longitude (-25, 68, 180) and with Antartica as America (arbitrary, irrelevant for palms)
# --> our coding is coarser but broadly fits the continent level of TDWG


data <- read.table("Uses_by_region_Figs1Band3_Realtable_wContinent.txt", header=T, stringsAsFactors = T, sep="\t")

# keep only the 5 most palm-rich regions in each continent, color by continent
# subset by continent
data_Am <- subset(data, data$Continent == "Americas")
data_Af <- subset(data, data$Continent == "AfricaWAsia")
data_As <- subset(data, data$Continent == "EAsiaPacific")
# order
data_Am_o <- data_Am[ order(-data_Am$LEVEL3_SPE_speNum), ]
data_Af_o <- data_Af[ order(-data_Af$LEVEL3_SPE_speNum), ]
data_As_o <- data_As[ order(-data_As$LEVEL3_SPE_speNum), ]
# get top 5
data_Am_o_T5 <- data_Am_o[1:5,]
data_Af_o_T5 <- data_Af_o[1:5,]
data_As_o_T5 <- data_As_o[1:5,]
data_Am_o_T10 <- data_Am_o[1:10,]
data_Af_o_T10 <- data_Af_o[1:10,]
data_As_o_T10 <- data_As_o[1:10,]
# cbind all back
data_T5 <- rbind(data_Am_o_T5, data_Af_o_T5, data_As_o_T5)
data_T10 <- rbind(data_Am_o_T10, data_Af_o_T10, data_As_o_T10)



data_labs <- as.data.frame(data_T10$LEVEL3_NAM)
data_labs$LEVEL3_SPE_FD2_B <- data_T10$LEVEL3_SPE_FD2_B
names(data_labs)[1] <- "LEVEL3_NAM"

data_labs_W <- as.data.frame(data_T10$LEVEL3_NAM)
data_labs_W$LEVEL3_SPE_FD2_W <- data_T10$LEVEL3_SPE_FD2_W
names(data_labs_W)[1] <- "LEVEL3_NAM"


# plot, Best-case scenario, color by continent
pB <- ggparcoord(data_T10, 
                 columns = c(10,11,9,5,6,7,8), 
                 groupColumn = 57,
                 scale = "globalminmax",
                 showPoints = T,
                 #title = "Region ranking based on proportion of HUF species that are threatened",
                 alphaLines = 1) + 
  scale_color_manual(values = c("#0D0887FF", "#CC4678FF", "#F89441FF"))+
  geom_text_repel(data = data_labs, aes(x = 1, y = LEVEL3_SPE_FD2_B, label = LEVEL3_NAM), size=3,
  inherit.aes = FALSE, xlim = c(NA, 1), max.overlaps = 20, min.segment.length = 0)+
  theme_classic()+
  theme(legend.position = "top") +
  labs(x= "", y = "Percentage of threatened species")


# plot, Worst-case scenario, color by continent

pW <- ggparcoord(data_T10, 
                 columns = c(18,19,17,13,14,15,16), 
                 groupColumn = 57,
                 scale = "globalminmax",
                 showPoints = T,
                 #title = "Region ranking based on proportion of HUF species that are threatened",
                 alphaLines = 1) + 
  scale_color_manual(values = c("#0D0887FF", "#CC4678FF", "#F89441FF"))+
  geom_text_repel(data = data_labs_W, aes(x = 1, y = LEVEL3_SPE_FD2_W, label = LEVEL3_NAM), size=3,
                  inherit.aes = FALSE, xlim = c(NA, 1), max.overlaps = 20, min.segment.length = 0)+
  theme_classic()+
  theme(legend.position = "top") +
  labs(x= "", y = "Percentage of threatened species")


pdf("Fig1B_Bump_uses_T10_2.pdf", 8, 24) 
grid.arrange(pB, pW, nrow=2)
dev.off()

setEPS()
postscript("Fig1B_Bump_uses_T10.eps", height = 24, width = 8)
grid.arrange(pB, pW, nrow=2)
dev.off() 



# plot all regions, Best-case scenario, color by species richness
pB <- ggparcoord(data, 
columns = c(10,11,9,5,6,7,8), 
groupColumn = 56,
scale = "globalminmax",
showPoints = T,
#title = "Region ranking based on proportion of HUF species that are threatened",
alphaLines = 1) + 
#scale_color_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF"))+
theme_classic()+
theme(legend.position = "top") +
labs(x= "", y = "Percentage of threatened species")


# plot all regions, Worst-case scenario, color by species richness

pW <- ggparcoord(data, 
                columns = c(18,19,17,13,14,15,16), 
                groupColumn = 56,
                scale = "globalminmax",
                showPoints = T,
                #title = "Region ranking based on proportion of HUF species that are threatened",
                alphaLines = 1) + 
  #scale_color_manual(values = c("#0D0887FF", "#7E03A8FF", "#CC4678FF"))+
  theme_classic()+
  theme(legend.position = "top") +
  labs(x= "", y = "Percentage of threatened species")


pdf("FigS2C_Bump_uses.pdf", 8, 24) 
grid.arrange(pB, pW, nrow=2)
dev.off()
 
setEPS()
postscript("FigS2C_Bump_uses.eps", height = 24, width = 8)
grid.arrange(pB, pW, nrow=2)
dev.off() 

##################################




#######################################
# plot uses on consensus tree (Fig. S4 - part 1)
#######################################

# import old table with tree name equivalent to species name
d <- read.table("Phylogeny_NoCon_Henderson_Tree_1_taxa_forPy_CONS-CONT-TRAITSALL-USES_lessCol.txt", header=T, stringsAsFactors = F, sep="\t")

# import a maximum clade credibility tree made with Tree Annotator on trees 251-1000 of Phylogeny_NoCon_Henderson_AllTreesCorrected
tree_ed <- read.nexus("Phylogeny_NoCon_Henderson_AllTreesCorrected_NEX_TA.nex")
tree_ed2 <- root(tree_ed, node=2489) # root on Calamoideae
tree_ed3 <- ladderize(tree_ed2)

# ensure taxon names are same as in tree
d$Tree_name <- gsub(" ", "_", d$Tree_name) # taxon labels have to be the first column of the attached data! (does not matter below because we do not use %<+%)

# make a subset with the uses only (easier later)
d_sub <- d[,c(36, 38, 39,41)]
rownames(d_sub) <- d$Tree_name

# code each column differently so that can color them differently (1 to 8, odd corresponding to no use, even to use), with always black for odd/no use

# replace NA by 0 for the uses
NA_rep_0 <- function(table) {
  for ( x in 1:length(table[,1])) {
    for (y in 1:length(names(table))) {
      if (is.na(table[x,y])) {
        table[x,y] <- 0
      }
    }
  }
  return(table)
}


d_sub2 <- NA_rep_0(d_sub)


# recode to have different colors
d_sub2$HumanFood <- mapvalues(as.factor(d_sub2$HumanFood), from = levels(as.factor(d_sub2$HumanFood)), to = c("1", "2"))
d_sub2$Materials <- mapvalues(as.factor(d_sub2$Materials), from = levels(as.factor(d_sub2$Materials)), to = c("3", "4"))
d_sub2$Medicines <- mapvalues(as.factor(d_sub2$Medicines), from = levels(as.factor(d_sub2$Medicines)), to = c("5", "6"))
d_sub2$SocialUses <- mapvalues(as.factor(d_sub2$SocialUses), from = levels(as.factor(d_sub2$SocialUses)), to = c("7", "8"))


# transform d_sub into a character df
rn <- rownames(d_sub2)
d_sub2 <- as.data.frame(sapply(d_sub2, as.character))
rownames(d_sub2) <- rn

# "#f0f0f0", "#004529", "#f0f0f0", "#238443", "#f0f0f0", "#78c679", "#f0f0f0", "#fed976"
##fcc5c0

heatmap.colours <- c("gray96", "#49006a", "gray96", "#ae017e", "gray96", "#f768a1", "gray96", "salmon")
names(heatmap.colours) <- 1:8

p <- ggtree(tree_ed3, layout='circular', branch.length = "none", size=0.1)

p2 <- p + geom_tiplab2(align=T, linetype="dotted", size=0.2, hjust=-0.02)

ph <- gheatmap(p, d_sub2, offset = 0, width=0.3, color=NULL, colnames=F)+
  #colnames_position="top", 
  # colnames_angle=90, colnames_offset_y = 1, 
  #hjust=0, font.size=2) +
  scale_fill_manual(values=heatmap.colours,breaks=1:8)

ph2 <- gheatmap(p2, d_sub2, offset = 1, width=0.15, color=NULL, colnames=F)+
  #colnames_position="top", 
  # colnames_angle=90, colnames_offset_y = 1, 
  #hjust=0, font.size=2) +
  scale_fill_manual(values=heatmap.colours,breaks=1:8)


# pdf("FigX_Uses_on_tree.pdf", 10, 10) 
# ph
# dev.off()
# 
# setEPS()
# postscript("FigX_Uses_on_tree.eps", height = 10, width = 10)
# ph
# dev.off()

pdf("FigS4_Uses_on_tree_wTipLabs.pdf", 20, 20) 
ph2
dev.off()

setEPS()
postscript("FigS4_Uses_on_tree_wTipLabs.eps", height = 20, width = 20)
ph2
dev.off()

######################################







#################################
# Calculate probability of being a good substitute for each species, for each use category
#################################

# First make a list with PD dictionaries for all species, but using accepted names

# get distance table
pd1 <- read.table("PhyloDist_Average100Trees.txt", header=T, stringsAsFactors = F, sep="\t")

# logtransform after adding 1 to all (so that 0 does not give -Inf), then rescale all to be between 0 and 1 (the logtr is mostly useful for FD which are quite skewed, but doing it for both)
pd2 <- pd1 + 1
pd3 <- log(pd2)
pd <- scale(pd3, center = FALSE, scale = rep(max(pd3), length(names(pd3))))
pd <- as.data.frame(pd)

# get tree names and accepted names link
d <- read.table("Phylogeny_NoCon_Henderson_Tree_1_taxa_forPy_CONS-CONT-TRAITSALL-USES_lessCol.txt", header=T, stringsAsFactors = F, sep="\t")
#tax_h <- hash(keys=d$Tree_name, values=d$Acc_name)

# get all accepted species and make a list where we will store a dictionary of PD for each species
# SLOW! NO NEED TO DO IT AGAIN!!! 
# If need to regenerate it, start from the table exported below (PhyloDist_AccNames_scaled.txt) and build a list of hashes where one hash is a row of the df (commands below)

d2 <- read.table("Best_case_Worst_case_table.txt", header=T, stringsAsFactors = F, sep="\t")
PD_l <- vector(mode = "list", length = length(d2$species))
names(PD_l) <- d2$species

# fill the dictionaries
# distances may be NA if the accepted species was not in the tree, even under another name

for (s in 1:length(PD_l)) { # for each accepted species, make a dictionary of distances with other accepted species

  print(s)
  occs <- grep(paste("^",names(PD_l)[s],"$",sep=""), d$Acc_name)
  
  if (length(occs) == 0) { #if the accepted species is not in the tree
    PD_l[[s]] <- NA
  } else { 
    
    if (length(occs) == 1) { # if it is, and only once, create a dictionary of its distance to all other species
      
      h <- hash(keys=names(PD_l), values=rep(NA, length(names(PD_l))))
      
      r <- which(rownames(pd) == gsub(" ", "_", d$Tree_name[occs[1]])) # the row of pd with the distances for that species
      
      for (t in 1: length(names(pd))) { # for each species in the tree
        t_name <- names(pd)[t]
        a_name <- d$Acc_name[grep(paste("^",gsub("_", " ", t_name),"$",sep=""), d$Tree_name)]
        #a_name <- unique(a_name)
        if (!is.na(a_name)) { # if the species has an accepted name (if it does not, we will not include it in the dictionary)
          if (!a_name %in% names(PD_l)[s]) {# if the species is not the one for which we are making the dictionary (if it is, we will put NA as distance)
          p <- grep(paste("^",a_name,"$",sep=""), d$Acc_name)
          if (length(p) == 1) { # if it is the only one to have this accepted name
            h[a_name] <- pd[r,t]
          } else { # if many species have this accepted name, average the distances
            occs2 <- gsub(" ", "_", d$Tree_name[p])
            dist <- c()
            for (o in 1:length(occs2)) {
              dist <- c(dist, pd[r,which(names(pd) == occs2[o])])
            }
            h[a_name] <- mean(dist)
          }
          } else {
            h[a_name] <- NA
          }
        }
      }
      
      PD_l[[s]] <- h
      
    } else { # if the accepted species is more than once in the tree (via synonyms), make a dictionary of distances for each synonym and average their values
    
      H_l <- vector(mode = "list", length=length(occs)) # make a list of temporary dictionaries that will later be averaged
      
      for (r in 1:length(occs)) {
      
        h <- hash(keys=names(PD_l), values=rep(NA, length(names(PD_l))))
        
        r2 <- which(rownames(pd) == gsub(" ", "_", d$Tree_name[occs[r]]))
        
        for (t in 1: length(names(pd))) { # for each species in the tree
          t_name <- names(pd)[t]
          a_name <- d$Acc_name[grep(paste("^",gsub("_", " ", t_name),"$",sep=""), d$Tree_name)]
          #a_name <- unique(a_name)
          if (!is.na(a_name)) { # if the species has an accepted name (if it does not, we will not include it in the dictionary)
            if (!a_name %in% names(PD_l)[s]) {# if the species is not the one for which we are making the dictionary (if it is, we will put NA as distance)
              p <- grep(paste("^",a_name,"$",sep=""), d$Acc_name)
              if (length(p) == 1) { # if it is the only one to have this accepted name
                h[a_name] <- pd[r2,t]
              } else { # if many species have this accepted name, average the distances
                occs2 <- gsub(" ", "_", d$Tree_name[p])
                dist <- c()
                for (o in 1:length(occs2)) {
                  dist <- c(dist, pd[r2,which(names(pd) == occs2[o])])
                }
                h[a_name] <- mean(dist)
              }
            } else {
              h[a_name] <- NA
            }
          }
        }
        H_l[[r]] <- h
      }
    
      h2 <- hash(keys=names(PD_l), values=rep(NA, length(names(PD_l))))
      for (x in 1:length(names(PD_l))) {
        v <- c()
        for (z in 1:length(H_l)) {
          v <- c(v, H_l[[z]][[names(PD_l)[x]]])
        }
        k_mean <- mean(v, na.rm=T)
        h2[names(PD_l)[x]] <- k_mean
      }
      
      PD_l[[s]] <- h2
    }
  }
}


# export the list as a dataframe
PD_m <- matrix(rep(NA,length(names(PD_l))*length(names(PD_l))), ncol=length(names(PD_l)))
PD_df <- as.data.frame(PD_m, row.names = names(PD_l))
names(PD_df) <- keys(PD_l[[1]])

for (x in 1:length(PD_l)) { # super inefficient
  print(x)
  sp <- rownames(PD_df)[which(rownames(PD_df) == names(PD_l)[x])]
  if (!is.na(PD_l[[x]])) {
  for (y in 1:length(names(PD_df))) {
    val <- PD_l[[x]][[names(PD_df)[y]]]
    PD_df[sp,y] <- val
  }
  }
}

write.table(PD_df, "PhyloDist_Average100Trees_AccNames_Scaled.txt", sep="\t", row.names = T)

# Regenerate from the dataframe if needed:
pd2 <- read.table("PhyloDist_Average100Trees_AccNames_Scaled.txt", sep="\t", header=T, row.names = 1)
names(pd2) <- gsub("\\.", " ", names(pd2))
names(pd2)[grep("-", rownames(pd2))] <- c("Bactris ana-juliae","Calamus novae-georgii","Calamus toli-toliensis",
                                          "Chamaedorea castillo-montii","Chamaedorea ernesti-augusti",
                                          "Gaussia gomez-pompae","Raphia palma-pinus","Wallaceodoxa raja-ampat")
# (the hyphen was removed which creates problems below)

# replace the NA of the diagonal by 0
for (n in 1:length(rownames(pd2))) {
  pd2[n,n] <- 0
}

PD_l <- vector(mode = "list", length = length(rownames(pd2)))
names(PD_l) <- rownames(pd2)

for (s in 1:length(rownames(pd2))) {
  print(s)
  h <- hash(keys=names(pd2), values=pd2[s,])
  PD_l[[s]] <- h
}



# For species absent from the tree (i.e. with NA everywhere), get the median distances that their congeneric species have with other species, leave NA if no congeneric species

# save a copy of PD_l just in case
PD_l_orig <- vector(mode="list", length=length(PD_l))
for (x in 1:length(PD_l)) {
  PD_l_orig[[x]] <- copy(PD_l[[x]])
}
names(PD_l_orig) <- names(PD_l)

# and another one
# PD_l_2 is the one we will modify, but we made another copy of PD_l (PD_l_orig) just in case
PD_l_2 <- vector(mode="list", length=length(PD_l))
for (x in 1:length(PD_l)) {
  PD_l_2[[x]] <- copy(PD_l[[x]])
}
names(PD_l_2) <- names(PD_l)

# calculate and input medians [!!!!!!! TOOK ONE WEEK TO RUN BECAUSE OF CALAMUS]
for (s in 1:length(names(PD_l))) {
  if (is.na(PD_l[[names(PD_l)[s]]][[names(PD_l)[1]]])) { # check if NA for first species (will be the case for any species with NAs but the first, which has been checked and does not have NAs)

    print(paste(names(PD_l)[s], "NAspecies"))
    genus <- paste("^", str_split(names(PD_l)[s], " ")[[1]][1], " ", sep = "")
    cogeneric_id <- grep(genus, names(PD_l))
    cogeneric_id <- cogeneric_id[which(!cogeneric_id ==s)] # removing the species itself so that does not generate a 0 distance when s2= s below
    
    for (s2 in 1:length(names(PD_l))) { # for each species, get the median distance from the congeneric species to it, and append that in the dictionary of the species that has the NAs
      
      cogeneric_id3 <- cogeneric_id[which(!cogeneric_id ==s2)] # remove s2 from the congeneric of s to avoid distances of 0, which should only happen for truly identical species
      
      if (length(cogeneric_id3) < 1) { # if there is only s2 among the congeneric of s, then cogeneric_id will be empty. In that case there is no way to get a median distance between s and s2, so it will be (stay) NA
        
        PD_l_2[[names(PD_l_2)[s]]][[names(PD_l_2)[s2]]] <- NA
        PD_l_2[[names(PD_l_2)[s2]]][[names(PD_l_2)[s]]] <- NA
        
      } else {
        
        dist_v <- c()
        for (k in 1:length(cogeneric_id3)) {
          dist_v <- c(dist_v, PD_l[[cogeneric_id3[k]]][[names(PD_l)[s2]]])
        }
        new_dist <- median(dist_v, na.rm=T)
      
        if (!is.na(new_dist)) {
      
          PD_l_2[[names(PD_l_2)[s]]][[names(PD_l_2)[s2]]] <- new_dist
          PD_l_2[[names(PD_l_2)[s2]]][[names(PD_l_2)[s]]] <- new_dist
      
        } else { # it can happen that the new distance is NA if the s2 was a species with NAs too (so any distance from a cogeneric of s to s2 was NA). In that case, we try to get the median distance between the congeneric species of s and the congeneric species of s2
        
          print(paste(names(PD_l)[s], names(PD_l)[s2], "TWO NAs?"))
          genus2 <- paste("^", str_split(names(PD_l)[s2], " ")[[1]][1], " ", sep = "")
          cogeneric_id2 <- grep(genus2, names(PD_l))
          cogeneric_id2 <- cogeneric_id2[which(!cogeneric_id2 ==s)] # removing the species under focus from the congeneric vector, otherwise will have some 0 distance biasing median if s and s2 are in same genus
          cogeneric_id2 <- cogeneric_id2[which(!cogeneric_id2 ==s2)]
          
          if (length(cogeneric_id2) < 1) { # if s2 does not have congeneric species, then cogeneric_id2 will be empty. In that case there is no way to get a median distance between s and s2, so it will be (stay) NA
            
            print(paste(names(PD_l)[s], names(PD_l)[s2], "NO CONGENERIC-ONE"))
            PD_l_2[[names(PD_l_2)[s]]][[names(PD_l_2)[s2]]] <- NA
            PD_l_2[[names(PD_l_2)[s2]]][[names(PD_l_2)[s]]] <- NA
            
          } else {
          
            dist_v2 <- c()
            for (k in 1:length(cogeneric_id3)) {
              for (z in 1:length(cogeneric_id2)) {
                dist_v2 <- c(dist_v2, PD_l[[cogeneric_id3[k]]][[names(PD_l)[cogeneric_id2[z]]]])
              }
            }
            new_dist2 <- median(dist_v2, na.rm = T)
        
            if (!is.na(new_dist2)) { # will see new_dist2 as NA if dist_v2 was only made of NAs, even if using na.rm =T, so will allow to detect cases where no distance can be calculated
          
              PD_l_2[[names(PD_l_2)[s]]][[names(PD_l_2)[s2]]] <- new_dist2
              PD_l_2[[names(PD_l_2)[s2]]][[names(PD_l_2)[s]]] <- new_dist2
          
            } else { # it can happen that the new distance is still NA if there are no congeneric species for both species # should be covered above though
          
              print(paste(names(PD_l)[s], names(PD_l)[s2], "NO CONGENERIC_TWO"))
            }
          }
        }
        
      }
      
    }
  }
}



# export the list as a dataframe
PD_m_2 <- matrix(rep(NA,length(names(PD_l_2))*length(names(PD_l_2))), ncol=length(names(PD_l_2)))
PD_df_2 <- as.data.frame(PD_m_2, row.names = names(PD_l_2))
names(PD_df_2) <- keys(PD_l_2[[1]])

for (x in 1:length(PD_l_2)) { # super inefficient
  print(x)
  sp <- rownames(PD_df_2)[which(rownames(PD_df_2) == names(PD_l_2)[x])]
  if (!is.na(PD_l_2[[x]])) {
    for (y in 1:length(names(PD_df_2))) {
      val <- PD_l_2[[x]][[names(PD_df_2)[y]]]
      PD_df_2[sp,y] <- val
    }
  }
}

write.table(PD_df_2, "PhyloDist_Average100Trees_AccNames_Scaled_2.txt", sep="\t", row.names = T)











# Calculate the first component of the score (here called score 1, PU in the paper)
# The score is the mean of ratios distance/median, on previously scaled matrices
# We add a little trick so that all bad scores (i.e. scores < 0 become 0, this allows later to multiply score 1 and score 2 and make sure two negative (bad) scores do not end up giving a positive final score)
# The trick is s1 = (s1 + |s1|) / 2, so when s1 is negative, s1 becomes 0, while when s1 is positive it remains s1
# So score 1 = ((1- (((PDu/pd_med)+(FDu/fd_med))/2)) + abs((1- (((PDu/pd_med)+(FDu/fd_med))/2))))/2


# Obtain PDs (scaled and logtransformed, simulations show that logtr does not make much difference but just in case)

pd2 <- read.table("PhyloDist_Average100Trees_AccNames_Scaled_2.txt", sep="\t", header=T, row.names = 1)
names(pd2) <- gsub("\\.", " ", names(pd2))
names(pd2)[grep("-", rownames(pd2))] <- c("Bactris ana-juliae","Calamus novae-georgii","Calamus toli-toliensis",
                                          "Chamaedorea castillo-montii","Chamaedorea ernesti-augusti",
                                          "Gaussia gomez-pompae","Raphia palma-pinus","Wallaceodoxa raja-ampat")
# (the hyphen was removed which creates problems below)

# make a dictionary of distances for each species
PD_l <- vector(mode = "list", length = length(rownames(pd2)))
names(PD_l) <- rownames(pd2)

for (s in 1:length(rownames(pd2))) {
  print(s)
  h <- hash(keys=colnames(pd2), values=pd2[s,])
  PD_l[[s]] <- h
}



# scale and logtr FDs
fd1 <- read.table("FunDist_3traits.txt", header=T, stringsAsFactors = F, sep="\t")

fd2 <- fd1 + 1
fd3 <- log(fd2)
fd <- scale(fd3, center = FALSE, scale = rep(max(fd3), length(names(fd3))))
fd <- as.data.frame(fd)


fd2 <- fd
names(fd2) <- gsub("\\.", " ", names(fd2))
names(fd2)[grep("-", rownames(fd2))] <- c("Bactris ana-juliae","Calamus novae-georgii","Calamus toli-toliensis",
                                          "Chamaedorea castillo-montii","Chamaedorea ernesti-augusti",
                                          "Gaussia gomez-pompae","Raphia palma-pinus","Wallaceodoxa raja-ampat")

# make a dictionary of distances for each species
FD_l <- vector(mode = "list", length = length(rownames(fd2)))
names(FD_l) <- rownames(fd2)

for (s in 1:length(rownames(fd2))) {
  print(s)
  h <- hash(keys=names(fd2), values=fd2[s,])
  FD_l[[s]] <- h
}

PD_l2 <- PD_l
FD_l2 <- FD_l

# 
# Make vectors of the distance matrices (only lower triangle)
fd2_v <- c()
z <- 0
for (c in 1:length(names(fd2))) {
  z <- z + 1
  fd2_v <- c(fd2_v, fd2[z:length(fd2[,1]),c])
}

pd2_v <- c()
z <- 0
for (c in 1:length(names(pd2))) {
  z <- z + 1
  pd2_v <- c(pd2_v, pd2[z:length(pd2[,1]),c])
}

# check distribution of the distances
pd2_df <- as.data.frame(pd2_v)
fd2_df <- as.data.frame(fd2_v)


hist_PD <- ggplot(pd2_df, aes(x=pd2_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PD", x="PD", y = "Count")

hist_FD <- ggplot(fd2_df, aes(x=fd2_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FD", x="FD", y = "Count")


pdf("PD_FD_ScaledLogTr_Distribution.pdf", 8, 10) 
grid.arrange(hist_PD, hist_FD, nrow = 2)
dev.off()








# Do the same but only scale the distances, do not logtransform them

# Obtain PDs (scaled and logtransformed, simulations show that logtr does not make much difference but just in case)

pd3 <- read.table("PhyloDist_Average100Trees_AccNames_Scaled_2.txt", sep="\t", header=T, row.names = 1)
names(pd3) <- gsub("\\.", " ", names(pd3))
names(pd3)[grep("-", rownames(pd3))] <- c("Bactris ana-juliae","Calamus novae-georgii","Calamus toli-toliensis",
                                          "Chamaedorea castillo-montii","Chamaedorea ernesti-augusti",
                                          "Gaussia gomez-pompae","Raphia palma-pinus","Wallaceodoxa raja-ampat")
# (the hyphen was removed which creates problems below)


# "detransform" pd3
pd3_revsc <- scale(pd3, center = FALSE, scale = rep(0.1867555, length(names(pd3)))) # the 0.1867555 value equals 1/5.354596, which is max(pd3) obtained from above when we first scaled
pd3_revsc2 <- exp(pd3_revsc)
pd3 <- pd3_revsc2 - 1

# scale it again, without logtransform
pd3 <- scale(pd3, center = FALSE, scale = rep(max(pd3, na.rm = T), length(colnames(pd3))))


# make a dictionary of distances for each species
PD_l3 <- vector(mode = "list", length = length(rownames(pd3)))
names(PD_l3) <- rownames(pd3)

for (s in 1:length(rownames(pd3))) {
  print(s)
  h <- hash(keys=colnames(pd3), values=pd3[s,])
  PD_l3[[s]] <- h
}



# scale FDs
fd1 <- read.table("FunDist_3traits.txt", header=T, stringsAsFactors = F, sep="\t")

fd3 <- scale(fd1, center = FALSE, scale = rep(max(fd1), length(names(fd1))))
fd3 <- as.data.frame(fd3)

names(fd3) <- gsub("\\.", " ", names(fd3))
names(fd3)[grep("-", rownames(fd3))] <- c("Bactris ana-juliae","Calamus novae-georgii","Calamus toli-toliensis",
                                          "Chamaedorea castillo-montii","Chamaedorea ernesti-augusti",
                                          "Gaussia gomez-pompae","Raphia palma-pinus","Wallaceodoxa raja-ampat")

# make a dictionary of distances for each species
FD_l3 <- vector(mode = "list", length = length(rownames(fd3)))
names(FD_l3) <- rownames(fd3)

for (s in 1:length(rownames(fd3))) {
  print(s)
  h <- hash(keys=names(fd3), values=fd3[s,])
  FD_l3[[s]] <- h
}


# 
# Make vectors of the distance matrices (only lower triangle)
fd3_v <- c()
z <- 0
for (c in 1:length(names(fd3))) {
  z <- z + 1
  fd3_v <- c(fd3_v, fd3[z:length(fd3[,1]),c])
}

pd3_v <- c()
z <- 0
for (c in 1:length(colnames(pd3))) {
  z <- z + 1
  pd3_v <- c(pd3_v, pd3[z:length(pd3[,1]),c])
}

# check distribution of the distances
pd3_df <- as.data.frame(pd3_v)
fd3_df <- as.data.frame(fd3_v)


hist_PD <- ggplot(pd3_df, aes(x=pd3_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PD", x="PD", y = "Count")

hist_FD <- ggplot(fd3_df, aes(x=fd3_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FD", x="FD", y = "Count")


pdf("PD_FD_ScaledOnly_Distribution.pdf", 8, 10) 
grid.arrange(hist_PD, hist_FD, nrow = 2)
dev.off()








# THEN CALCULATE SCORE
# score1 is PU from the paper, scoreTwo is PSU from the paper, and score2/PSA is SP from the paper
# Do it with scaled only (pd3, fd3, PD_l3, FD_l3) and scaled and logtransformed (pd2,fd2, PD_l2, FD_l2) distances


PD_l <- PD_l2 # change for log/not log !!!!!!!!!
FD_l <- FD_l2

# calculate medians
pd_median <- median(as.matrix(pd2), na.rm=T)   # change for log/not log !!!!!!!!!
fd_median <- median(as.matrix(fd2))

# > median(as.matrix(pd2), na.rm=T)
# [1] 0.9644245
# > median(as.matrix(fd2))
# [1] 0.06758556

# > median(as.matrix(pd3), na.rm=T)
# [1] 0.8257284
# > median(as.matrix(fd3))
# [1] 0.0440437

# Make use dictionaries
use_h_HUF <- hash(keys=d2$species, values=d2$Human_use_food)
use_h_HUMa <- hash(keys=d2$species, values=d2$Human_use_materials)
use_h_HUMe <- hash(keys=d2$species, values=d2$Human_use_medicine)
use_h_HUC <- hash(keys=d2$species, values=d2$Human_use_cultural)

# Make status dictionaries
status_h_BS <- hash(keys=d2$species, values=d2$RLr_MLpr_BS)
status_h_WS <- hash(keys=d2$species, values=d2$RLr_MLpr_WS)


# Make a table where we will fill the scores for each use (it will be NA if the species is threatened)

score_tab <- as.data.frame(d2$species)
names(score_tab)[1] <- "Species"
score_tab$s1_HUF_BS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUF_BS_pd <- rep(NA, length(score_tab$Species)) # we also fill the distances, in case we want to plot species according to them
score_tab$s1_HUF_BS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMa_BS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMa_BS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMa_BS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMe_BS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMe_BS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMe_BS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUC_BS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUC_BS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUC_BS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUF_WS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUF_WS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUF_WS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMa_WS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMa_WS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMa_WS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMe_WS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMe_WS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUMe_WS_fd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUC_WS <- rep(NA, length(score_tab$Species))
score_tab$s1_HUC_WS_pd <- rep(NA, length(score_tab$Species))
score_tab$s1_HUC_WS_fd <- rep(NA, length(score_tab$Species))


# Calculate scores (score 1: Potential for Use)
# For the moment, if min_pd is NA (because the species is not in the tree and not used itself), 
# the score remains NA
# only a few species (see output of below)

for (s in 1:length(PD_l)) {
  if (!is.na(status_h_BS[[names(PD_l)[s]]])) {
    if (status_h_BS[[names(PD_l)[s]]] == "LC") { # only calculate a score if non threatened, else score will be NA (since will not use it below anyway)
      
      h <- PD_l[[s]]
      h2 <- FD_l[[s]]
      h_df <- as.data.frame(names(PD_l))
      names(h_df)[1] <- "Species"
      h_df$PD <- rep(NA, length(h_df$Species))
      h_df$FD <- rep(NA, length(h_df$Species))
      h_df$HUF <- rep(NA, length(h_df$Species))
      h_df$HUMa <- rep(NA, length(h_df$Species))
      h_df$HUMe <- rep(NA, length(h_df$Species))
      h_df$HUC <- rep(NA, length(h_df$Species))
      for (k in 1:length(h_df$Species)) {
        h_df$PD[k] <- h[[h_df$Species[k]]]
        h_df$FD[k] <- h2[[h_df$Species[k]]]
        h_df$HUF[k] <- use_h_HUF[[h_df$Species[k]]]
        h_df$HUMa[k] <- use_h_HUMa[[h_df$Species[k]]]
        h_df$HUMe[k] <- use_h_HUMe[[h_df$Species[k]]]
        h_df$HUC[k] <- use_h_HUC[[h_df$Species[k]]]
      }
      
      #for HUF, get PD with closest used species (if that species is itself, will be 0 so will give a maximal score of 1)
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUF == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUF_BS"))} 
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUF == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUF_BS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUF_BS_pd[s] <- min_pd 
      score_tab$s1_HUF_BS_fd[s] <- min_fd
      
      # same for HUMa
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUMa == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUMa_BS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUMa == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUMa_BS[s] <-((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUMa_BS_pd[s] <- min_pd
      score_tab$s1_HUMa_BS_fd[s] <- min_fd
      
      #same for HUMe
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUMe == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUMe_BS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUMe == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUMe_BS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUMe_BS_pd[s] <- min_pd
      score_tab$s1_HUMe_BS_fd[s] <- min_fd
      
      #same for HUC
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUC == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUC_BS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUC == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUC_BS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUC_BS_pd[s] <- min_pd
      score_tab$s1_HUC_BS_fd[s] <- min_fd
      
    } 
  }
  
  if (!is.na(status_h_WS[[names(PD_l)[s]]])) {
    if (status_h_WS[[names(PD_l)[s]]] == "LC") { # only calculate a score if non threatened according to WS, else score will be NA
      
      h <- PD_l[[s]]
      h2 <- FD_l[[s]]
      h_df <- as.data.frame(names(PD_l))
      names(h_df)[1] <- "Species"
      h_df$PD <- rep(NA, length(h_df$Species))
      h_df$FD <- rep(NA, length(h_df$Species))
      h_df$HUF <- rep(NA, length(h_df$Species))
      h_df$HUMa <- rep(NA, length(h_df$Species))
      h_df$HUMe <- rep(NA, length(h_df$Species))
      h_df$HUC <- rep(NA, length(h_df$Species))
      for (k in 1:length(h_df$Species)) {
        h_df$PD[k] <- h[[h_df$Species[k]]]
        h_df$FD[k] <- h2[[h_df$Species[k]]]
        h_df$HUF[k] <- use_h_HUF[[h_df$Species[k]]]
        h_df$HUMa[k] <- use_h_HUMa[[h_df$Species[k]]]
        h_df$HUMe[k] <- use_h_HUMe[[h_df$Species[k]]]
        h_df$HUC[k] <- use_h_HUC[[h_df$Species[k]]]
      }
      
      #get PD with closest used species
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUF == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUF_WS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUF == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUF_WS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUF_WS_pd[s] <- min_pd
      score_tab$s1_HUF_WS_fd[s] <- min_fd
      
      #get PD with closest used species
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUMa == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUMa_WS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUMa == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUMa_WS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUMa_WS_pd[s] <- min_pd
      score_tab$s1_HUMa_WS_fd[s] <- min_fd
      
      #get PD with closest used species
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUMe == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUMe_WS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUMe == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUMe_WS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUMe_WS_pd[s] <- min_pd
      score_tab$s1_HUMe_WS_fd[s] <- min_fd
      
      #get PD with closest used species
      h_df2 <- h_df[order(h_df$PD),]
      h_df3 <- subset(h_df2, h_df2$HUC == "YES")
      min_pd <- h_df3$PD[1]
      min_pd_s <- min_pd/pd_median
      if (is.na(min_pd)) {print(paste(names(PD_l)[s], "HUC_WS"))}
      #get FD with closest used species
      h_df2 <- h_df[order(h_df$FD),]
      h_df3 <- subset(h_df2, h_df2$HUC == "YES")
      min_fd <- h_df3$FD[1]
      min_fd_s <- min_fd/fd_median
      score_tab$s1_HUC_WS[s] <- ((1-((min_pd_s+min_fd_s)/2)) + abs((1-((min_pd_s+min_fd_s)/2))))/2
      score_tab$s1_HUC_WS_pd[s] <- min_pd
      score_tab$s1_HUC_WS_fd[s] <- min_fd
      
    }
  }
}

write.table(score_tab, "Score1_table_final_ScaledLogTr.txt", sep="\t", row.names = F)    # change for log/not log !!!!!!!!!



# plot the distribution of score 1 depending on pd and fd
s1 <- read.table("Score1_table_final_ScaledLogTr.txt", sep="\t", header=T)      # change for log/not log !!!!!!!!!
names(s1)
s1_HUF_BS <- s1[,c(2,3,4)]
s1_HUMa_BS <- s1[,c(5,6,7)]
s1_HUMe_BS <- s1[,c(8,9,10)]
s1_HUC_BS <- s1[,c(11,12,13)]
names(s1_HUF_BS) <- c("Score", "PD", "FD")
names(s1_HUMa_BS) <- c("Score", "PD", "FD")
names(s1_HUMe_BS) <- c("Score", "PD", "FD")
names(s1_HUC_BS) <- c("Score", "PD", "FD")
s4_b <- rbind(s1_HUF_BS, s1_HUMa_BS, s1_HUMe_BS, s1_HUC_BS)

s4_plot <- ggplot(s4_b, aes(x=FD, y=PD, color=Score)) +
  geom_point() +
  labs(title="Score 1",
       x="Functional distance to closest utilized species", 
       y = "Phylogenetic distance to closest utilized species")

s_hist <- ggplot(s4_b, aes(x=Score)) + 
geom_histogram(binwidth=0.01)

pdf("PD_FD_score1_final_ScaledLogTr.pdf", 8, 10)        # change for log/not log !!!!!!!!!
grid.arrange(s4_plot, s_hist, nrow = 2)
dev.off()








# dictionary of score1 for each use (needed below)
# use BS scores, as same as WS but more species are scored
score1_tab <- read.table("Score1_table_final_ScaledLogTr.txt", header=T, sep="\t")       # change log/not log!!!!!!!!
score1_h_l <- vector(mode="list", length=4)
names(score1_h_l) <- c("HUF", "HUMa", "HUMe", "HUC")

score1_h_l[["HUF"]] <- hash(keys=score1_tab$Species, values=score1_tab$s1_HUF_BS)
score1_h_l[["HUMa"]] <- hash(keys=score1_tab$Species, values=score1_tab$s1_HUMa_BS)
score1_h_l[["HUMe"]] <- hash(keys=score1_tab$Species, values=score1_tab$s1_HUMe_BS)
score1_h_l[["HUC"]] <- hash(keys=score1_tab$Species, values=score1_tab$s1_HUC_BS)



# dictionary of species for each region (will be needed below)
# warnings ok (happen when check if region has already species or not: when it has, it is a vector, hence the warning, but fine because NA only indicates truly empty value, as vector of species can never start with NA - checked with grep("TRUE", is.na(df12$species)), which gives integer(0))
regions_sp_h <- hash(keys=tdwg3$LEVEL3_COD, values=rep(NA, length(tdwg3$LEVEL3_COD)))
for (x in 1:length(df12$species)) {
  if (is.na(regions_sp_h[[df12$area_code_l3[x]]])) {
    regions_sp_h[[df12$area_code_l3[x]]] <- c(df12$species[x])
  } else {
    regions_sp_h[[df12$area_code_l3[x]]] <- c(regions_sp_h[[df12$area_code_l3[x]]], df12$species[x])
  }
}

# Make status dictionary (will need it too)
status_d <- read.table("Species_status_dictionary.txt", sep="\t", header=T, stringsAsFactors = T)
status_d$Status_BSWS2 <- mapvalues(status_d$Status_BSWS, from=levels(status_d$Status_BSWS), to=c("None", "WSonly", "Unknown", "BSandWS"))
status_d_h <- hash(keys=status_d$Species, values=as.character(status_d$Status_BSWS2))

# Calculate score 2 of co-occuring non-threatened species for each UT species, for each region, for each use, for each scenario
# What we call score2 below is what we call later PSA and in the paper SP, i.e. the Substitution Potential. The intermediate Potential for Similar Use is called scoreTwo below.
# BS: only consider UT species under BS and only calculate score 2 for co-occuring species that are NT according to BS (using s1 BS from above)
# WS: only consider UT species under WS and only calculate score 2 for co-occuring species that are NT according to WS (using s1 WS from above)

# Make an empty vector to store the scores to be able to study their distribution
scores_v <- c()


# Make a list of results tables per region with as many tables as regions and names of tables being region names
reg_l <- vector(mode = "list", length = length(keys(regions_sp_h)))
names(reg_l) <- keys(regions_sp_h)

uses <- c("HUF", "HUMa", "HUMe", "HUC") # should be in the same order as names(score1_h_l) !!!

for (r in 1:length(names(reg_l))) {
 if (!is.na(regions_sp_h[[names(reg_l)[r]]][1])) { # if there are palms in the region (this will generate a warning as it is a vector of the species in the region, fine)
 tab <- as.data.frame(regions_sp_h[[names(reg_l)[r]]])
 tab[,2:102] <- as.data.frame(matrix(data=NA, nrow = length(tab[,1]), ncol = 101))
 names(tab) <- c("Species", "Status", uses, 
                 paste("HUF", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep=""),
                 paste("HUMa", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep=""),
                 paste("HUMe", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep=""),
                 paste("HUC", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep="")
                 )

 
 for (s in 1:length(tab$Species)) {
   tab$Status[s] <- status_d_h[[tab$Species[s]]]
   tab$HUF[s] <- use_h_HUF[[tab$Species[s]]]
   tab$HUMa[s] <- use_h_HUMa[[tab$Species[s]]]
   tab$HUMe[s] <- use_h_HUMe[[tab$Species[s]]]
   tab$HUC[s] <- use_h_HUC[[tab$Species[s]]]
   u2 <- 0
   
   if (tab$Status[s] == "WSonly" || tab$Status[s] == "BSandWS") { # if the species is threatened (in any scenario)
   for (u in 1:length(uses)) {
     if (!is.na(tab[s,2+u])) { # if something is known about the species use  # The 2 is hard-coded, change if more than two columns before the uses columns in tab
     if (tab[s,2+u] == "YES") { # if the species is used for the considered use
       
       # build vector of alternative species with scores, for each scenario
       
       species_C1_BS <- c()
       species_C2_BS <- c()
       species_C3_BS <- c()
       species_C4_BS <- c()
       species_C1_WS <- c()
       species_C2_WS <- c()
       species_C3_WS <- c()
       species_C4_WS <- c()
       
       for (s2 in 1:length(tab$Species)) { # for each co-occurring species
         
         if (!tab$Species[s2] == tab$Species[s]) { # excluding the species to replace
           
            
              if (status_d_h[[tab$Species[s2]]] == "None" || status_d_h[[tab$Species[s2]]] == "WSonly") { # we want species that are known to be not threatened under the BS
                
                score1 <- score1_h_l[[uses[u]]][[tab$Species[s2]]] # score calculation is repeated here and below because score1 is NA for threatened species so would produce errors if did it for any species before entering the loop
                PD <- PD_l[[tab$Species[s]]][[tab$Species[s2]]]
                FD <- FD_l[[tab$Species[s]]][[tab$Species[s2]]]
                if (!is.na(PD) && !is.na(FD) && !is.na(score1)) {
                  scoreTwo <- ( (  (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) + abs( (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) )   )/2)
                  score2 <- score1 * scoreTwo
                } else {
                  score2 <- NA
                }
                score_all <- paste(PD, FD, score1, scoreTwo, score2, sep="\t")
                scores_v <- c(scores_v, score_all)
                
                if (!is.na(score2)) { # if it is NA, the species will not count as alternative, it can happen if there is no congeneric species: then cannot calculate median distance (see PD calculations above)
                if (score2 <= 0.4) {
                species_C1_BS <- c(species_C1_BS, paste(tab$Species[s2], round(score2, 3), sep="___")) # input species in the right alternatives vector depending on score2 category
                } else {
                  if (score2 <=0.6) {
                    species_C2_BS <- c(species_C2_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                  } else {
                    if (score2 <=0.8) {
                      species_C3_BS <- c(species_C3_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                    } else {
                      species_C4_BS <- c(species_C4_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                    }
                  }
                }
               } else {print(paste(names(reg_l)[r], tab$Species[s], tab$Species[s2], "BS", "NO CONGENERIC!!!!!!!!!!!!!!!!!!!!!!!!!"))} # just to check how many are missed
              }

 
                            
              if (status_d_h[[tab$Species[s2]]] == "None") { # we want species that are known to be not threatened even under the BS
                
                score1 <- score1_h_l[[uses[u]]][[tab$Species[s2]]]
                PD <- PD_l[[tab$Species[s]]][[tab$Species[s2]]]
                FD <- FD_l[[tab$Species[s]]][[tab$Species[s2]]]
                if (!is.na(PD) && !is.na(FD) && !is.na(score1)) {
                  scoreTwo <- ( (  (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) + abs( (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) )   )/2)
                  score2 <- score1 * scoreTwo
                } else {
                  score2 <- NA
                }

                
                if (!is.na(score2)) { # if it is NA, the species will not count as alternative, it can happen if there is no congeneric species: then cannot calculate median distance (see PD calculations above)
                if (score2 <= 0.4) {
                  species_C1_WS <- c(species_C1_WS, paste(tab$Species[s2], round(score2, 3), sep="___")) # input species in the right alternatives vector depending on score2 category
                } else {
                  if (score2 <=0.6) {
                    species_C2_WS <- c(species_C2_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                  } else {
                    if (score2 <=0.8) {
                      species_C3_WS <- c(species_C3_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                    } else {
                      species_C4_WS <- c(species_C4_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                    }
                  }
                }
               } else {print(paste(names(reg_l)[r], tab$Species[s], tab$Species[s2], "WS", "NO CONGENERIC!!!!!!!!!!!!!!!!!!!!!!!!!"))} # just to check how many are missed
              }
    
          }      
       }
       
      # fill the table
       vec_l <- vector(mode="list", length=8)
       if (is.null(species_C1_BS)) {vec_l[[1]] <- NA} else {vec_l[[1]] <- species_C1_BS} # make a list of vectors for each alternative category, with NA if no alternatives
       if (is.null(species_C2_BS)) {vec_l[[2]] <- NA} else {vec_l[[2]] <- species_C2_BS}
       if (is.null(species_C3_BS)) {vec_l[[3]] <- NA} else {vec_l[[3]] <- species_C3_BS}
       if (is.null(species_C4_BS)) {vec_l[[4]] <- NA} else {vec_l[[4]] <- species_C4_BS}
       if (is.null(species_C1_WS)) {vec_l[[5]] <- NA} else {vec_l[[5]] <- species_C1_WS}
       if (is.null(species_C2_WS)) {vec_l[[6]] <- NA} else {vec_l[[6]] <- species_C2_WS}
       if (is.null(species_C3_WS)) {vec_l[[7]] <- NA} else {vec_l[[7]] <- species_C3_WS}
       if (is.null(species_C4_WS)) {vec_l[[8]] <- NA} else {vec_l[[8]] <- species_C4_WS}
       
       for (v in 1:length(vec_l)) { # transform vectors into strings and fill them in tab # 6, 8 and 16 are hard-coded, change if changing the tab structure
         if (!is.na(vec_l[[v]])) {
          tab[s,(6+u2+v)] <- toString(vec_l[[v]]) # alternatives themselves
          tab[s,(6+u2+8+v)] <- length(vec_l[[v]]) # number of alternatives for that category
          if (length(vec_l[[v]]) > 0) { # binary coding of presence (1) / absence (0) of alternatives for that category
            tab[s,(6+u2+16+v)] <- 1
          } else {tab[s,(6+u2+16+v)] <- 0}
         } else {
          tab[s,(6+u2+v)] <- NA
          tab[s,(6+u2+8+v)] <- 0
          tab[s,(6+u2+16+v)] <- 0
         }
       }
       
       
     } else { # if the species is not used, put NA for alternatives
       tab[s,(6+u2+1):(6+u2+24)] <- rep(NA, 24)      # The numbers 6, 24, and 96 or 7:102 below are hard coded, change if change structure of tab
     }
       
     } else { # if nothing is known about the species use, put NA for alternatives
       tab[s,(6+u2+1):(6+u2+24)] <- rep(NA, 24)
     }
     
     u2 <- u2+24
   }
   
   } else { # if the species is not threatened, put NA for alternatives
     tab[s,7:102] <- rep(NA, 96)
   }
   
 }
 
 reg_l[[r]] <- tab
 
 } else {
   reg_l[[r]] <- NA # region table will be NA if no palms in the region
 }
 
 }


# export vector of scores as a table
scores_df <- as.data.frame(scores_v)
names(scores_df)[1] <- "PD\tFD\tScore1\tScore2\tPSA" 
# In the table/header, Score1 is PU from the paper (score1 above), score2 is PSU from the paper (scoreTwo above), and PSA is SP from the paper (score2 above)
# PD and FD are the distances used to calculate PSU
write.table(scores_df, "All_scores_final_ScaledLogTr.txt", row.names = F)        # change for log/not log !!!!!!!!!
# Remove quotes manually in Notepad, then can reimport and plot anytime !!!!!!!!!!!!!!!!!!!!

# export region tables

final_table <- reg_l[[2]][1,] # use first region with a table to create a template, CHANGE AS NEEDED!!!
final_table$Region <- names(reg_l)[2] # use first region with a table to create a template, CHANGE AS NEEDED!!!

for (f in 1:length(reg_l)) {
  if (!is.na(reg_l[[f]])) {
    t <- reg_l[[f]]
    t$Region <- rep(names(reg_l)[f], length(t$Species))
    final_table <- rbind(final_table, t)
  }
}

final_table <- final_table[-1,]

write.table(final_table, "Alternatives_per_region_final_ScaledLogTr.txt", sep="\t")   # change for log/not log !!!!!!!!!




# create lists of tables, with as many tables as regions and names of tables being region names
reg_l_noNA <- reg_l[!is.na(reg_l)]
regP_l <- vector(mode = "list", length = length(reg_l_noNA))
names(regP_l) <- names(reg_l_noNA)

uses <- c("HUF", "HUMa", "HUMe", "HUC")
cats <- c("C1", "C2", "C3", "C4")

final_tab <- as.data.frame(matrix(data=NA, nrow = 1, ncol = 7))
names(final_tab) <- c("Region", "Scenario", "Use", "Category", "PercentageA", "NumberA", "NumberUT")

for (r in 1:length(names(regP_l))) {
  tab <- as.data.frame(matrix(data=NA, nrow = 32, ncol = 7))
  names(tab) <- c("Region", "Scenario", "Use", "Category", "PercentageA", "NumberA", "NumberUT")
  
  region <- names(regP_l)[r]
  
  t <- reg_l_noNA[[r]]
  t_BS <- subset(t, t$Status == "BSandWS")
  t_WS <- subset(t, t$Status == "BSandWS" | t$Status == "WSonly") # WS in this case is when species to replace are threatened under BS and WS or even just WS, whereas BS is when they are threatened under BS, and therefore also WS
  numUT_BS <- length(t_BS[,1]) # not really numUT, in fact just numT
  numUT_WS <- length(t_WS[,1])
  
  numUT_BS_v <- c()
  numUT_WS_v <- c()
  
  for (u in 1:length(uses)) {
    use <- uses[u]
    numUT_BS_v[u]<- length(which(t_BS[,2+u] == "YES"))
    numUT_WS_v[u]<- length(which(t_WS[,2+u] == "YES"))
  }
  
  z <- 0
  
  for (u in 1:length(uses)) {
    use <- uses[u]
    for (c in 1:length(cats)) {
      z <- z+1
      
      category <- cats[c]
      BS_col <- paste(use,"Bin_", category, "_BS", sep="") # will depend on u and c and original table ordering
      WS_col <- paste(use,"Bin_", category, "_WS", sep="") # will depend on u and c and original table ordering
      
      numA_BS <- sum(t_BS[,which(names(t_BS) == BS_col)], na.rm = T) 
      perA_BS <- (numA_BS*100)/numUT_BS_v[u]
      numA_WS <- sum(t_WS[,which(names(t_WS) == WS_col)], na.rm = T) 
      perA_WS <- (numA_WS*100)/numUT_WS_v[u]
      
      resBS <- c(region, "BS", use, category, perA_BS, numA_BS, numUT_BS_v[u])
      resWS <- c(region, "WS", use, category, perA_WS, numA_WS, numUT_WS_v[u])
      
      tab[z,1:7] <- resBS
      tab[z+1,1:7] <- resWS
      z <- z+1
    }
  }
  
  
  final_tab <- rbind(final_tab, tab)
}

final_tab <- final_tab[-1,]




write.table(final_tab, "Alternatives_Percentages_per_region_final_ScaledLogTr.txt", sep="\t")          # change for log/not log !!!!!!!!!








# do the same thing but cumulatively, i.e. alternative species of category 4 also count as cat3, 2 and 1, alt. species of cat3 also count as 2 and 1 etc

# Calculate score 2 of co-occuring non-threatened species for each UT species, for each region, for each use, for each scenario
# BS: only consider UT species under BS and only calculate score 2 for co-occuring species that are NT according to BS (using s1 BS from above)
# WS: only consider UT species under WS and only calculate score 2 for co-occuring species that are NT according to WS (using s1 WS from above)

# Make a list of results tables per region with as many tables as regions and names of tables being region names
reg_l <- vector(mode = "list", length = length(keys(regions_sp_h)))
names(reg_l) <- keys(regions_sp_h)

uses <- c("HUF", "HUMa", "HUMe", "HUC") # should be in the same order as names(score1_h_l) !!!

for (r in 1:length(names(reg_l))) {
  if (!is.na(regions_sp_h[[names(reg_l)[r]]][1])) { # if there are palms in the region (this will generate a warning as it is a vector of the species in the region, fine)
    tab <- as.data.frame(regions_sp_h[[names(reg_l)[r]]])
    tab[,2:102] <- as.data.frame(matrix(data=NA, nrow = length(tab[,1]), ncol = 101))
    names(tab) <- c("Species", "Status", uses, 
                    paste("HUF", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                   paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                   paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep=""),
                    paste("HUMa", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                    paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                    paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep=""),
                    paste("HUMe", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                    paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                    paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep=""),
                    paste("HUC", c(paste("Spe_C", c(1:4), "_BS", sep=""), paste("Spe_C", c(1:4), "_WS", sep=""),
                                   paste("Num_C", c(1:4), "_BS", sep=""), paste("Num_C", c(1:4), "_WS", sep=""),
                                   paste("Bin_C", c(1:4), "_BS", sep=""), paste("Bin_C", c(1:4), "_WS", sep="")), sep="")
    )
    
    
    for (s in 1:length(tab$Species)) {
      tab$Status[s] <- status_d_h[[tab$Species[s]]]
      tab$HUF[s] <- use_h_HUF[[tab$Species[s]]]
      tab$HUMa[s] <- use_h_HUMa[[tab$Species[s]]]
      tab$HUMe[s] <- use_h_HUMe[[tab$Species[s]]]
      tab$HUC[s] <- use_h_HUC[[tab$Species[s]]]
      u2 <- 0
      
      if (tab$Status[s] == "WSonly" || tab$Status[s] == "BSandWS") { # if the species is threatened (in any scenario)
        for (u in 1:length(uses)) {
          if (!is.na(tab[s,2+u])) { # if something is known about the species use  # The 2 is hard-coded, change if more than two columns before the uses columns in tab
            if (tab[s,2+u] == "YES") { # if the species is used for the considered use
              
              # build vector of alternative species with scores, for each scenario
              
              species_C1_BS <- c()
              species_C2_BS <- c()
              species_C3_BS <- c()
              species_C4_BS <- c()
              species_C1_WS <- c()
              species_C2_WS <- c()
              species_C3_WS <- c()
              species_C4_WS <- c()
              
              for (s2 in 1:length(tab$Species)) { # for each co-occurring species
                
                if (!tab$Species[s2] == tab$Species[s]) { # excluding the species to replace
                  
                  
                  if (status_d_h[[tab$Species[s2]]] == "None" || status_d_h[[tab$Species[s2]]] == "WSonly") { # we want species that are known to be not threatened under the BS
                    
                    score1 <- score1_h_l[[uses[u]]][[tab$Species[s2]]] # score calculation is repeated here and below because score1 is NA for threatened species so would produce errors if did it for any species before entering the loop
                    PD <- PD_l[[tab$Species[s]]][[tab$Species[s2]]]
                    FD <- FD_l[[tab$Species[s]]][[tab$Species[s2]]]
                    if (!is.na(PD) && !is.na(FD) && !is.na(score1)) {
                      score2 <- score1 * ( (  (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) + abs( (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) )   )/2)
                    } else {
                      score2 <- NA
                    }
                    
                    if (!is.na(score2)) { # if it is NA, the species will not count as alternative, it can happen if there is no congeneric species: then cannot calculate median distance (see PD calculations above)
                      if (score2 <= 0.4) {
                        species_C1_BS <- c(species_C1_BS, paste(tab$Species[s2], round(score2, 3), sep="___")) # input species in the right alternatives vector depending on score2 category
                      } else {
                        if (score2 <=0.6) {
                          species_C2_BS <- c(species_C2_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                          species_C1_BS <- c(species_C1_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                        } else {
                          if (score2 <=0.8) {
                            species_C3_BS <- c(species_C3_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C1_BS <- c(species_C1_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C2_BS <- c(species_C2_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                          } else {
                            species_C4_BS <- c(species_C4_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C1_BS <- c(species_C1_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C2_BS <- c(species_C2_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C3_BS <- c(species_C3_BS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                          }
                        }
                      }
                    } else {print(paste(names(reg_l)[r], tab$Species[s], tab$Species[s2], "BS", "NO CONGENERIC!!!!!!!!!!!!!!!!!!!!!!!!!"))} # just to check how many are missed
                  }
                  
                  
                  
                  if (status_d_h[[tab$Species[s2]]] == "None") { # we want species that are known to be not threatened even under the BS
                    
                    score1 <- score1_h_l[[uses[u]]][[tab$Species[s2]]]
                    PD <- PD_l[[tab$Species[s]]][[tab$Species[s2]]]
                    FD <- FD_l[[tab$Species[s]]][[tab$Species[s2]]]
                    if (!is.na(PD) && !is.na(FD) && !is.na(score1)) {
                      score2 <- score1 * ( (  (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) + abs( (1-( ( (PD/pd_median) + (FD/fd_median) ) /2)) )   )/2)
                    } else {
                      score2 <- NA
                    }
                    
                    if (!is.na(score2)) { # if it is NA, the species will not count as alternative, it can happen if there is no congeneric species: then cannot calculate median distance (see PD calculations above)
                      if (score2 <= 0.4) {
                        species_C1_WS <- c(species_C1_WS, paste(tab$Species[s2], round(score2, 3), sep="___")) # input species in the right alternatives vector depending on score2 category
                      } else {
                        if (score2 <=0.6) {
                          species_C2_WS <- c(species_C2_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                          species_C1_WS <- c(species_C1_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                        } else {
                          if (score2 <=0.8) {
                            species_C3_WS <- c(species_C3_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C2_WS <- c(species_C2_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C1_WS <- c(species_C1_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                          } else {
                            species_C4_WS <- c(species_C4_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C3_WS <- c(species_C3_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C2_WS <- c(species_C2_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                            species_C1_WS <- c(species_C1_WS, paste(tab$Species[s2], round(score2, 3), sep="___"))
                          }
                        }
                      }
                    } else {print(paste(names(reg_l)[r], tab$Species[s], tab$Species[s2], "WS", "NO CONGENERIC!!!!!!!!!!!!!!!!!!!!!!!!!"))} # just to check how many are missed
                  }
                  
                }      
              }
              
              # fill the table
              vec_l <- vector(mode="list", length=8)
              if (is.null(species_C1_BS)) {vec_l[[1]] <- NA} else {vec_l[[1]] <- species_C1_BS} # make a list of vectors for each alternative category, with NA if no alternatives
              if (is.null(species_C2_BS)) {vec_l[[2]] <- NA} else {vec_l[[2]] <- species_C2_BS}
              if (is.null(species_C3_BS)) {vec_l[[3]] <- NA} else {vec_l[[3]] <- species_C3_BS}
              if (is.null(species_C4_BS)) {vec_l[[4]] <- NA} else {vec_l[[4]] <- species_C4_BS}
              if (is.null(species_C1_WS)) {vec_l[[5]] <- NA} else {vec_l[[5]] <- species_C1_WS}
              if (is.null(species_C2_WS)) {vec_l[[6]] <- NA} else {vec_l[[6]] <- species_C2_WS}
              if (is.null(species_C3_WS)) {vec_l[[7]] <- NA} else {vec_l[[7]] <- species_C3_WS}
              if (is.null(species_C4_WS)) {vec_l[[8]] <- NA} else {vec_l[[8]] <- species_C4_WS}
              
              for (v in 1:length(vec_l)) { # transform vectors into strings and fill them in tab # 6, 8 and 16 are hard-coded, change if changing the tab structure
                if (!is.na(vec_l[[v]])) {
                  tab[s,(6+u2+v)] <- toString(vec_l[[v]]) # alternatives themselves
                  tab[s,(6+u2+8+v)] <- length(vec_l[[v]]) # number of alternatives for that category
                  if (length(vec_l[[v]]) > 0) { # binary coding of presence (1) / absence (0) of alternatives for that category
                    tab[s,(6+u2+16+v)] <- 1
                  } else {tab[s,(6+u2+16+v)] <- 0}
                } else {
                  tab[s,(6+u2+v)] <- NA
                  tab[s,(6+u2+8+v)] <- 0
                  tab[s,(6+u2+16+v)] <- 0
                }
              }
              
              
            } else { # if the species is not used, put NA for alternatives
              tab[s,(6+u2+1):(6+u2+24)] <- rep(NA, 24)      # The numbers 6, 24, and 96 or 7:102 below are hard coded, change if change structure of tab
            }
            
          } else { # if nothing is known about the species use, put NA for alternatives
            tab[s,(6+u2+1):(6+u2+24)] <- rep(NA, 24)
          }
          
          u2 <- u2+24
        }
        
      } else { # if the species is not threatened, put NA for alternatives
        tab[s,7:102] <- rep(NA, 96)
      }
      
    }
    
    reg_l[[r]] <- tab
    
  } else {
    reg_l[[r]] <- NA # region table will be NA if no palms in the region
  }
  
}





# export region tables

final_table <- reg_l[[2]][1,] # use first region with a table to create a template, CHANGE AS NEEDED!!!
final_table$Region <- names(reg_l)[2] # use first region with a table to create a template, CHANGE AS NEEDED!!!

for (f in 1:length(reg_l)) {
  if (!is.na(reg_l[[f]])) {
    t <- reg_l[[f]]
    t$Region <- rep(names(reg_l)[f], length(t$Species))
    final_table <- rbind(final_table, t)
  }
}

final_table <- final_table[-1,]

write.table(final_table, "Alternatives_per_region_cumulative_final_ScaledLogTr.txt", sep="\t")        # change for log/not log !!!!!!!!!



# for each region, for each scenario, for each use, get % of UT species with alternatives in each s2 category

# create lists of tables, with as many tables as regions and names of tables being region names
reg_l_noNA <- reg_l[!is.na(reg_l)]
regP_l <- vector(mode = "list", length = length(reg_l_noNA))
names(regP_l) <- names(reg_l_noNA)

uses <- c("HUF", "HUMa", "HUMe", "HUC")
cats <- c("C1", "C2", "C3", "C4")

final_tab <- as.data.frame(matrix(data=NA, nrow = 1, ncol = 7))
names(final_tab) <- c("Region", "Scenario", "Use", "Category", "PercentageA", "NumberA", "NumberUT")

for (r in 1:length(names(regP_l))) {
  tab <- as.data.frame(matrix(data=NA, nrow = 32, ncol = 7))
  names(tab) <- c("Region", "Scenario", "Use", "Category", "PercentageA", "NumberA", "NumberUT")

  region <- names(regP_l)[r]
  
  t <- reg_l_noNA[[r]]
  t_BS <- subset(t, t$Status == "BSandWS")
  t_WS <- subset(t, t$Status == "BSandWS" | t$Status == "WSonly") # WS in this case is when species to replace are threatened under BS and WS or even just WS, whereas BS is when they are threatened under BS, and therefore also WS
  numUT_BS <- length(t_BS[,1]) # not really numUT, in fact just numT
  numUT_WS <- length(t_WS[,1])
  
  numUT_BS_v <- c()
  numUT_WS_v <- c()
  
  for (u in 1:length(uses)) {
    use <- uses[u]
    numUT_BS_v[u]<- length(which(t_BS[,2+u] == "YES"))
    numUT_WS_v[u]<- length(which(t_WS[,2+u] == "YES"))
  }
  
  z <- 0
  
  for (u in 1:length(uses)) {
    use <- uses[u]
    for (c in 1:length(cats)) {
      z <- z+1
      
      category <- cats[c]
      BS_col <- paste(use,"Bin_", category, "_BS", sep="") # will depend on u and c and original table ordering
      WS_col <- paste(use,"Bin_", category, "_WS", sep="") # will depend on u and c and original table ordering
      
      numA_BS <- sum(t_BS[,which(names(t_BS) == BS_col)], na.rm = T) 
      perA_BS <- (numA_BS*100)/numUT_BS_v[u]
      numA_WS <- sum(t_WS[,which(names(t_WS) == WS_col)], na.rm = T) 
      perA_WS <- (numA_WS*100)/numUT_WS_v[u]
      
      resBS <- c(region, "BS", use, category, perA_BS, numA_BS, numUT_BS_v[u])
      resWS <- c(region, "WS", use, category, perA_WS, numA_WS, numUT_WS_v[u])
      
      tab[z,1:7] <- resBS
      tab[z+1,1:7] <- resWS
      z <- z+1
    }
  }
  
  
  final_tab <- rbind(final_tab, tab)
}

final_tab <- final_tab[-1,]




write.table(final_tab, "Alternatives_Percentages_per_region_cumulative_final_ScaledLogTr.txt", sep="\t")         # change for log/not log !!!!!!!!!



# NB: when "NumberUT" is 0, NumberA is 0 but this is meaningless as if there is no species to substitute we do not look how many of them have alternatives


##########################################



########################
# Examine observed PSA behavior
########################

# Remove quotes manually from All_scores_final* in Notepad, then can reimport and plot anytime !!!!!!!!!!!!!!!!!!!!
scores_all <- read.table("All_scores_final_ScaledLogTr.txt", header=T)       # change for log/not log !!!!!!!!

hist_PDt <- ggplot(scores_all, aes(x=PD)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PDt", x="PDt", y = "Count")
hist_FDt <- ggplot(scores_all, aes(x=FD)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FDt", x="FDt", y = "Count")
hist_s1 <- ggplot(scores_all, aes(x=Score1)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="Score 1", x="Score 1", y = "Count")
hist_s2 <- ggplot(scores_all, aes(x=Score2)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="Score 2", x="Score 2", y = "Count")
hist_PSA <- ggplot(scores_all, aes(x=PSA)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PSA", x="PSA", y = "Count")


pdf("PD_FD_Observed_Scores_distributions_ScaledLogTr.pdf", 10, 3)      # change for log/not log !!!!!!!!
grid.arrange(hist_PDt, hist_FDt, hist_s1, hist_s2, hist_PSA, nrow = 1)
dev.off()


# check score 2
plot_s2 <- ggplot(scores_all, aes(x=PD, y=FD, color=cut(Score2, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="Score 2")+
  theme(legend.position="bottom")+
  labs(title="Score 2",
       x="PDt", 
       y = "FDt")

plot_PSA <- ggplot(scores_all, aes(x=Score1, y=Score2, color=cut(PSA, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="PSA")+
  theme(legend.position="bottom")+
  labs(title="PSA",
       x="Score 1", 
       y = "Score 2")



pdf("PD_FD_Observed_Scores_plots_ScaledLogTr.pdf", 18, 8)          # change for log/not log !!!!!!!!
grid.arrange(plot_s2, plot_PSA, nrow = 1)
dev.off()



# Seems that scaled only makes more sense when looking at plots
# In terms of categories, it makes sense to use <0.6, [0.6-0.8], and >0.8 (but will also use <0.4 to keep 4 categories since coded like that)







########################
# Examine simulated score behavior
########################
# simulate some scores 1 and 2 and PSA depending on simulated (scaled) FDu, PDu, FD and PD


# > median(as.matrix(pd2), na.rm=T)
# [1] 0.9644245
# > median(as.matrix(fd2))
# [1] 0.06758556

# > median(as.matrix(pd3), na.rm=T)
# [1] 0.8257284
# > median(as.matrix(fd3))
# [1] 0.0440437


# First for ScaledLogTr
# min and max of fd and pd are 0 and 1
fd_med <- median(as.matrix(fd2))
pd_med <- median(as.matrix(pd2), na.rm=T)

# Make vectors of the distance matrices (only lower triangle)
fd2_v <- c()
z <- 0
for (c in 1:length(names(fd2))) {
  z <- z + 1
  fd2_v <- c(fd2_v, fd2[z:length(fd2[,1]),c])
}

pd2_v <- c()
z <- 0
for (c in 1:length(names(pd2))) {
  z <- z + 1
  pd2_v <- c(pd2_v, pd2[z:length(pd2[,1]),c])
}

# remove NAs from the PD vector (the NAs are due to Sabinaria and Wallaceodoxa - see above - and we do not want to sample them)
pd2_v <- pd2_v[!is.na(pd2_v)]

# sample randomly the vectors to obtain some PDu, FDu, PDt and FDt
PDu_v <- sample(pd2_v, 100000, replace = T)
PDt_v <- sample(pd2_v, 100000, replace = T)
FDu_v <- sample(fd2_v, 100000, replace = T)
FDt_v <- sample(fd2_v, 100000, replace = T)

# Put that in a dataframe, so that one row is a potential species to score, with a given PDu, FDu, PDt and FDt (the latter to an hypothetical species to replace)
simulated_s <- as.data.frame(PDu_v)
simulated_s$PDt_v <- PDt_v
simulated_s$FDu_v <- FDu_v
simulated_s$FDt_v <- FDt_v
simulated_s$score1 <- rep(NA, length(simulated_s$PDu_v))
simulated_s$score2 <- rep(NA, length(simulated_s$PDu_v))
simulated_s$PSA <- rep(NA, length(simulated_s$PDu_v))

# Calculate the scores for each row
for (x in 1:length(simulated_s$PDu_v)) {
  simulated_s$score1[x] <- ((1- (((simulated_s$PDu_v[x]/pd_med)+(simulated_s$FDu_v[x]/fd_med))/2)) + abs((1- (((simulated_s$PDu_v[x]/pd_med)+(simulated_s$FDu_v[x]/fd_med))/2))))/2
  simulated_s$score2[x] <- ((1- (((simulated_s$PDt_v[x]/pd_med)+(simulated_s$FDt_v[x]/fd_med))/2)) + abs((1- (((simulated_s$PDt_v[x]/pd_med)+(simulated_s$FDt_v[x]/fd_med))/2))))/2
  simulated_s$PSA[x] <- simulated_s$score1[x] * simulated_s$score2[x]
}

# Export the scores
write.table(simulated_s, "simulated_PSA_ScaledlogTr.txt", sep="\t")

max(simulated_s$score1)
# [1] 0.9859333
max(simulated_s$score2)
# [1] 0.9837785
max(simulated_s$PSA)
# [1] 0.3682433
# > 
min(simulated_s$score1)
# [1] 0
min(simulated_s$score2)
# [1] 0
min(simulated_s$PSA)
# [1] 0


# Check the distribution of all distances and scores
sim_hist_PDu <- ggplot(simulated_s, aes(x=PDu_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PDu", x="PDu", y = "Count")
sim_hist_FDu <- ggplot(simulated_s, aes(x=FDu_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FDu", x="FDu", y = "Count")
sim_hist_PDt <- ggplot(simulated_s, aes(x=PDt_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PDt", x="PDt", y = "Count")
sim_hist_FDt <- ggplot(simulated_s, aes(x=FDt_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FDt", x="FDt", y = "Count")
sim_hist_s1 <- ggplot(simulated_s, aes(x=score1)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="Score 1", x="Score 1", y = "Count")
sim_hist_s2 <- ggplot(simulated_s, aes(x=score2)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="Score 2", x="Score 2", y = "Count")
sim_hist_PSA <- ggplot(simulated_s, aes(x=PSA)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PSA", x="PSA", y = "Count")


# check score 1
sim_plot_s1 <- ggplot(simulated_s, aes(x=PDu_v, y=FDu_v, color=cut(score1, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="Score 1")+
  theme(legend.position="bottom")+
  labs(title="Score 1",
       x="PDu", 
       y = "FDu")

# check score 2
sim_plot_s2 <- ggplot(simulated_s, aes(x=PDt_v, y=FDt_v, color=cut(score2, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="Score 2")+
  theme(legend.position="bottom")+
  labs(title="Score 2",
       x="PDt", 
       y = "FDt")

sim_plot_PSA <- ggplot(simulated_s, aes(x=score1, y=score2, color=cut(PSA, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="PSA")+
  theme(legend.position="bottom")+
  labs(title="PSA",
       x="Score 1", 
       y = "Score 2")



pdf("PD_FD_Simulated_Scores_distributions_ScaledLogTr.pdf", 10, 5) 
grid.arrange(sim_hist_PDu, sim_hist_FDu, sim_hist_PDt, sim_hist_FDt, 
             sim_hist_s1, sim_hist_s2, sim_hist_PSA, nrow = 2)
dev.off()

pdf("PD_FD_Simulated_Scores_plots_ScaledLogTr.pdf", 26, 8) 
grid.arrange(sim_plot_s1, sim_plot_s2, sim_plot_PSA, nrow = 1)
dev.off()



# Then for scaled Only:

# min and max of fd and pd are 0 and 1
fd_med <- median(as.matrix(fd3))
pd_med <- median(as.matrix(pd3), na.rm=T)

# Make vectors of the distance matrices (only lower triangle)
fd3_v <- c()
z <- 0
for (c in 1:length(names(fd3))) {
  z <- z + 1
  fd3_v <- c(fd3_v, fd3[z:length(fd3[,1]),c])
}

pd3_v <- c()
z <- 0
for (c in 1:length(names(pd3))) {
  z <- z + 1
  pd3_v <- c(pd3_v, pd3[z:length(pd3[,1]),c])
}

# remove NAs from the PD vector (the NAs are due to Sabinaria and Wallaceodoxa - see above - and we do not want to sample them)
pd3_v <- pd3_v[!is.na(pd3_v)]

# sample randomly the vectors to obtain some PDu, FDu, PDt and FDt
PDu_v <- sample(pd3_v, 100000, replace = T)
PDt_v <- sample(pd3_v, 100000, replace = T)
FDu_v <- sample(fd3_v, 100000, replace = T)
FDt_v <- sample(fd3_v, 100000, replace = T)

# Put that in a dataframe, so that one row is a potential species to score, with a given PDu, FDu, PDt and FDt (the latter to an hypothetical species to replace)
simulated_s <- as.data.frame(PDu_v)
simulated_s$PDt_v <- PDt_v
simulated_s$FDu_v <- FDu_v
simulated_s$FDt_v <- FDt_v
simulated_s$score1 <- rep(NA, length(simulated_s$PDu_v))
simulated_s$score2 <- rep(NA, length(simulated_s$PDu_v))
simulated_s$PSA <- rep(NA, length(simulated_s$PDu_v))

# Calculate the scores for each row
for (x in 1:length(simulated_s$PDu_v)) {
  simulated_s$score1[x] <- ((1- (((simulated_s$PDu_v[x]/pd_med)+(simulated_s$FDu_v[x]/fd_med))/2)) + abs((1- (((simulated_s$PDu_v[x]/pd_med)+(simulated_s$FDu_v[x]/fd_med))/2))))/2
  simulated_s$score2[x] <- ((1- (((simulated_s$PDt_v[x]/pd_med)+(simulated_s$FDt_v[x]/fd_med))/2)) + abs((1- (((simulated_s$PDt_v[x]/pd_med)+(simulated_s$FDt_v[x]/fd_med))/2))))/2
  simulated_s$PSA[x] <- simulated_s$score1[x] * simulated_s$score2[x]
}

# Export the scores
write.table(simulated_s, "simulated_PSA_ScaledOnly.txt", sep="\t")

max(simulated_s$score1)
# [1] 0.986129
max(simulated_s$score2)
# [1] 0.972905
max(simulated_s$PSA)
# [1] 0.7323346
# > 
min(simulated_s$score1)
# [1] 0
min(simulated_s$score2)
# [1] 0
min(simulated_s$PSA)
# [1] 0


# Check the distribution of all distances and scores
sim_hist_PDu <- ggplot(simulated_s, aes(x=PDu_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PDu", x="PDu", y = "Count")
sim_hist_FDu <- ggplot(simulated_s, aes(x=FDu_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FDu", x="FDu", y = "Count")
sim_hist_PDt <- ggplot(simulated_s, aes(x=PDt_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PDt", x="PDt", y = "Count")
sim_hist_FDt <- ggplot(simulated_s, aes(x=FDt_v)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="FDt", x="FDt", y = "Count")
sim_hist_s1 <- ggplot(simulated_s, aes(x=score1)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="Score 1", x="Score 1", y = "Count")
sim_hist_s2 <- ggplot(simulated_s, aes(x=score2)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="Score 2", x="Score 2", y = "Count")
sim_hist_PSA <- ggplot(simulated_s, aes(x=PSA)) + 
  geom_histogram(binwidth=0.01)+
  theme_classic()+
  labs(title="PSA", x="PSA", y = "Count")


# check score 1
sim_plot_s1 <- ggplot(simulated_s, aes(x=PDu_v, y=FDu_v, color=cut(score1, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="Score 1")+
  theme(legend.position="bottom")+
  labs(title="Score 1",
       x="PDu", 
       y = "FDu")

# check score 2
sim_plot_s2 <- ggplot(simulated_s, aes(x=PDt_v, y=FDt_v, color=cut(score2, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="Score 2")+
  theme(legend.position="bottom")+
  labs(title="Score 2",
       x="PDt", 
       y = "FDt")

sim_plot_PSA <- ggplot(simulated_s, aes(x=score1, y=score2, color=cut(PSA, c(0,0.2,0.4,0.6,0.8,1), include.lowest=T))) +
  geom_point() +
  coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  scale_color_brewer(type="seq",palette = "YlGn", name="PSA")+
  theme(legend.position="bottom")+
  labs(title="PSA",
       x="Score 1", 
       y = "Score 2")



pdf("PD_FD_Simulated_Scores_distributions_ScaledOnly.pdf", 10, 5) 
grid.arrange(sim_hist_PDu, sim_hist_FDu, sim_hist_PDt, sim_hist_FDt, 
             sim_hist_s1, sim_hist_s2, sim_hist_PSA, nrow = 2)
dev.off()

pdf("PD_FD_Simulated_Scores_plots_ScaledOnly.pdf", 26, 8) 
grid.arrange(sim_plot_s1, sim_plot_s2, sim_plot_PSA, nrow = 1)
dev.off()



#########################







####################################
# Make bar plots for all regions
# Separate uses
# Once with PSA>0.6 (C3), once with PSA>0.8 (C4)
#####################################

# import proportions of species with alternatives according to phylogeny and according to functional traits, and merge them
p <- read.table("Alternatives_Percentages_per_region_cumulative_final_ScaledOnly.txt", sep="\t", header=T)

# make a discrete classification of the number of species

p$NumberUT_cat <- p$NumberUT
for (x in 1:length(p$Region)) {
  if(p$NumberUT[x] > 4) {
    p$NumberUT_cat[x] <- "5+"
  } else {
      p$NumberUT_cat[x] <- "<5"
  }
}

# make a discrete classification of the percentage of alternatives
p$PercentageA_cat <- p$PercentageA
p$PercentageA_cat0 <- p$PercentageA
p$PercentageA_cat25 <- p$PercentageA
p$PercentageA_cat50 <- p$PercentageA
p$PercentageA_cat75 <- p$PercentageA

for (x in 1:length(p$Region)) {
  if (!is.na(p$PercentageA[x])){
    if(p$PercentageA[x] >= 75) {
      p$PercentageA_cat[x] <- "[75-100]"
      p$PercentageA_cat0[x] <- "YES"
      p$PercentageA_cat25[x] <- "YES"
      p$PercentageA_cat50[x] <- "YES"
      p$PercentageA_cat75[x] <- "YES"
    } else {
      if (p$PercentageA[x] >= 50) {
        p$PercentageA_cat[x] <- "[50-75["
        p$PercentageA_cat0[x] <- "YES"
        p$PercentageA_cat25[x] <- "YES"
        p$PercentageA_cat50[x] <- "YES"
        p$PercentageA_cat75[x] <- "NO"
      } else {
        if (p$PercentageA[x] >= 25) {
          p$PercentageA_cat[x] <- "[25-50["
          p$PercentageA_cat0[x] <- "YES"
          p$PercentageA_cat25[x] <- "YES"
          p$PercentageA_cat50[x] <- "NO"
          p$PercentageA_cat75[x] <- "NO"
        } else {
          p$PercentageA_cat[x] <- "[0-25["
          p$PercentageA_cat0[x] <- "YES"
          p$PercentageA_cat25[x] <- "NO"
          p$PercentageA_cat50[x] <- "NO"
          p$PercentageA_cat75[x] <- "NO"
        }
      }
    }
  } else {
    p$PercentageA_cat[x] <- NA
    p$PercentageA_cat0[x] <- NA
    p$PercentageA_cat25[x] <- NA
    p$PercentageA_cat50[x] <- NA
    p$PercentageA_cat75[x] <- NA
  }
}


# Subset p depending on category (focus first on species qualifying for categories 3 and 4, i.e. PSA>0.6)
p60 <- subset(p, p$Category == "C3")

# Subset by scenario (otherwise will be messy to add numbers on plots)
p60_BS <- subset(p60, p60$Scenario == "BS")
p60_WS <- subset(p60, p60$Scenario == "WS")

# Subset by use
p60_BS_HUF <- subset(p60_BS, p60_BS$Use == "HUF")
p60_BS_HUMa <- subset(p60_BS, p60_BS$Use == "HUMa")
p60_BS_HUMe <- subset(p60_BS, p60_BS$Use == "HUMe")
p60_BS_HUC <- subset(p60_BS, p60_BS$Use == "HUC")

p60_WS_HUF <- subset(p60_WS, p60_WS$Use == "HUF")
p60_WS_HUMa <- subset(p60_WS, p60_WS$Use == "HUMa")
p60_WS_HUMe <- subset(p60_WS, p60_WS$Use == "HUMe")
p60_WS_HUC <- subset(p60_WS, p60_WS$Use == "HUC")

# Remove unnecessary columns
p60_BS_HUF_2 <- p60_BS_HUF[,c(8,9)]
p60_BS_HUMa_2 <- p60_BS_HUMa[,c(8,9)]
p60_BS_HUMe_2 <- p60_BS_HUMe[,c(8,9)]
p60_BS_HUC_2 <- p60_BS_HUC[,c(8,9)]

p60_WS_HUF_2 <- p60_WS_HUF[,c(8,9)]
p60_WS_HUMa_2 <- p60_WS_HUMa[,c(8,9)]
p60_WS_HUMe_2 <- p60_WS_HUMe[,c(8,9)]
p60_WS_HUC_2 <- p60_WS_HUC[,c(8,9)]

# rearrange data to be able to plot numbers (make a function first)

reformat <- function(tab) {

  tab_sum <- as.data.frame(table(tab))

  tab_sum_s <- arrange(tab_sum, NumberUT_cat, desc(PercentageA_cat))

  tab_sum_s_prop <- ddply(tab_sum_s, "NumberUT_cat", transform, prop=proportions(Freq))

  tab_sum_s_prop$perc <- round(tab_sum_s_prop$prop*100)

  tab_sum_s_cumper <- ddply(tab_sum_s_prop, "NumberUT_cat", transform, laby_pos=cumsum(perc))

  tab_sum_s_cumper$laby_lab <- paste(tab_sum_s_cumper$perc, "%", " (", tab_sum_s_cumper$Freq , ")", sep="")

return(tab_sum_s_cumper)

}

p60_BS_HUF_3 <- reformat(p60_BS_HUF_2)
p60_BS_HUMa_3 <- reformat(p60_BS_HUMa_2)
p60_BS_HUMe_3 <- reformat(p60_BS_HUMe_2)
p60_BS_HUC_3 <- reformat(p60_BS_HUC_2)

p60_WS_HUF_3 <- reformat(p60_WS_HUF_2)
p60_WS_HUMa_3 <- reformat(p60_WS_HUMa_2)
p60_WS_HUMe_3 <- reformat(p60_WS_HUMe_2)
p60_WS_HUC_3 <- reformat(p60_WS_HUC_2)


# Combine BS and WS
p60_BS_HUF_3$Scenario <- rep("Best", length(p60_BS_HUF_3$NumberUT_cat))
p60_BS_HUMa_3$Scenario <- rep("Best", length(p60_BS_HUMa_3$NumberUT_cat))
p60_BS_HUMe_3$Scenario <- rep("Best", length(p60_BS_HUMe_3$NumberUT_cat))
p60_BS_HUC_3$Scenario <- rep("Best", length(p60_BS_HUC_3$NumberUT_cat))

p60_WS_HUF_3$Scenario <- rep("Worst", length(p60_WS_HUF_3$NumberUT_cat))
p60_WS_HUMa_3$Scenario <- rep("Worst", length(p60_WS_HUMa_3$NumberUT_cat))
p60_WS_HUMe_3$Scenario <- rep("Worst", length(p60_WS_HUMe_3$NumberUT_cat))
p60_WS_HUC_3$Scenario <- rep("Worst", length(p60_WS_HUC_3$NumberUT_cat))

p60_BSWS_HUF <- rbind(p60_BS_HUF_3, p60_WS_HUF_3)
p60_BSWS_HUMa <- rbind(p60_BS_HUMa_3, p60_WS_HUMa_3)
p60_BSWS_HUMe <- rbind(p60_BS_HUMe_3, p60_WS_HUMe_3)
p60_BSWS_HUC <- rbind(p60_BS_HUC_3, p60_WS_HUC_3)

# plot
# CAREFUL colors order needs customisation because changes if all classes are not in the first bar
# ALSO CAREFUL because even when color is customised, order of some numbers is wrong: needs manual correction in correl

p_HUF <- ggplot(data=p60_BSWS_HUF, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#7B2382FF", "#410F75FF", "#0C0927FF"), name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Food")

p_HUMa <- ggplot(data=p60_BSWS_HUMa, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#410F75FF", "#0C0927FF", "#7B2382FF"), name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Utensils, tools and construction")

p_HUMe <- ggplot(data=p60_BSWS_HUMe, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#410F75FF", "#0C0927FF", "#7B2382FF") , name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Medicine")

p_HUC <- ggplot(data=p60_BSWS_HUC, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#410F75FF", "#0C0927FF", "#7B2382FF"), name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Cultural uses")


# Export plots
pdf("Fig3A_Percentage_regions_PSA60.pdf", 3, 4) 
grid.arrange(p_HUF, p_HUMa, p_HUMe, p_HUC, nrow = 2)
dev.off()



# same for PSA>0.8

# Subset p depending on category
p60 <- subset(p, p$Category == "C4") # do not bother changing df name so stay p60

# Subset by scenario (otherwise will be messy to add numbers on plots)
p60_BS <- subset(p60, p60$Scenario == "BS")
p60_WS <- subset(p60, p60$Scenario == "WS")

# Subset by use
p60_BS_HUF <- subset(p60_BS, p60_BS$Use == "HUF")
p60_BS_HUMa <- subset(p60_BS, p60_BS$Use == "HUMa")
p60_BS_HUMe <- subset(p60_BS, p60_BS$Use == "HUMe")
p60_BS_HUC <- subset(p60_BS, p60_BS$Use == "HUC")

p60_WS_HUF <- subset(p60_WS, p60_WS$Use == "HUF")
p60_WS_HUMa <- subset(p60_WS, p60_WS$Use == "HUMa")
p60_WS_HUMe <- subset(p60_WS, p60_WS$Use == "HUMe")
p60_WS_HUC <- subset(p60_WS, p60_WS$Use == "HUC")

# Remove unnecessary columns
p60_BS_HUF_2 <- p60_BS_HUF[,c(8,9)]
p60_BS_HUMa_2 <- p60_BS_HUMa[,c(8,9)]
p60_BS_HUMe_2 <- p60_BS_HUMe[,c(8,9)]
p60_BS_HUC_2 <- p60_BS_HUC[,c(8,9)]

p60_WS_HUF_2 <- p60_WS_HUF[,c(8,9)]
p60_WS_HUMa_2 <- p60_WS_HUMa[,c(8,9)]
p60_WS_HUMe_2 <- p60_WS_HUMe[,c(8,9)]
p60_WS_HUC_2 <- p60_WS_HUC[,c(8,9)]

# rearrange data to be able to plot numbers (make a function first)

p60_BS_HUF_3 <- reformat(p60_BS_HUF_2)
p60_BS_HUMa_3 <- reformat(p60_BS_HUMa_2)
p60_BS_HUMe_3 <- reformat(p60_BS_HUMe_2)
p60_BS_HUC_3 <- reformat(p60_BS_HUC_2)

p60_WS_HUF_3 <- reformat(p60_WS_HUF_2)
p60_WS_HUMa_3 <- reformat(p60_WS_HUMa_2)
p60_WS_HUMe_3 <- reformat(p60_WS_HUMe_2)
p60_WS_HUC_3 <- reformat(p60_WS_HUC_2)


# Combine BS and WS
p60_BS_HUF_3$Scenario <- rep("Best", length(p60_BS_HUF_3$NumberUT_cat))
p60_BS_HUMa_3$Scenario <- rep("Best", length(p60_BS_HUMa_3$NumberUT_cat))
p60_BS_HUMe_3$Scenario <- rep("Best", length(p60_BS_HUMe_3$NumberUT_cat))
p60_BS_HUC_3$Scenario <- rep("Best", length(p60_BS_HUC_3$NumberUT_cat))

p60_WS_HUF_3$Scenario <- rep("Worst", length(p60_WS_HUF_3$NumberUT_cat))
p60_WS_HUMa_3$Scenario <- rep("Worst", length(p60_WS_HUMa_3$NumberUT_cat))
p60_WS_HUMe_3$Scenario <- rep("Worst", length(p60_WS_HUMe_3$NumberUT_cat))
p60_WS_HUC_3$Scenario <- rep("Worst", length(p60_WS_HUC_3$NumberUT_cat))

p60_BSWS_HUF <- rbind(p60_BS_HUF_3, p60_WS_HUF_3)
p60_BSWS_HUMa <- rbind(p60_BS_HUMa_3, p60_WS_HUMa_3)
p60_BSWS_HUMe <- rbind(p60_BS_HUMe_3, p60_WS_HUMe_3)
p60_BSWS_HUC <- rbind(p60_BS_HUC_3, p60_WS_HUC_3)

# plot
# CAREFUL colors order needs customisation because changes if all classes are not in the first bar
# ALSO CAREFUL because even when color is customised, order of numbers is wrong: needs manual correction in correl

p_HUF <- ggplot(data=p60_BSWS_HUF, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#7B2382FF", "#410F75FF", "#0C0927FF"), name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Food")

p_HUMa <- ggplot(data=p60_BSWS_HUMa, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#7B2382FF", "#410F75FF", "#0C0927FF"), name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Utensils, tools and construction")

p_HUMe <- ggplot(data=p60_BSWS_HUMe, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#7B2382FF", "#410F75FF", "#0C0927FF") , name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Medicine")

p_HUC <- ggplot(data=p60_BSWS_HUC, aes(x=Scenario, y=perc, fill=PercentageA_cat)) +
  geom_bar(stat="identity") +
  facet_grid(~NumberUT_cat)+
  geom_text(aes(y =laby_pos, label=laby_lab), vjust=1.5, color="grey", size=1)+
  scale_fill_manual(values=c("#D3436EFF", "#7B2382FF", "#410F75FF", "#0C0927FF"), name = "Percentage of utilized threatened species \nwith at least one alternative")+
  #scale_y_continuous(breaks=seq(0,50,10), limits=c(0,50))+
  theme_classic() +
  theme(legend.position="bottom", axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.x=element_text(angle = 0, hjust = 0.5, vjust = 0.25, size=5), axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=5), legend.text=element_text(size=5), legend.title=element_text(size=5),
        legend.key.size = unit(0.5,"line")) +
  labs(x = "", y = "Percentage of regions", title="Cultural uses")





# Export plots
pdf("Fig3A_Percentage_regions_PSA80.pdf", 3, 4) 
grid.arrange(p_HUF, p_HUMa, p_HUMe, p_HUC, nrow = 2)
dev.off()
#########################################






#############################
# Make a summary plot, faceted by species richness
# Build subtables separately for ease and then bring them together
###################################

# import proportions of species with alternatives according to phylogeny and according to functional traits, and merge them
p <- read.table("Alternatives_Percentages_per_region_cumulative_final_ScaledOnly.txt", sep="\t", header=T)

# make a discrete classification of the number of species

p$NumberUT_cat <- p$NumberUT
for (x in 1:length(p$Region)) {
  if(p$NumberUT[x] > 4) {
    p$NumberUT_cat[x] <- "5+"
  } else {
    p$NumberUT_cat[x] <- "<5"
  }
}

# make a discrete classification of the percentage of alternatives
p$PercentageA_cat <- p$PercentageA

for (x in 1:length(p$Region)) {
  if (!is.na(p$PercentageA[x])){
    if(p$PercentageA[x] <50) {
      p$PercentageA_cat[x] <- "[0-50["
    } else {
      p$PercentageA_cat[x] <- "[50-100]"
    } 
  } else {
    p$PercentageA_cat[x] <- NA
  }
}


# Subset p depending on category (focus first on species qualifying for categories 3 and 4, i.e. PSA>0.6)
p60 <- subset(p, p$Category == "C3")

# Subset by scenario
p60_BS <- subset(p60, p60$Scenario == "BS")
p60_WS <- subset(p60, p60$Scenario == "WS")

# Subset by use
p60_BS_HUF <- subset(p60_BS, p60_BS$Use == "HUF")
p60_BS_HUMa <- subset(p60_BS, p60_BS$Use == "HUMa")
p60_BS_HUMe <- subset(p60_BS, p60_BS$Use == "HUMe")
p60_BS_HUC <- subset(p60_BS, p60_BS$Use == "HUC")

p60_WS_HUF <- subset(p60_WS, p60_WS$Use == "HUF")
p60_WS_HUMa <- subset(p60_WS, p60_WS$Use == "HUMa")
p60_WS_HUMe <- subset(p60_WS, p60_WS$Use == "HUMe")
p60_WS_HUC <- subset(p60_WS, p60_WS$Use == "HUC")

# Remove unnecessary columns
p60_BS_HUF_2 <- p60_BS_HUF[,c(8,9)]
p60_BS_HUMa_2 <- p60_BS_HUMa[,c(8,9)]
p60_BS_HUMe_2 <- p60_BS_HUMe[,c(8,9)]
p60_BS_HUC_2 <- p60_BS_HUC[,c(8,9)]

p60_WS_HUF_2 <- p60_WS_HUF[,c(8,9)]
p60_WS_HUMa_2 <- p60_WS_HUMa[,c(8,9)]
p60_WS_HUMe_2 <- p60_WS_HUMe[,c(8,9)]
p60_WS_HUC_2 <- p60_WS_HUC[,c(8,9)]

# rearrange data to be able to plot numbers (make a function first)

reformat <- function(tab) {
  
  tab_sum <- as.data.frame(table(tab))
  
  tab_sum_s <- arrange(tab_sum, NumberUT_cat, desc(PercentageA_cat))
  
  tab_sum_s_prop <- ddply(tab_sum_s, "NumberUT_cat", transform, prop=proportions(Freq))
  
  tab_sum_s_prop$perc <- round(tab_sum_s_prop$prop*100)
  
  tab_sum_s_cumper <- ddply(tab_sum_s_prop, "NumberUT_cat", transform, laby_pos=cumsum(perc))
  
  tab_sum_s_cumper$laby_lab <- paste(tab_sum_s_cumper$perc, "%", " (", tab_sum_s_cumper$Freq , ")", sep="")
  
  return(tab_sum_s_cumper)
  
}

p60_BS_HUF_3 <- reformat(p60_BS_HUF_2)
p60_BS_HUMa_3 <- reformat(p60_BS_HUMa_2)
p60_BS_HUMe_3 <- reformat(p60_BS_HUMe_2)
p60_BS_HUC_3 <- reformat(p60_BS_HUC_2)

p60_WS_HUF_3 <- reformat(p60_WS_HUF_2)
p60_WS_HUMa_3 <- reformat(p60_WS_HUMa_2)
p60_WS_HUMe_3 <- reformat(p60_WS_HUMe_2)
p60_WS_HUC_3 <- reformat(p60_WS_HUC_2)


# Combine BS and WS
p60_BS_HUF_3$Scenario <- rep("Best", length(p60_BS_HUF_3$NumberUT_cat))
p60_BS_HUMa_3$Scenario <- rep("Best", length(p60_BS_HUMa_3$NumberUT_cat))
p60_BS_HUMe_3$Scenario <- rep("Best", length(p60_BS_HUMe_3$NumberUT_cat))
p60_BS_HUC_3$Scenario <- rep("Best", length(p60_BS_HUC_3$NumberUT_cat))

p60_WS_HUF_3$Scenario <- rep("Worst", length(p60_WS_HUF_3$NumberUT_cat))
p60_WS_HUMa_3$Scenario <- rep("Worst", length(p60_WS_HUMa_3$NumberUT_cat))
p60_WS_HUMe_3$Scenario <- rep("Worst", length(p60_WS_HUMe_3$NumberUT_cat))
p60_WS_HUC_3$Scenario <- rep("Worst", length(p60_WS_HUC_3$NumberUT_cat))

p60_BSWS_HUF <- rbind(p60_BS_HUF_3, p60_WS_HUF_3)
p60_BSWS_HUMa <- rbind(p60_BS_HUMa_3, p60_WS_HUMa_3)
p60_BSWS_HUMe <- rbind(p60_BS_HUMe_3, p60_WS_HUMe_3)
p60_BSWS_HUC <- rbind(p60_BS_HUC_3, p60_WS_HUC_3)

# Combine Uses
p60_BSWS_HUF$Use <- rep("Food", length(p60_BSWS_HUF$NumberUT_cat))
p60_BSWS_HUMa$Use <- rep("Utensils", length(p60_BSWS_HUMa$NumberUT_cat))
p60_BSWS_HUMe$Use <- rep("Medicine", length(p60_BSWS_HUMe$NumberUT_cat))
p60_BSWS_HUC$Use <- rep("Culture", length(p60_BSWS_HUC$NumberUT_cat))

p60_BSWS_AllU <- rbind(p60_BSWS_HUF, p60_BSWS_HUMa, p60_BSWS_HUMe, p60_BSWS_HUC)

p60_BSWS_AllU$Score <- rep("06", length(p60_BSWS_AllU$NumberUT_cat))

# Do the same for PSA>0.8 (C4)

# Subset p depending on category 
p80 <- subset(p, p$Category == "C4")

# Subset by scenario
p80_BS <- subset(p80, p80$Scenario == "BS")
p80_WS <- subset(p80, p80$Scenario == "WS")

# Subset by use
p80_BS_HUF <- subset(p80_BS, p80_BS$Use == "HUF")
p80_BS_HUMa <- subset(p80_BS, p80_BS$Use == "HUMa")
p80_BS_HUMe <- subset(p80_BS, p80_BS$Use == "HUMe")
p80_BS_HUC <- subset(p80_BS, p80_BS$Use == "HUC")

p80_WS_HUF <- subset(p80_WS, p80_WS$Use == "HUF")
p80_WS_HUMa <- subset(p80_WS, p80_WS$Use == "HUMa")
p80_WS_HUMe <- subset(p80_WS, p80_WS$Use == "HUMe")
p80_WS_HUC <- subset(p80_WS, p80_WS$Use == "HUC")

# Remove unnecessary columns
p80_BS_HUF_2 <- p80_BS_HUF[,c(8,9)]
p80_BS_HUMa_2 <- p80_BS_HUMa[,c(8,9)]
p80_BS_HUMe_2 <- p80_BS_HUMe[,c(8,9)]
p80_BS_HUC_2 <- p80_BS_HUC[,c(8,9)]

p80_WS_HUF_2 <- p80_WS_HUF[,c(8,9)]
p80_WS_HUMa_2 <- p80_WS_HUMa[,c(8,9)]
p80_WS_HUMe_2 <- p80_WS_HUMe[,c(8,9)]
p80_WS_HUC_2 <- p80_WS_HUC[,c(8,9)]

# rearrange data to be able to plot numbers (make a function first)

reformat <- function(tab) {
  
  tab_sum <- as.data.frame(table(tab))
  
  tab_sum_s <- arrange(tab_sum, NumberUT_cat, desc(PercentageA_cat))
  
  tab_sum_s_prop <- ddply(tab_sum_s, "NumberUT_cat", transform, prop=proportions(Freq))
  
  tab_sum_s_prop$perc <- round(tab_sum_s_prop$prop*100)
  
  tab_sum_s_cumper <- ddply(tab_sum_s_prop, "NumberUT_cat", transform, laby_pos=cumsum(perc))
  
  tab_sum_s_cumper$laby_lab <- paste(tab_sum_s_cumper$perc, "%", " (", tab_sum_s_cumper$Freq , ")", sep="")
  
  return(tab_sum_s_cumper)
  
}

p80_BS_HUF_3 <- reformat(p80_BS_HUF_2)
p80_BS_HUMa_3 <- reformat(p80_BS_HUMa_2)
p80_BS_HUMe_3 <- reformat(p80_BS_HUMe_2)
p80_BS_HUC_3 <- reformat(p80_BS_HUC_2)

p80_WS_HUF_3 <- reformat(p80_WS_HUF_2)
p80_WS_HUMa_3 <- reformat(p80_WS_HUMa_2)
p80_WS_HUMe_3 <- reformat(p80_WS_HUMe_2)
p80_WS_HUC_3 <- reformat(p80_WS_HUC_2)


# Combine BS and WS
p80_BS_HUF_3$Scenario <- rep("Best", length(p80_BS_HUF_3$NumberUT_cat))
p80_BS_HUMa_3$Scenario <- rep("Best", length(p80_BS_HUMa_3$NumberUT_cat))
p80_BS_HUMe_3$Scenario <- rep("Best", length(p80_BS_HUMe_3$NumberUT_cat))
p80_BS_HUC_3$Scenario <- rep("Best", length(p80_BS_HUC_3$NumberUT_cat))

p80_WS_HUF_3$Scenario <- rep("Worst", length(p80_WS_HUF_3$NumberUT_cat))
p80_WS_HUMa_3$Scenario <- rep("Worst", length(p80_WS_HUMa_3$NumberUT_cat))
p80_WS_HUMe_3$Scenario <- rep("Worst", length(p80_WS_HUMe_3$NumberUT_cat))
p80_WS_HUC_3$Scenario <- rep("Worst", length(p80_WS_HUC_3$NumberUT_cat))

p80_BSWS_HUF <- rbind(p80_BS_HUF_3, p80_WS_HUF_3)
p80_BSWS_HUMa <- rbind(p80_BS_HUMa_3, p80_WS_HUMa_3)
p80_BSWS_HUMe <- rbind(p80_BS_HUMe_3, p80_WS_HUMe_3)
p80_BSWS_HUC <- rbind(p80_BS_HUC_3, p80_WS_HUC_3)

# Combine Uses
p80_BSWS_HUF$Use <- rep("Food", length(p80_BSWS_HUF$NumberUT_cat))
p80_BSWS_HUMa$Use <- rep("Utensils", length(p80_BSWS_HUMa$NumberUT_cat))
p80_BSWS_HUMe$Use <- rep("Medicine", length(p80_BSWS_HUMe$NumberUT_cat))
p80_BSWS_HUC$Use <- rep("Culture", length(p80_BSWS_HUC$NumberUT_cat))

p80_BSWS_AllU <- rbind(p80_BSWS_HUF, p80_BSWS_HUMa, p80_BSWS_HUMe, p80_BSWS_HUC)

p80_BSWS_AllU$Score <- rep("08", length(p80_BSWS_AllU$NumberUT_cat))


# combine both p60 and p80
p6080_BSWS_AllU <- rbind(p60_BSWS_AllU, p80_BSWS_AllU)

# keep only [0-50] and necessary columns

p6080_BSWS_AllU_0 <- subset(p6080_BSWS_AllU, p6080_BSWS_AllU$PercentageA_cat %in% "[0-50[")
p6080_BSWS_AllU_0$ScenarioScore <- paste(p6080_BSWS_AllU_0$Scenario, p6080_BSWS_AllU_0$Score, sep="")

p6080_BSWS_AllU_0_2 <- p6080_BSWS_AllU_0[,c(1,3,5,9,10,11)]

# create columns to draw segments between some BS and WS of a same score, category and use
p6080_BSWS_AllU_0_2$Startx <- rep(NA, length(p6080_BSWS_AllU_0_2$NumberUT_cat))
p6080_BSWS_AllU_0_2$Starty <- rep(NA, length(p6080_BSWS_AllU_0_2$NumberUT_cat))
p6080_BSWS_AllU_0_2$Endx <- rep(NA, length(p6080_BSWS_AllU_0_2$NumberUT_cat))
p6080_BSWS_AllU_0_2$Endy <- rep(NA, length(p6080_BSWS_AllU_0_2$NumberUT_cat))

# Fill column

for (x in 1:length(p6080_BSWS_AllU_0_2$NumberUT_cat)) {
  if ( length(grep("Best", p6080_BSWS_AllU_0_2$ScenarioScore[x]))>0) {
    if ( length(grep("Medicine", p6080_BSWS_AllU_0_2$Use[x]))>0) {
  p6080_BSWS_AllU_0_2$Startx[x] <- p6080_BSWS_AllU_0_2$Freq[x]
  p6080_BSWS_AllU_0_2$Starty[x] <- p6080_BSWS_AllU_0_2$perc[x]
  p6080_BSWS_AllU_0_2$Endx[x] <- p6080_BSWS_AllU_0_2$Freq[x+1]
  p6080_BSWS_AllU_0_2$Endy[x] <- p6080_BSWS_AllU_0_2$perc[x+1]
    } else {
      p6080_BSWS_AllU_0_2$Startx[x] <- p6080_BSWS_AllU_0_2$Freq[x]
      p6080_BSWS_AllU_0_2$Starty[x] <- p6080_BSWS_AllU_0_2$perc[x]
      p6080_BSWS_AllU_0_2$Endx[x] <- p6080_BSWS_AllU_0_2$Freq[x+2]
      p6080_BSWS_AllU_0_2$Endy[x] <- p6080_BSWS_AllU_0_2$perc[x+2]      
    }
  }
}


# plot

  

write.table(p6080_BSWS_AllU_0_2, "Regions_digest_for_plot.txt", sep="\t", row.names = F)

p6080_BSWS_AllU_0_3 <- read.table("Regions_digest_for_plot.txt", header = T)

plot <- ggplot(p6080_BSWS_AllU_0_3, aes(x=Freq, y=perc, color=ScenarioScore)) +
  geom_point( size=8, shape = 19 )+
  scale_color_manual(values=c("#0072B2", "#56B4E9", "#D55E00", "#F0E442" ))+
  #scale_shape_manual(values=c("\U0001F9E0", "\U0001f958", "\U0001f48A", "\U0001F6E0" ))+   # do not differentiate uses - see below
  geom_segment( aes(x=Startx, xend=Endx, y=Starty, yend=Endy), color="grey", size=1) +
  facet_grid(~ NumberUT_cat) +
  theme_classic() +
  theme(legend.position="bottom")+
  xlab("Number of regions with lower resilience") +
  ylab("Percentage of regions with lower resilience")


pdf("Figure3B_Summary.pdf", 7, 5) 
plot
dev.off()

# Try to use unicodes but not working well and not sure about rights of use so just use it to check where each use falls on the plot and will replace by nice images in correl
# Export that one in svg from the plot window
plot2 <- ggplot(p6080_BSWS_AllU_0_3, aes(x=Freq, y=perc, color=ScenarioScore, shape=Use)) +
  geom_point( size=8 )+
  scale_color_manual(values=c("#0072B2", "#56B4E9", "#D55E00", "#F0E442" ))+
  scale_shape_manual(values=c("\U0001F9E0", "\U0001f958", "\U0001f48A", "\U0001F6E0" ))+
  geom_segment( aes(x=Startx, xend=Endx, y=Starty, yend=Endy), color="grey", size=1) +
  facet_grid(~ NumberUT_cat) +
  theme_classic() +
  xlab("Number of regions with lower resilience") +
  ylab("Percentage of regions with lower resilience")






















##################################
# Make maps of percentage of UT species with at lest 1 alternative
# Only BS (can see on plot above that WS is similar)
# Do PSA 0.6 and PSA 0.8 (the latter for OSM)
# Separate for each use
# Diverging colors according to species richness category
##################################

# Get region table from above 

p60 <- subset(p, p$Category == "C3")

# Subset by scenario
p60_BS <- subset(p60, p60$Scenario == "BS")

# Subset by use
p60_BS_HUF <- subset(p60_BS, p60_BS$Use == "HUF")
p60_BS_HUMa <- subset(p60_BS, p60_BS$Use == "HUMa")
p60_BS_HUMe <- subset(p60_BS, p60_BS$Use == "HUMe")
p60_BS_HUC <- subset(p60_BS, p60_BS$Use == "HUC")

# Code a new variable combining UT species richness and the percentage of alternatives
p60_BS_HUF$NumUT_PerA <- paste(p60_BS_HUF$NumberUT_cat, p60_BS_HUF$PercentageA_cat, sep = "_")
p60_BS_HUMa$NumUT_PerA <- paste(p60_BS_HUMa$NumberUT_cat, p60_BS_HUMa$PercentageA_cat, sep = "_")
p60_BS_HUMe$NumUT_PerA <- paste(p60_BS_HUMe$NumberUT_cat, p60_BS_HUMe$PercentageA_cat, sep = "_")
p60_BS_HUC$NumUT_PerA <- paste(p60_BS_HUC$NumberUT_cat, p60_BS_HUC$PercentageA_cat, sep = "_")

# Recode the ones with 0 UT species to be simply NA
for (x in 1:length(p60_BS_HUF$Region)) {
  if (p60_BS_HUF$NumberUT[x] == 0) {
    p60_BS_HUF$NumUT_PerA[x] <- NA
  }
}
for (x in 1:length(p60_BS_HUMa$Region)) {
  if (p60_BS_HUMa$NumberUT[x] == 0) {
    p60_BS_HUMa$NumUT_PerA[x] <- NA
  }
}
for (x in 1:length(p60_BS_HUMe$Region)) {
  if (p60_BS_HUMe$NumberUT[x] == 0) {
    p60_BS_HUMe$NumUT_PerA[x] <- NA
  }
}
for (x in 1:length(p60_BS_HUC$Region)) {
  if (p60_BS_HUC$NumberUT[x] == 0) {
    p60_BS_HUC$NumUT_PerA[x] <- NA
  }
}

# make dictionaries

p60_BS_HUF_h <- hash(keys=p60_BS_HUF$Region, values=p60_BS_HUF$NumUT_PerA)
p60_BS_HUMa_h <- hash(keys=p60_BS_HUMa$Region, values=p60_BS_HUMa$NumUT_PerA)
p60_BS_HUMe_h <- hash(keys=p60_BS_HUMe$Region, values=p60_BS_HUMe$NumUT_PerA)
p60_BS_HUC_h <- hash(keys=p60_BS_HUC$Region, values=p60_BS_HUC$NumUT_PerA)

# Map
# create column in the raster
tdwg3$LEVEL3_SPE_HUF_ALT60BS <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_ALT60BS <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_ALT60BS <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_ALT60BS <- rep(NA, length(tdwg3$LEVEL3_NAM))

# populate columns 
for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUF_h)) {
    tdwg3$LEVEL3_SPE_HUF_ALT60BS[x] <- p60_BS_HUF_h[[tdwg3$LEVEL3_COD[x]]]
  }
}

for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUMa_h)) {
    tdwg3$LEVEL3_SPE_HUMa_ALT60BS[x] <- p60_BS_HUMa_h[[tdwg3$LEVEL3_COD[x]]]
  }
}

for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUMe_h)) {
    tdwg3$LEVEL3_SPE_HUMe_ALT60BS[x] <- p60_BS_HUMe_h[[tdwg3$LEVEL3_COD[x]]]
  }
}

for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUC_h)) {
    tdwg3$LEVEL3_SPE_HUC_ALT60BS[x] <- p60_BS_HUC_h[[tdwg3$LEVEL3_COD[x]]]
  }
}




# plot

# viridis colors for 5 categories: "#0D0887FF" "#7E03A8FF" "#CC4678FF" "#F89441FF" "#F0F921FF"
# need to change them manually below depending on what categories are in the data
# check the categories:

summary(as.factor(tdwg3$LEVEL3_SPE_HUF_ALT60BS))
summary(as.factor(tdwg3$LEVEL3_SPE_HUMa_ALT60BS))
summary(as.factor(tdwg3$LEVEL3_SPE_HUMe_ALT60BS))
summary(as.factor(tdwg3$LEVEL3_SPE_HUC_ALT60BS))

# c("#b35806", "#e08214", "#fdb863", "#fee0b6", 
# "#d8daeb", "#b2abd2", "#8073ac", "#542788")



# plot the map for best case scenario (legends are plotted in separate plot below but first try with legend to 
# make sure the colors are right!!!!!)

map_HUF_ALT60BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUF_ALT60BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#e08214", "#fdb863", "#fee0b6", "#b2abd2", "#d8daeb"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of \"interesting\" species known or predicted to be threatened")+ 
  ggtitle("Food")+
  theme(legend.position="none")

map_HUMa_ALT60BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUMa_ALT60BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#fdb863", "#fee0b6", "#b2abd2", "#d8daeb"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of \"interesting\" species known or predicted to be threatened")+ 
  ggtitle("Utensils, tools and construction")+
  theme(legend.position="none")

map_HUMe_ALT60BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUMe_ALT60BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#fdb863", "#fee0b6"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of \"interesting\" species known or predicted to be threatened")+ 
  ggtitle("Medicine")+
  theme(legend.position="none")

map_HUC_ALT60BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUC_ALT60BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#fdb863", "#fee0b6", "#d8daeb"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of \"interesting\" species known or predicted to be threatened")+ 
  ggtitle("Cultural uses")+
  theme(legend.position="none")



pdf("Fig3_Maps_p60.pdf", 16, 8) 
grid.arrange(map_HUF_ALT60BS, map_HUMa_ALT60BS, map_HUMe_ALT60BS, map_HUC_ALT60BS, nrow=2)
dev.off()


### Same for PSA>0.8 (C4; keep it called p60 below for simplicity)


# Get region table from above 

p60 <- subset(p, p$Category == "C4")

# Subset by scenario
p60_BS <- subset(p60, p60$Scenario == "BS")

# Subset by use
p60_BS_HUF <- subset(p60_BS, p60_BS$Use == "HUF")
p60_BS_HUMa <- subset(p60_BS, p60_BS$Use == "HUMa")
p60_BS_HUMe <- subset(p60_BS, p60_BS$Use == "HUMe")
p60_BS_HUC <- subset(p60_BS, p60_BS$Use == "HUC")

# Code a new variable combining UT species richness and the percentage of alternatives
p60_BS_HUF$NumUT_PerA <- paste(p60_BS_HUF$NumberUT_cat, p60_BS_HUF$PercentageA_cat, sep = "_")
p60_BS_HUMa$NumUT_PerA <- paste(p60_BS_HUMa$NumberUT_cat, p60_BS_HUMa$PercentageA_cat, sep = "_")
p60_BS_HUMe$NumUT_PerA <- paste(p60_BS_HUMe$NumberUT_cat, p60_BS_HUMe$PercentageA_cat, sep = "_")
p60_BS_HUC$NumUT_PerA <- paste(p60_BS_HUC$NumberUT_cat, p60_BS_HUC$PercentageA_cat, sep = "_")

# Recode the ones with 0 UT species to be simply NA
for (x in 1:length(p60_BS_HUF$Region)) {
  if (p60_BS_HUF$NumberUT[x] == 0) {
    p60_BS_HUF$NumUT_PerA[x] <- NA
  }
}
for (x in 1:length(p60_BS_HUMa$Region)) {
  if (p60_BS_HUMa$NumberUT[x] == 0) {
    p60_BS_HUMa$NumUT_PerA[x] <- NA
  }
}
for (x in 1:length(p60_BS_HUMe$Region)) {
  if (p60_BS_HUMe$NumberUT[x] == 0) {
    p60_BS_HUMe$NumUT_PerA[x] <- NA
  }
}
for (x in 1:length(p60_BS_HUC$Region)) {
  if (p60_BS_HUC$NumberUT[x] == 0) {
    p60_BS_HUC$NumUT_PerA[x] <- NA
  }
}

# make dictionaries

p60_BS_HUF_h <- hash(keys=p60_BS_HUF$Region, values=p60_BS_HUF$NumUT_PerA)
p60_BS_HUMa_h <- hash(keys=p60_BS_HUMa$Region, values=p60_BS_HUMa$NumUT_PerA)
p60_BS_HUMe_h <- hash(keys=p60_BS_HUMe$Region, values=p60_BS_HUMe$NumUT_PerA)
p60_BS_HUC_h <- hash(keys=p60_BS_HUC$Region, values=p60_BS_HUC$NumUT_PerA)

# Map

# create column in the raster (CHANGED NAME FOR A DIFFERENT PSA, to not overwrite PSA>0.6)
tdwg3$LEVEL3_SPE_HUF_ALT80BS <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMa_ALT80BS <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUMe_ALT80BS <- rep(NA, length(tdwg3$LEVEL3_NAM))
tdwg3$LEVEL3_SPE_HUC_ALT80BS <- rep(NA, length(tdwg3$LEVEL3_NAM))

# populate columns (changed here too depending on PSA)
for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUF_h)) {
    tdwg3$LEVEL3_SPE_HUF_ALT80BS[x] <- p60_BS_HUF_h[[tdwg3$LEVEL3_COD[x]]]
  }
}

for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUMa_h)) {
    tdwg3$LEVEL3_SPE_HUMa_ALT80BS[x] <- p60_BS_HUMa_h[[tdwg3$LEVEL3_COD[x]]]
  }
}

for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUMe_h)) {
    tdwg3$LEVEL3_SPE_HUMe_ALT80BS[x] <- p60_BS_HUMe_h[[tdwg3$LEVEL3_COD[x]]]
  }
}

for (x in 1:length(tdwg3$LEVEL3_NAM)){
  if (tdwg3$LEVEL3_COD[x] %in% keys(p60_BS_HUC_h)) {
    tdwg3$LEVEL3_SPE_HUC_ALT80BS[x] <- p60_BS_HUC_h[[tdwg3$LEVEL3_COD[x]]]
  }
}




# plot

# viridis colors for 5 categories: "#0D0887FF" "#7E03A8FF" "#CC4678FF" "#F89441FF" "#F0F921FF"
# need to change them manually below depending on what categories are in the data
# check the categories:

# select right column here and below depending on PSA
summary(as.factor(tdwg3$LEVEL3_SPE_HUF_ALT80BS))
summary(as.factor(tdwg3$LEVEL3_SPE_HUMa_ALT80BS))
summary(as.factor(tdwg3$LEVEL3_SPE_HUMe_ALT80BS))
summary(as.factor(tdwg3$LEVEL3_SPE_HUC_ALT80BS))

# c("#b35806", "#e08214", "#fdb863", "#fee0b6", 
# "#d8daeb", "#b2abd2", "#8073ac", "#542788")


# plot the map for best case scenario (legends are plotted in separate plot below but first try with legend to 
# make sure the colors are right!!!!!)
map_HUF_ALT80BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUF_ALT80BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#e08214", "#fdb863", "#fee0b6", "#542788", "#b2abd2", "#d8daeb"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of threatened utilized species\nwith alternatives")+ 
  ggtitle("Food")+
  theme(legend.position="none")

map_HUMa_ALT80BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUMa_ALT80BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#e08214", "#fdb863", "#fee0b6", "#542788", "#8073ac", "#b2abd2", "#d8daeb"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of threatened utilized species\nwith alternatives")+ 
  ggtitle("Utensils, tools and construction")+
  theme(legend.position="none")

map_HUMe_ALT80BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUMe_ALT80BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#e08214", "#fdb863", "#fee0b6"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of threatened utilized species\nwith alternatives")+ 
  ggtitle("Medicine")+
  theme(legend.position="none")

map_HUC_ALT80BS <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUC_ALT80BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#e08214", "#fdb863", "#fee0b6", "#b2abd2"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of threatened utilized species\nwith alternatives")+ 
  ggtitle("Cultural uses")+
  theme(legend.position="none")


pdf("Fig3_Maps_p80.pdf", 16, 8) 
grid.arrange(map_HUF_ALT80BS, map_HUMa_ALT80BS, map_HUMe_ALT80BS, map_HUC_ALT80BS, nrow=2)
dev.off()


# make a legend plot
map_HUMa_ALT80BS_legend <- ggplot(tdwg3) +
  geom_sf(aes(fill = LEVEL3_SPE_HUMa_ALT80BS), color=NA) + # remove color=NA if want country borders
  scale_fill_manual(values = c("#b35806", "#e08214", "#fdb863", "#fee0b6", "#542788", "#8073ac", "#b2abd2", "#d8daeb"), na.value="grey") +
  theme_classic() + 
  labs(fill = "Percentage of threatened utilized species\nwith alternatives")+ 
  ggtitle("")+
  theme(legend.position="bottom")

pdf("Fig3_Maps_legend.pdf", 16, 8) 
map_HUMa_ALT80BS_legend
dev.off()


##################################






################################
# calculate phylogenetic signal using the Fritz & Purvis D stat
# https://conbio.onlinelibrary.wiley.com/doi/epdf/10.1111/j.1523-1739.2010.01455.x
################################

# get data
d <- read.table("Phylogeny_NoCon_Henderson_Tree_1_taxa_forPy_CONS-CONT-TRAITSALL-USES_lessCol.txt", header=T, stringsAsFactors = F, sep="\t")

# ensure taxon names are same as in tree
d$Tree_name <- gsub(" ", "_", d$Tree_name)

# get tree
tree_ed <- read.nexus("Phylogeny_NoCon_Henderson_AllTreesCorrected_NEX_TA.nex")
tree_ed2 <- root(tree_ed, node=2489) # root on Calamoideae
#tree_ed3 <- ladderize(tree_ed2)

tree_ed4 <- tree_ed
tree_ed4$edge.length[tree_ed4$edge.length<=0] <- quantile(tree_ed4$edge.length,0.1)*0.1

tree_ed5 <- tree_ed2
tree_ed5$edge.length[tree_ed5$edge.length<=0] <- quantile(tree_ed5$edge.length,0.1)*0.1

# 0: Brownian phylo structure; 1: random
D_HUF <- phylo.d(d, tree_ed4, Tree_name, HumanFood, permut = 1000, rnd.bias=NULL)
D_HUF_r <- phylo.d(d, tree_ed5, Tree_name, HumanFood, permut = 1000, rnd.bias=NULL)
D_HUMa <- phylo.d(d, tree_ed4, Tree_name, Materials, permut = 1000, rnd.bias=NULL)
D_HUMa_r <- phylo.d(d, tree_ed5, Tree_name, Materials, permut = 1000, rnd.bias=NULL)
D_HUMe <- phylo.d(d, tree_ed4, Tree_name, Medicines, permut = 1000, rnd.bias=NULL)
D_HUMe_r <- phylo.d(d, tree_ed5, Tree_name, Medicines, permut = 1000, rnd.bias=NULL)
D_HUC <- phylo.d(d, tree_ed4, Tree_name, SocialUses, permut = 1000, rnd.bias=NULL)
D_HUC_r <- phylo.d(d, tree_ed5, Tree_name, SocialUses, permut = 1000, rnd.bias=NULL)

Use_cat <- as.data.frame(c("HumanFood", "Materials", "Medicine", "Cultural"))
names(Use_cat)[1] <- "Use"
Use_cat$D <- c(D_HUF$DEstimate, D_HUMa$DEstimate, D_HUMe$DEstimate, D_HUC$DEstimate)
Use_cat$pval_BM <- c(D_HUF$Pval0, D_HUMa$Pval0, D_HUMe$Pval0, D_HUC$Pval0)
Use_cat$pval_Random <- c(D_HUF$Pval1, D_HUMa$Pval1, D_HUMe$Pval1, D_HUC$Pval1)
Use_cat$D_r <- c(D_HUF_r$DEstimate, D_HUMa_r$DEstimate, D_HUMe_r$DEstimate, D_HUC_r$DEstimate)
Use_cat$pval_BM_r <- c(D_HUF_r$Pval0, D_HUMa_r$Pval0, D_HUMe_r$Pval0, D_HUC_r$Pval0)
Use_cat$pval_Random_r <- c(D_HUF_r$Pval1, D_HUMa_r$Pval1, D_HUMe_r$Pval1, D_HUC_r$Pval1)

# see Fritz&Purvis paper for interpretation
# result of pval0 (pval_BM) = how significantly different is D from BM phylo (ie from 0) --> <0.05: significantly different from BM phylo, >0.05: not significantly different from BM phylo
# result of pval1 (pval_Random) = how significantly different is D from random (ie from 1) --> <0.05: significantly different from random, >0.05: not significantly different from random
# As long as D < 1 and pval1 < 0.05, there is some degree of phylo signal (if D>1 and pval1 is <0.05, the trait is "overdispersed", i.e. even more dispersed than randomly)
# if pval0 > 0.05, the trait follows BM phylo

# in our case, D is always < 1
# and pvalRandom=pval1 is always <0.05, so D is always significantly different from random, so there is some degree of phylo signal
# in addition, pvalBM=pval0 is always <0.05, so D is always significantly different from BM (although almost not for Materials: D=0.034 and D_r = 0.046, so Materials almost follow BM)

write.table(Use_cat, "Uses_phylogenetic_signal_D.txt", sep="\t")

# is the use with signal also the most documented? Kind of: (careful here some redundant because tree names, but real numbers are similar, see below)
# > sum(d$HumanFood, na.rm = T)
# [1] 368
# > sum(d$Materials, na.rm = T)
# [1] 632
# > sum(d$Medicines, na.rm = T)
# [1] 296
# > sum(d$SocialUses, na.rm = T)
# [1] 216



# calculate global numbers of threatened used species based on each scenario and RLr
d <- read.table("Best_case_Worst_case_table.txt", header=T, stringsAsFactors = F, sep="\t")

d_HUF <- subset(d, d$Human_use_food == "YES")
d_HUMa <- subset(d, d$Human_use_materials == "YES")
d_HUMe <- subset(d, d$Human_use_medicine == "YES")
d_HUC <- subset(d, d$Human_use_cultural == "YES")


Use_cat <- as.data.frame(c("HumanFood", "Materials", "Medicine", "Cultural"))
names(Use_cat)[1] <- "Use"
Use_cat$SpeNum <- c(length(d_HUF[,1]), length(d_HUMa[,1]), length(d_HUMe[,1]), length(d_HUC[,1]))
Use_cat$RLr_Tnum <- c(sum(d_HUF$RLr == "nonLC", na.rm = T), sum(d_HUMa$RLr == "nonLC", na.rm = T), sum(d_HUMe$RLr == "nonLC", na.rm = T), sum(d_HUC$RLr == "nonLC", na.rm = T))
Use_cat$BS_Tnum <- c(sum(d_HUF$RLr_MLpr_BS == "nonLC", na.rm = T), sum(d_HUMa$RLr_MLpr_BS == "nonLC", na.rm = T), sum(d_HUMe$RLr_MLpr_BS == "nonLC", na.rm = T), sum(d_HUC$RLr_MLpr_BS == "nonLC", na.rm = T))
Use_cat$WS_Tnum <- c(sum(d_HUF$RLr_MLpr_WS == "nonLC", na.rm = T), sum(d_HUMa$RLr_MLpr_WS == "nonLC", na.rm = T), sum(d_HUMe$RLr_MLpr_WS == "nonLC", na.rm = T), sum(d_HUC$RLr_MLpr_WS == "nonLC", na.rm = T))
Use_cat$RLr_NTnum <- c(sum(d_HUF$RLr == "LC", na.rm = T), sum(d_HUMa$RLr == "LC", na.rm = T), sum(d_HUMe$RLr == "LC", na.rm = T), sum(d_HUC$RLr == "LC", na.rm = T))
Use_cat$BS_NTnum <- c(sum(d_HUF$RLr_MLpr_BS == "LC", na.rm = T), sum(d_HUMa$RLr_MLpr_BS == "LC", na.rm = T), sum(d_HUMe$RLr_MLpr_BS == "LC", na.rm = T), sum(d_HUC$RLr_MLpr_BS == "LC", na.rm = T))
Use_cat$WS_NTnum <- c(sum(d_HUF$RLr_MLpr_WS == "LC", na.rm = T), sum(d_HUMa$RLr_MLpr_WS == "LC", na.rm = T), sum(d_HUMe$RLr_MLpr_WS == "LC", na.rm = T), sum(d_HUC$RLr_MLpr_WS == "LC", na.rm = T))
Use_cat$RLr_SpeNumKnown <- c(sum(Use_cat$RLr_Tnum[1], Use_cat$RLr_NTnum[1]), sum(Use_cat$RLr_Tnum[2], Use_cat$RLr_NTnum[2]), sum(Use_cat$RLr_Tnum[3], Use_cat$RLr_NTnum[3]), sum(Use_cat$RLr_Tnum[4], Use_cat$RLr_NTnum[4]))
Use_cat$BS_SpeNumKnown <- c(sum(Use_cat$BS_Tnum[1], Use_cat$BS_NTnum[1]), sum(Use_cat$BS_Tnum[2], Use_cat$BS_NTnum[2]), sum(Use_cat$BS_Tnum[3], Use_cat$BS_NTnum[3]), sum(Use_cat$BS_Tnum[4], Use_cat$BS_NTnum[4]))
Use_cat$WS_SpeNumKnown <- c(sum(Use_cat$WS_Tnum[1], Use_cat$WS_NTnum[1]), sum(Use_cat$WS_Tnum[2], Use_cat$WS_NTnum[2]), sum(Use_cat$WS_Tnum[3], Use_cat$WS_NTnum[3]), sum(Use_cat$WS_Tnum[4], Use_cat$WS_NTnum[4]))
Use_cat$RLr_Tper <- c((Use_cat$RLr_Tnum[1]*100/Use_cat$RLr_SpeNumKnown[1]), (Use_cat$RLr_Tnum[2]*100/Use_cat$RLr_SpeNumKnown[2]), (Use_cat$RLr_Tnum[3]*100/Use_cat$RLr_SpeNumKnown[3]), (Use_cat$RLr_Tnum[4]*100/Use_cat$RLr_SpeNumKnown[4]))
Use_cat$BS_Tper <- c((Use_cat$BS_Tnum[1]*100/Use_cat$BS_SpeNumKnown[1]), (Use_cat$BS_Tnum[2]*100/Use_cat$BS_SpeNumKnown[2]), (Use_cat$BS_Tnum[3]*100/Use_cat$BS_SpeNumKnown[3]), (Use_cat$BS_Tnum[4]*100/Use_cat$BS_SpeNumKnown[4]))
Use_cat$WS_Tper <- c((Use_cat$WS_Tnum[1]*100/Use_cat$WS_SpeNumKnown[1]), (Use_cat$WS_Tnum[2]*100/Use_cat$WS_SpeNumKnown[2]), (Use_cat$WS_Tnum[3]*100/Use_cat$WS_SpeNumKnown[3]), (Use_cat$WS_Tnum[4]*100/Use_cat$WS_SpeNumKnown[4]))


write.table(Use_cat, "Uses_global_props.txt", sep="\t")


####################################



#########################
# calculate functional distinctiveness for each species - done once, can reimport it when needed
# stored in Functional_distinctiveness.txt
########################

# format data to fit funrar requirements

# reimport traits table

traits <- read.table("All_species_cons_and_traits_newly_combined_Prec.txt", header=T,sep="\t")

# recode the categorical variables as binary ones
traits$ClimbingYes <- rep(NA, length(traits$Climbing))
#traits$ClimbingNo <- rep(NA, length(traits$Climbing))
traits$AcaulescentYes <- rep(NA, length(traits$Climbing))
#traits$AcaulescentNo <- rep(NA, length(traits$Climbing))
traits$ErectYes <- rep(NA, length(traits$Climbing))
#traits$ErectNo <- rep(NA, length(traits$Climbing))
traits$SolitaryYes <- rep(NA, length(traits$Climbing))
traits$ClusteringYes <- rep(NA, length(traits$Climbing))
traits$UnderstoreyYes <- rep(NA, length(traits$Climbing))
traits$CanopyYes <- rep(NA, length(traits$Climbing))
traits$ConspicuousYes <- rep(NA, length(traits$Climbing))

for (x in 1:length(traits$Climbing)) {
  traits$ClimbingYes[x] <- ifelse(as.numeric(traits$Climbing[x])>0, 1, 0)
  #traits$ClimbingNo[x] <- ifelse(as.numeric(traits$Climbing[x])>0, 0, 1)
  traits$AcaulescentYes[x] <- ifelse(as.numeric(traits$Acaulescent[x])>0, 1, 0)
  #traits$AcaulescentNo[x] <- ifelse(as.numeric(traits$Acaulescent[x])>0, 0, 1)
  traits$ErectYes[x] <- ifelse(as.numeric(traits$Erect[x])>0, 1, 0)
  #traits$ErectNo[x] <- ifelse(as.numeric(traits$Erect[x])>0, 0, 1)
  traits$SolitaryYes[x] <- ifelse(as.numeric(traits$StemSolitary[x]) >0, 1, 0)
  traits$ClusteringYes[x] <- ifelse(as.numeric(traits$StemSolitary[x]) < 1 || as.numeric(traits$StemSolitary[x]) > 1 , 1, 0)
  traits$UnderstoreyYes[x] <- ifelse(as.character(traits$UnderstoreyCanopy[x]) == "understorey" || as.character(traits$UnderstoreyCanopy[x]) == "both", 1, 0)
  traits$CanopyYes[x] <- ifelse(as.character(traits$UnderstoreyCanopy[x]) == "canopy" || as.character(traits$UnderstoreyCanopy[x]) == "both", 1, 0)
  traits$ConspicuousYes[x] <- ifelse(as.character(traits$Conspicuousness[x]) == "conspicuous", 1, 0)
}



# keep only interesting columns
to_keep <- c(names(traits)[11:12], names(traits)[14:15], names(traits)[51], names(traits)[55:57], names(traits)[73:80])

traits2 <- traits[,names(traits) %in% to_keep]
rownames(traits2) <- traits$species
traits <- traits2


# logtr and then rescale all variables by making them range from 0 to 1

traits_sc <- traits
for (i in 1:length(names(traits))) {
  traits_sc[,i] <- log(as.numeric(traits[,i]))
  traits_sc[,i] <- rescale(traits[,i])
}


# Imputte NAs using the weighted average of 11-NN
traits_knn <- knnImputation(traits_sc, k=11, scale = F) # we already scaled ; 11 to avoid ties with categorical vars
rownames(traits_knn) <- rownames(traits)

# for binary variables, round the weighted average so that 0 or 1
traits_knn2 <- traits_knn

for (i in c(5, 9, 10, 11, 12, 13, 14 , 15, 16)) {
  print(names(traits_knn)[i])
  print(which(traits_knn[,i] == 0.5)) # not sure this works, if not, 0.5 will have gone un-noticed and rounded to 0 (hopefully no 0.5 because used an odd k above)
  traits_knn2[,i] <- round(traits_knn[,i])
}

# make a df with binary as factors to use with Gower distance
traits_knn3 <- traits_knn2
for (x in c(5, 9, 10, 11, 12, 13, 14 , 15, 16)) {
  traits_knn3[,x] <- as.factor(traits_knn2[,x])
}


# compute distance matrix (all vars should be numeric) - try euclidean and Gower but will use Gower
trait_dist_e <- compute_dist_matrix(traits_knn2, metric="euclidean")
trait_dist_g <- compute_dist_matrix(traits_knn2, metric="gower")
trait_dist_g2 <- compute_dist_matrix(traits_knn3, metric="gower") # we will go for this one

# make site_species table (with only one site)
site_spe = matrix(1, ncol = nrow(trait_dist_g2))
colnames(site_spe) = colnames(trait_dist_g2)
row.names(site_spe) = c("global")

# compute functional distinctiveness
#FD_e = distinctiveness(site_spe, trait_dist_e)
FD_g = distinctiveness(site_spe, trait_dist_g2)

# export results to be joined later with other data
FD_g2 <- as.data.frame(t(FD_g))
colnames(FD_g2)[1] <- "FD_G"
FD_g2$species <- rownames(FD_g2)

write.csv(FD_g2, "Functional_distinctiveness.txt", row.names = F)

#############################


##############################
# Calculate functional distances between species, for a matrix of 3 traits 
# Will export them to calculate scores 1 and 2
################################

# impute missing data as above, same script
# format data to fit funrar requirements

# reimport traits table

traits <- read.table("All_species_cons_and_traits_newly_combined_Prec.txt", header=T,sep="\t")

# recode the categorical variables as binary ones
traits$ClimbingYes <- rep(NA, length(traits$Climbing))
#traits$ClimbingNo <- rep(NA, length(traits$Climbing))
traits$AcaulescentYes <- rep(NA, length(traits$Climbing))
#traits$AcaulescentNo <- rep(NA, length(traits$Climbing))
traits$ErectYes <- rep(NA, length(traits$Climbing))
#traits$ErectNo <- rep(NA, length(traits$Climbing))
traits$SolitaryYes <- rep(NA, length(traits$Climbing))
traits$ClusteringYes <- rep(NA, length(traits$Climbing))
traits$UnderstoreyYes <- rep(NA, length(traits$Climbing))
traits$CanopyYes <- rep(NA, length(traits$Climbing))
traits$ConspicuousYes <- rep(NA, length(traits$Climbing))

for (x in 1:length(traits$Climbing)) {
  traits$ClimbingYes[x] <- ifelse(as.numeric(traits$Climbing[x])>0, 1, 0)
  #traits$ClimbingNo[x] <- ifelse(as.numeric(traits$Climbing[x])>0, 0, 1)
  traits$AcaulescentYes[x] <- ifelse(as.numeric(traits$Acaulescent[x])>0, 1, 0)
  #traits$AcaulescentNo[x] <- ifelse(as.numeric(traits$Acaulescent[x])>0, 0, 1)
  traits$ErectYes[x] <- ifelse(as.numeric(traits$Erect[x])>0, 1, 0)
  #traits$ErectNo[x] <- ifelse(as.numeric(traits$Erect[x])>0, 0, 1)
  traits$SolitaryYes[x] <- ifelse(as.numeric(traits$StemSolitary[x]) >0, 1, 0)
  traits$ClusteringYes[x] <- ifelse(as.numeric(traits$StemSolitary[x]) < 1 || as.numeric(traits$StemSolitary[x]) > 1 , 1, 0)
  traits$UnderstoreyYes[x] <- ifelse(as.character(traits$UnderstoreyCanopy[x]) == "understorey" || as.character(traits$UnderstoreyCanopy[x]) == "both", 1, 0)
  traits$CanopyYes[x] <- ifelse(as.character(traits$UnderstoreyCanopy[x]) == "canopy" || as.character(traits$UnderstoreyCanopy[x]) == "both", 1, 0)
  traits$ConspicuousYes[x] <- ifelse(as.character(traits$Conspicuousness[x]) == "conspicuous", 1, 0)
}



# keep only interesting columns
traits2 <- traits[,c(11,12,14,15,51,55,56,57,73:80)]
rownames(traits2) <- traits$species
traits <- traits2

# create a new trait stem volume
traits$MaxStemVolume_cm3 <- pi * ((traits$MaxStemDia_cm/2)*(traits$MaxStemDia_cm/2)) * traits$MaxStemHeight_m * 100


# logtr and then rescale all variables by making them range from 0 to 1

traits_sc <- traits
for (i in 1:length(names(traits))) {
  traits_sc[,i] <- log(as.numeric(traits[,i]))
  traits_sc[,i] <- rescale(traits[,i])
}


# Imputte NAs using the weighted average of 11-NN
traits_knn <- knnImputation(traits_sc, k=11, scale = F) # we already scaled ; 11 to avoid ties with categorical vars
rownames(traits_knn) <- rownames(traits)

# for binary variables, round the weighted average so that 0 or 1
traits_knn2 <- traits_knn

for (i in c(5, 9, 10, 11, 12, 13, 14 , 15, 16)) {
  print(names(traits_knn)[i])
  print(which(traits_knn[,i] == 0.5)) # not sure this works, if not, 0.5 will have gone un-noticed and rounded to 0 (hopefully no 0.5 because used an odd k above)
  traits_knn2[,i] <- round(traits_knn[,i])
}

# make a df with binary as factors to use with Gower distance
traits_knn3 <- traits_knn2
for (x in c(5, 9, 10, 11, 12, 13, 14 , 15, 16)) {
  traits_knn3[,x] <- as.factor(traits_knn2[,x])
}

# keep only interesting traits, i.e. the ones in the Camara-Leret et al 2017 paper, armed, and climbing
traits_knn3_2 <- traits_knn3[,c(4, 5, 6, 9, 17)]

# compute a distance matrix from the traits and export it
traits_knn3_2_dist <- compute_dist_matrix(traits_knn3_2, metric="gower")


write.table(traits_knn3_2_dist, "FunDist_5traits.txt", row.names = T, sep="\t")


# Make a pca of the traits (needs first to make the binary numeric again)
traits_knn4 <- traits_knn3_2
for (x in c(2,4)) {
  traits_knn4[,x] <- as.numeric(traits_knn3_2[,x])
}

traits_knn4_pca <- prcomp(traits_knn4)

summary(traits_knn4_pca)

pca_plot <- fviz_pca_biplot(traits_knn4_pca,label="var")


# same but without binary because weird results
# keep only interesting traits, i.e. the ones in the Camara-Leret et al 2017 paper
traits_knn5 <- traits_knn4[,c(1,3,5)]

# compute a distance matrix from the traits and export it
traits_knn5_dist <- compute_dist_matrix(traits_knn5, metric="euclidean")

write.table(traits_knn5_dist, "FunDist_3traits.txt", row.names = T, sep="\t")

# Make a pca of the traits 
traits_knn5_pca <- prcomp(traits_knn5)

summary(traits_knn5_pca)

pca_plot <- fviz_pca_biplot(traits_knn5_pca,label="var")


pdf("PCA_3traits.pdf", 30, 30) 
pca_plot
dev.off()




##################################




################################
# calculate ED
# time consuming, done once, reimport when needed
# stored in Evolutionary_distinctiveness_Average.txt
# and Evolutionary_distinctiveness_allTrees.txt
################################

# import old table with tree name equivalent to species name
d <- read.table("Phylogeny_NoCon_Henderson_Tree_1_taxa_forPy_CONS-CONT-TRAITSALL-USES_lessCol.txt", header=T, stringsAsFactors = F, sep="\t")
tree_names_to_keep <- d$Tree_name

# build future result table with correspondancy tree name, accepted name
d2 <- d[,1:2]

# import all ultrametric trees
trees_ed <- read.tree("Phylogeny_NoCon_Henderson_AllTreesCorrected.new")

# use labels of first tree to build a list of tips to drop because not in trees_names_to_keep
drops <- c()
z <- 0
for (x in 1:length(trees_ed[[1]]$tip.label)) {
  if (gsub("_", " ", trees_ed[[1]]$tip.label[x]) %in% tree_names_to_keep) {
    z = z+1
  } else {
    drops <- c(drops, trees_ed[[1]]$tip.label[x])
  }
}

# in fact useless, tips correspond to the tree names to keep so don't do the following
# Drop the tips in the drops list from all trees.
# trees_ed_d <- lapply(trees_ed,drop.tip,tip=drops$V1)

# check that they are all scaled the same 
# they are not, they reflect incertainty in ages, but I think that's ok because the branch lengths are divided among child tips
# root_height <- rep(0,length(trees_ed_d))
# for (i in 1:length(trees_ed_d)) {
#   root_height[i] <- get.rooted.tree.height(trees_ed_d[[i]])
# }
# min(root_height)
# max(root_height)

# calculate ED for the first tree (will use results to build final result table)
evoldis <- ed.calc(trees_ed[[1]], polytomy.cf="isaac")
ED <- evoldis$spp
names(ED)[2] <- 1 # corresponds to the results for the first tree

# calculate ED for all trees
for (i in 2:length(trees_ed)) {
  evoldis <- ed.calc(trees_ed[[i]], polytomy.cf="isaac")
  names(evoldis$spp)[2] <- i
  ED <- merge(ED, evoldis$spp,by="species")
}

# checking some species randomly, it does not seem that the first trees have much different values from the rest, so it seems that there is no need to remove a burnin fraction
# But it may be good to ask Soren to confirm
plot(as.numeric(ED[1,2:1000]))
plot(as.numeric(ED[10,2:1000]))
plot(as.numeric(ED[500,2:1000]))

# calculate average ED per species
# do it also for only the trees 251-1000 just in case later need to remove burnin

ED$AverageED <- rep(0,length(ED$species))
ED$AverageED_noBurnin <- rep(0,length(ED$species))
for (i in 1:length(ED$species)){
  ED$AverageED[i] <- mean(as.numeric(ED[i,2:1001]))
  ED$AverageED_noBurnin[i] <- mean(as.numeric(ED[i,252:1001]))
}

names(ED)[1] <- "Tree_name"
ED$Tree_name <- gsub("_", " ", ED$Tree_name)

ED2 <- ED[,c(1,1002, 1003)]

# export results to be joined later with other data
# as well as a copy of ED with all values

d3 <- full_join(d2, ED, by="Tree_name")
d4 <- full_join(d2, ED2, by="Tree_name")

#write.table(d3, "Evolutionary_distinctiveness_allTrees.txt", row.names = F, sep="\t")
write.table(d4, "Evolutionary_distinctiveness_Average.txt", row.names = F, sep="\t")

##############################################



###########################################
# Calculate phylogenetic distances between species
# Will export them to calculate scores 1 and 2
###########################################

# calculate phylogenetic distances among species and export them
trees_ed <- read.tree("Phylogeny_NoCon_Henderson_AllTreesCorrected.new") # trees are rooted already



samp <- sample(251:1000,100)
# samp
# [1]  536  353  389  996  556  846  272  759  504  657  398  737  889  775  339  520  766  608  558  785  252  948  823 1000  342  400  296  957  372  762  643  945  658  738  275
# [36]  451  744  641  523  950  743  467  406  613  703  998  524  697  490  538  529  411  367  733  568  919  764  638  413  758  753  644  391  989  729  390  481  653  682  813
# [71]  838  699  603  534  904  437  355  942  714  600  507  746  254  991  393  257  325  936  465  370  868  463  260  983  262  754  914  545  634  514

pdm_l <- vector(mode = "list", length = 100)
z <- 0

for (i in 1:length(samp)) { # only use 100 post burnin trees (in case burnin trees were included; unlikely though) because all trees is too memory consuming
  z <- z+1
  print(z)
  pdm <- cophenetic.phylo(trees_ed[[samp[i]]])
  pdm_l[[z]] <- pdm
}


# Tried this to see if this allowed to use 750 trees but too slow
# t <- as.data.frame(matrix(data=NA, nrow=750, ncol=length(trees_ed[[1]]$tip.label)))
# names(t) <- trees_ed[[1]]$tip.label
# l <- vector(mode="list", length=length(trees_ed[[1]]$tip.label))
# for (x in 1:length(l)) {
#   l[[x]] <- t
# }
# 
# h <- hash(keys=trees_ed[[1]]$tip.label, values=l)
# 
# for (i in 251:255) { # only use post burnin trees, in case burnin trees were included (unlikely though)
#   z <- z+1
#   print(z)
#   pdm <- cophenetic.phylo(trees_ed[[i]])
#   for (s1 in 1:length(rownames(pdm))) {
#     for (s2 in 1:length(rownames(pdm))) {
#       n1 <- rownames(pdm)[s1]
#       n2 <- rownames(pdm)[s2]
#       h[[n1]][z,which(names(h[[n1]]) == n2)] <- pdm[n1,n2]
#     }
#   }
# }




pdm_mean <- as.data.frame(matrix(data=NA,nrow=length(trees_ed[[1]]$tip.label), ncol=length(trees_ed[[1]]$tip.label)))
rownames(pdm_mean) <- trees_ed[[1]]$tip.label
names(pdm_mean) <- trees_ed[[1]]$tip.label

for (s1 in 1:length(trees_ed[[1]]$tip.label)) {
  print(s1)
  for (s2 in 1:length(trees_ed[[1]]$tip.label)) {
    n1 <- trees_ed[[1]]$tip.label[s1]
    n2 <- trees_ed[[1]]$tip.label[s2]
    v <- c()
    for (x in 1:length(pdm_l)) {
      di <- pdm_l[[x]][n1,n2]
      if (di<0) {
        print(paste(n1,n2,di))
      }
      v <- c(v, di)
    }
    pdm_mean[n1,n2] <- mean(v)
  }
}




write.table(pdm_mean, "PhyloDist_Average100Trees.txt", row.names = T, sep="\t")

###########################################






