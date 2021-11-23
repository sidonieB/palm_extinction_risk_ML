####################################
# Prepare data for machine learning
####################################

library(lattice)
library(ggplot2)
library(caret)
library(HistogramTools)
library(grid)
library(gridExtra)
library(mosaic)
library(proxy)
library(plyr)
library(scales)
library(UBL)

### import and format data

all <- read.table("MLspecies_all_preds_and_cons.txt", sep="\t", header=T)

# make sure all predictor variables are numeric
all2 <- all
all2$Nbe_unique_occ. <- as.numeric(all2$Nbe_unique_occ.)
all2$Nbe_subPop <- as.numeric(all2$Nbe_subPop)
all2$Nbe_loc <- as.numeric(all2$Nbe_loc)
all2$AOO2km <- as.numeric(all2$AOO2km)
all2$ECO <- as.numeric(all2$ECO)
all2$TDWG3 <- as.numeric(all2$TDWG3)

# change predictor names so that nicer

names(all2) <- gsub("Nbe_unique_occ.","Nocc",names(all2))
names(all2) <- gsub("Nbe_subPop","NsubPop",names(all2))
names(all2) <- gsub("Nbe_loc","Nloc",names(all2))
names(all2) <- gsub("AOO2km","AOO",names(all2))
names(all2) <- gsub("bio4","Tseas",names(all2))
names(all2) <- gsub("bio15","Pseas",names(all2))
names(all2) <- gsub("Tflex_mean","Tflex",names(all2))
names(all2) <- gsub("Pflex_mean","Pflex",names(all2))
names(all2) <- gsub("EOO_rCAT_conR2","EOO",names(all2))

# add a category group corresponding to assessed/non-assessed (effectively same as splitting by $preds but nicer labels for later)
all2$Subset <- rep("ASSESSED", length(all2$accepted_name))
for (s in 1:length(all2$Subset)) {
  if (all2$IUCN_TS_sum[s] =="NE") {
    all2$Subset[s] <- "NOT_ASSESSED"
  }
}

all2_sub <- split(all2, all2$Subset)


### Assess the geographic and phylogenetic representation of the different subsets
length(all2$accepted_name)

traits <-read.table("Continental_distribution.txt", 
                    header=T, row.names=NULL, sep="\t", stringsAsFactors = F, na.strings = F)

cons <- all2

head(traits)
head(cons)
cons_2 <- cons[,c(1, 21)] 
names(cons_2)[1] <- "species"
df <- full_join(traits, cons_2, by="species")

for (x in 1:length(df$species)){
  if(is.na(df$Subset[x])) {
    df$Subset[x] <- "NO_ML"
  } else if (df$Subset[x] == "NOT_ASSESSED") {
    df$Subset[x] <- "PREDICTED"
  }
}


# remove the species that will not be included in the ML
# because we care about the difference between what we use to train the model and what we will predict
# create a column genus, which will be our approximation for phylogenetic position

df_onlyML <- subset(df, df$Subset != "NO_ML")

df_onlyML$genus <- gsub(" .*", "", df_onlyML$species)


# continent bias
pC <- ggplot(data=df_onlyML, aes(x=Continent, fill=Subset)) + 
  geom_bar(position = "dodge")+ 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Continental representation")

# taxonomy bias
pT <- ggplot(data=df_onlyML, aes(x=genus, fill=Subset)) + 
  geom_bar(position = "dodge") + 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Taxonomic representation")


# export

pdf("Geographic_and_taxonomic_biases.pdf", 24, 6) 
grid.arrange(pC, pT, nrow = 1)
dev.off()

setEPS()
postscript("Geographic_and_taxonomic_biases.eps", height = 6, width = 24)
grid.arrange(pC, pT, nrow = 1)
dev.off()




### Prepare the data subsets

# The assessed species have to be subsampled so that they don't have strong geographic and taxonomic biases compared to 
# the species to predict
# Ideally, the proportion of species from continent x in the training set should be <= the proportion of species from this 
# continent in the species to predict, and same for genus
# We will then split the training and test set based on this new subset of assessed species or on the
# old one.

# Combine genus_continent as the category to equalize

df_onlyML$GC <- paste(df_onlyML$genus, df_onlyML$Continent, sep="_")

predicted <- subset(df_onlyML, df_onlyML$Subset == "PREDICTED")
assessed <- subset(df_onlyML, df_onlyML$Subset == "ASSESSED")

p_sum <- data.frame(table(predicted$GC))
p_sum$prop <- p_sum$Freq/length(predicted$species)

a_sum <- data.frame(table(assessed$GC))
a_sum$prop <- a_sum$Freq/length(assessed$species)

a_sum$p_res <- rep(NA, length(a_sum$Var1))

# calculate the proportion of the current sample that has to be sampled so that its new proportion equals the target proportion
# a = 32, p=12 so r= p/a
# if r < 1 perform random undersampling
# if r > 1, keep all
# plot again and see

for (x in 1:length(a_sum$Var1)) {
  if (as.character(a_sum$Var1[x]) %in% p_sum$Var1){
  a <- a_sum$prop[x]
  p <- p_sum$prop[which(p_sum$Var1 == as.character(a_sum$Var1[x]))]
  r <- p/a
  a_sum$p_res[x] <- r
  } else {
    if(a_sum$Freq[x] <= 3){
      a_sum$p_res[x] <- 1 # keep all species if no more than 3
    } else {
      a_sum$p_res[x] <- 0.25 # keep a fourth if more than 3 (range from 4 to 13 so will keep ca. 1 to 3 but not exactly: will depend on final numbers)
    }
    }
}

# prepare input for resampling function
# keep only classes to undersample (i.e. classes with p_res <1)
a_sum2 <- subset(a_sum, a_sum$p_res < 1)

a_sum_l <- vector(mode="list", length=length(a_sum2$Var1))
for (x in 1:length(a_sum2$Var1)) {
  a_sum_l[[x]] <- a_sum2$p_res[x]
  names(a_sum_l)[[x]] <- as.character(a_sum2$Var1[x])
}




# undersample the assessed species to reach the new proportions
assessed_us <- RandUnderClassif(GC~., assessed, C.perc = a_sum_l, repl = FALSE)


# plot the new assessed dataset in comparison to the predicted
df_onlyML_us <- rbind(predicted, assessed_us)


pCus <- ggplot(data=df_onlyML_us, aes(x=Continent, fill=Subset)) + 
  geom_bar(position = "dodge")+ 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Continental representation")

# taxonomy bias
pTus <- ggplot(data=df_onlyML_us, aes(x=genus, fill=Subset)) + 
  geom_bar(position = "dodge") + 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Taxonomic representation")


# export

pdf("Geographic_and_taxonomic_biases_BeforeAfterUS.pdf", 24, 12) 
grid.arrange(pC, pT, pCus, pTus, nrow = 2)
dev.off()

setEPS()
postscript("Geographic_and_taxonomic_biases_BeforeAfterUS.eps", height = 12, width = 24)
grid.arrange(pC, pT, pCus, pTus, nrow = 2)
dev.off()


# much better but so much lost, try to do the undersampling only on African plants so that all species from 
# other continents are left

a_sum$p_resAf <- rep(NA, length(a_sum$Var1))

for (x in 1:length(a_sum$Var1)) {
  if (length(grep("Africa", as.character(a_sum$Var1[x]))) > 0 ) {
  if (as.character(a_sum$Var1[x]) %in% p_sum$Var1){
    a <- a_sum$prop[x]
    p <- p_sum$prop[which(p_sum$Var1 == as.character(a_sum$Var1[x]))]
    r <- p/a
    a_sum$p_resAf[x] <- r
  } else {
    if(a_sum$Freq[x] <= 10){
      a_sum$p_resAf[x] <- 1 # keep all species if no more than 10
    } else {
      a_sum$p_resAf[x] <- 0.5 # keep a half if more than 10 (only one with 18 in Africa)
    }
  }
  } else {
    a_sum$p_resAf[x] <- 1
  }
  if (a_sum$p_resAf[x]*a_sum$Freq[x] < 1) { # do not resample if would bring species number to less than 1 # effectively will resample only Dypsis
    a_sum$p_resAf[x] <- 1
  }
  
  
}

# prepare input for resampling function
# keep only classes to undersample (i.e. classes with p_res <1)
a_sum2_Af <- subset(a_sum, a_sum$p_resAf < 1)

a_sum_l_Af <- vector(mode="list", length=length(a_sum2_Af$Var1))

for (x in 1:length(a_sum2_Af$Var1)) {
  a_sum_l_Af[[x]] <- a_sum2_Af$p_res[x]
  names(a_sum_l_Af)[[x]] <- as.character(a_sum2_Af$Var1[x])
}


# undersample the assessed species to reach the new proportions
assessed_usAf <- RandUnderClassif(GC~., assessed, C.perc = a_sum_l_Af, repl = FALSE)
# seed lacking here! Downstream reproducible if use same dataset splits, 
# i.e. assessed_usAf contains species in all2_sub[[1]] that are not in Unused_assessed.txt


# plot the new assessed dataset in comparison to the predicted
df_onlyML_usAf <- rbind(predicted, assessed_usAf)


pCus <- ggplot(data=df_onlyML_usAf, aes(x=Continent, fill=Subset)) + 
  geom_bar(position = "dodge")+ 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Continental representation")

# taxonomy bias
pTus <- ggplot(data=df_onlyML_usAf, aes(x=genus, fill=Subset)) + 
  geom_bar(position = "dodge") + 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Taxonomic representation")


# export

pdf("Geographic_and_taxonomic_biases_BeforeAfterUSAf.pdf", 24, 12) 
grid.arrange(pC, pT, pCus, pTus, nrow = 2)
dev.off()

pdf("Taxonomic_biases_BeforeAfterUSAf.pdf", 48, 12) 
grid.arrange(pT, pTus, nrow = 2)
dev.off()


setEPS()
postscript("Geographic_and_taxonomic_biases_BeforeAfterUSAf.eps", height = 12, width = 24)
grid.arrange(pC, pT, pCus, pTus, nrow = 2)
dev.off()

setEPS()
postscript("Ttaxonomic_biases_BeforeAfterUSAf.eps", height = 12, width = 48)
grid.arrange(pT, pTus, nrow = 2)
dev.off()


# The resampling of only African species is well balanced.
# Now our set of assessed species had 300 species, in assessed_usAf
# rebuild our global dataset to keep only these species in the TRAIN subset
# export the unused assessed species for use in downstream plots

assessed_usAf_all_data <- all2[which(all2$accepted_name %in% assessed_usAf$species),]
unused <- all2_sub[[1]][-which(all2_sub[[1]]$accepted_name %in% assessed_usAf$species),]
write.table(unused, "Unused_assessed.txt", row.names = F, sep="\t")

all2_usAf <- rbind(assessed_usAf_all_data, all2_sub[[2]])

all2_usAf_sub <- split(all2_usAf, all2_usAf$Subset) # kind of unnecessary but naming easier later


# split the assessed set randomly in a training set (75% of the assessed species = 225 sp) and a test set (25% of the assessed species = 75 sp)
# the training set will be used to optimize the model
# the test set will be used to test the performance of the final model on data that were not used to optimize the model
# ideally the training, test and to predict sets would have similar distributions for each of the predictors
# this would help to increase the model performance on the test set 
# but also to increase chances that the performance on the test set will be truely indicative of the performance on the set to predict
# which cannot be known (except later if some manual assessments are done on the predicted species)

# We know the test set has to be made of 75 species to leave ca. 75% of the assessed species for the training set, 
# so start building the test set with a random sample of 15 taken in the assessed set
# Then increment the test set with the most dissimilar samples that can be found in the assessed until reach 75 samples
# https://topepo.github.io/caret/data-splitting.html

# copy the data to keep only predictor and IUCN class column so that dissimilarity is not biased by irrelevant columns
drops <- c(names(all2_usAf)[1:2], names(all2_usAf)[16:17], names(all2_usAf)[19:21])
assessed_simple <- all2_usAf_sub[[1]][,!(names(all2_usAf_sub[[1]]) %in% drops)]

# create a normalized dataset so that some predictors don't have a higher impact when calculating the distances just because they range to higher values
# rescale to min 0 and max 1

assessed_simple_sc <- assessed_simple

for (i in 1:(length(names(assessed_simple))-1)) {
assessed_simple_sc[,i] <- rescale(assessed_simple[,i])
}

# check that it worked
for (i in 1:(length(names(assessed_simple))-1)) {
  print(names(assessed_simple)[i])
  print(min(assessed_simple[,i], na.rm=T))
  print(max(assessed_simple[,i], na.rm=T))
  print(min(assessed_simple_sc[,i], na.rm=T))
  print(max(assessed_simple_sc[,i], na.rm=T))
}


# Draw a random sample of 15 species from the simplified assessed set
set.seed(501)
startSet <- sample(1:length(assessed_simple[,1]), 15) # index of the samples, not the samples themselves
start <- assessed_simple[startSet,] # the 15 samples
# same from the scaled set (will compare the resulting sets)
start_sc <- assessed_simple_sc[startSet,]
# define the remaining species from the assessed set as a pool to get the remaining 60 samples from
samplePool <- assessed_simple[-startSet,]
samplePool_sc <- assessed_simple_sc[-startSet,]
# build the whole test set comprising the first 15 samples + 60 other samples as dissimilar as possible
# do it with both methods minDiss and sumDiss and will compare
# as well as using original and rescaled data
test_set_idx_m <- maxDissim(start, samplePool, n = 75, obj = minDiss)
test_set_idx_s <- maxDissim(start, samplePool, n = 75, obj = sumDiss)
test_set_idx_m_sc <- maxDissim(start_sc, samplePool_sc, n = 75, obj = minDiss)
test_set_idx_s_sc <- maxDissim(start_sc, samplePool_sc, n = 75, obj = sumDiss)


# Also generate a test set by drawing 75 samples among the assessed species, completely randomly
# we will compare it to the other sets, and decide what to use depending on what test set represents better the
# non assessed species, in terms of predictors
set.seed(501)
test_set_idx_r <- sample(1:length(assessed_simple[,1]), 75)

# generate the corresponding sets
test_set_m <-  all2_usAf_sub[[1]][test_set_idx_m,]
test_set_s <-  all2_usAf_sub[[1]][test_set_idx_s,]
test_set_m_sc <-  all2_usAf_sub[[1]][test_set_idx_m_sc,] # the data are not scaled, they are just chosen after scaling
test_set_s_sc <-  all2_usAf_sub[[1]][test_set_idx_s_sc,]
test_set_r <-  all2_usAf_sub[[1]][test_set_idx_r,]
train_set_m <- all2_usAf_sub[[1]][-test_set_idx_m,]
train_set_s <- all2_usAf_sub[[1]][-test_set_idx_s,]
train_set_m_sc <- all2_usAf_sub[[1]][-test_set_idx_m_sc,]
train_set_s_sc <- all2_usAf_sub[[1]][-test_set_idx_s_sc,]
train_set_r <- all2_usAf_sub[[1]][-test_set_idx_r,]
pred_set <- all2_usAf_sub[[2]]

# change the subset column of test and train to get them differentiated in histograms
test_set_m$Subset <- rep("TEST", length(test_set_m$Subset))
test_set_s$Subset <- rep("TEST", length(test_set_s$Subset))
train_set_m$Subset <- rep("TRAIN", length(train_set_m$Subset))
train_set_s$Subset <- rep("TRAIN", length(train_set_s$Subset))
test_set_m_sc$Subset <- rep("TEST", length(test_set_m_sc$Subset))
test_set_s_sc$Subset <- rep("TEST", length(test_set_s_sc$Subset))
train_set_m_sc$Subset <- rep("TRAIN", length(train_set_m_sc$Subset))
train_set_s_sc$Subset <- rep("TRAIN", length(train_set_s_sc$Subset))
test_set_r$Subset <- rep("TEST", length(test_set_r$Subset))
train_set_r$Subset <- rep("TRAIN", length(train_set_r$Subset))

# join all subsets
assessed_split_m <- rbind(test_set_m,train_set_m)
assessed_split_s <- rbind(test_set_s,train_set_s)
assessed_split_r <- rbind(test_set_r,train_set_r)
all3_m <- rbind(assessed_split_m,pred_set)
all3_s <- rbind(assessed_split_s,pred_set)
all3_r <- rbind(assessed_split_r,pred_set)

assessed_split_m_sc <- rbind(test_set_m_sc,train_set_m_sc)
assessed_split_s_sc <- rbind(test_set_s_sc,train_set_s_sc)
all3_m_sc <- rbind(assessed_split_m_sc,pred_set)
all3_s_sc <- rbind(assessed_split_s_sc,pred_set)

# export in case need later
write.table(all3_m, "Data_split_minDiss.txt", row.names = FALSE, sep="\t") 
write.table(all3_s, "Data_split_sumDiss.txt", row.names = FALSE, sep="\t")
write.table(all3_r, "Data_split_random.txt", row.names = FALSE, sep="\t")
write.table(all3_m_sc, "Data_split_minDiss_scaled.txt", row.names = FALSE, sep="\t") 
write.table(all3_s_sc, "Data_split_sumDiss_scaled.txt", row.names = FALSE, sep="\t")

# split by subset to plot
all3_sub_m <- split(all3_m, all3_m$Subset)
all3_sub_s <- split(all3_s, all3_s$Subset)
all3_sub_r <- split(all3_r, all3_r$Subset)

all3_sub_m_sc <- split(all3_m_sc, all3_m_sc$Subset)
all3_sub_s_sc <- split(all3_s_sc, all3_s_sc$Subset)


### Assess the adequacy of the splitting

# Compare the predictor distributions of the test data with these of the training and non-assessed
# Using histogram intersections (HI) and distributions
# we want high HI and distributions of similar shapes, especially for NA (not-assessed) and test, and NA and train

# load functions

plot_HI = function (data, data_sub, c) {
  
  min_b <- min(data[,c], na.rm = T)
  max_b <- max(data[,c], na.rm = T)
  bin_num <- 20
  m <- ceiling((max_b-min_b)/bin_num)
  
  min_b <- min_b - min_b %% m
  if ((max_b %% m) != 0) {
    max_b <- max_b + m - max_b %% m
  }
  
  if (max_b == 1) {
    bin_w <- (max_b - min_b) / bin_num
  } else {
    bin_w <- m
  }
  
  h1 <- hist(data_sub[[1]][,c], breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # pred
  h2 <- hist(c(data_sub[[2]][,c], data_sub[[3]][,c]), breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # test+train
  h3 <- hist(data_sub[[3]][,c], breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # train
  h4 <- hist(data_sub[[2]][,c], breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # test
  hi_pras <- round(intersect.dist(h2, h1), 2)
  hi_prtr <- round(intersect.dist(h3, h1), 2)
  hi_prte <- round(intersect.dist(h4, h1), 2)
  hi_trte <- round(intersect.dist(h4, h3), 2)
  
  HI_pras <- paste("HI_NA-A=", hi_pras, sep="")
  HI_prtr <- paste("HI_NA-train=", hi_prtr, sep="")
  HI_prte <- paste("HI_NA-test=", hi_prte, sep="")
  HI_trte <- paste("HI_train-test=", hi_trte, sep="")
  HI <- paste(HI_pras, HI_prtr, HI_prte, HI_trte, sep="\n")
  grob <- grobTree(textGrob(HI, x=0.3,  y=0.9, hjust=0,
                            gp=gpar(col="black", fontsize=8)))
  
  ggplot(data, aes(x=data[,c], fill=data[,21], color=data[,21])) +
    geom_histogram(alpha=0.7, position="identity", breaks=seq(from=min_b,to=max_b,by=bin_w), na.rm = T) + 
    theme_classic() +
    theme(legend.position="top") +
    scale_color_manual(values=c("grey30", "grey30", "grey30"))+
    scale_fill_manual(values=c("grey100", "grey30", "grey70")) +
    labs(fill = "", color="", x=names(data)[c], y = "Count") + 
    annotation_custom(grob)
  
}

# same but remove largest outliers (> third quartile) to be able to see better some of the distributions

plot_HI_zoom = function (data, data_sub, c) {
  
  #print(c)
  min_b <- min(data[,c], na.rm = T)
  max_b <- favstats(data[,c],data=data)[[4]] # Q3
  bin_num <- 20
  m <- ceiling((max_b-min_b)/bin_num)
  
  min_b <- min_b - min_b %% m
  if ((max_b %% m) != 0) {
    max_b <- max_b + m - max_b %% m
  }
  
  if (all.equal(max_b, 1) == TRUE) {
    bin_w <- (max_b - min_b) / bin_num
  } else {
    bin_w <- m
  }
  
  # need to feed only the values inside the range that we want to see
  h1 <- hist(data_sub[[1]][,c][data_sub[[1]][,c] >= min_b & data_sub[[1]][,c] <= max_b], breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # pred
  h2 <- hist(c(data_sub[[2]][,c][data_sub[[2]][,c] >= min_b & data_sub[[2]][,c] <= max_b], data_sub[[3]][,c][data_sub[[3]][,c] >= min_b & data_sub[[3]][,c] <= max_b]), breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # test+train
  h3 <- hist(data_sub[[3]][,c][data_sub[[3]][,c] >= min_b & data_sub[[3]][,c] <= max_b], breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # train
  h4 <- hist(data_sub[[2]][,c][data_sub[[2]][,c] >= min_b & data_sub[[2]][,c] <= max_b], breaks=seq(from=min_b,to=max_b,by=bin_w), plot=F) # test
  hi_pras <- round(intersect.dist(h2, h1), 2)
  hi_prtr <- round(intersect.dist(h3, h1), 2)
  hi_prte <- round(intersect.dist(h4, h1), 2)
  hi_trte <- round(intersect.dist(h4, h3), 2)
  
  HI_pras <- paste("HI_NA-A=", hi_pras, sep="")
  HI_prtr <- paste("HI_NA-train=", hi_prtr, sep="")
  HI_prte <- paste("HI_NA-test=", hi_prte, sep="")
  HI_trte <- paste("HI_train-test=", hi_trte, sep="")
  HI <- paste(HI_pras, HI_prtr, HI_prte, HI_trte, sep="\n")
  grob <- grobTree(textGrob(HI, x=0.3,  y=0.9, hjust=0,
                            gp=gpar(col="black", fontsize=8)))
  
  
  df <- data[data[,c] >= min_b & data[,c] <= max_b,] # may not be necessary, ggplot seems robust to values out of the breaks limits
  
  ggplot(df, aes(x=df[,c], fill=df[,21], color=df[,21])) +
    geom_histogram(alpha=0.7, position="identity", breaks=seq(from=min_b,to=max_b,by=bin_w), na.rm = T) + 
    theme_classic() +
    theme(legend.position="top") +
    scale_color_manual(values=c("grey30", "grey30", "grey30"))+
    scale_fill_manual(values=c("grey100", "grey30", "grey70")) +
    labs(fill = "", color="", x=names(data)[c], y = "Count") + 
    annotation_custom(grob)
  
}


# run the functions on the data using the m, s and r test and train sets
# all data, or zooming excluding the third quartile

plots_HI <- lapply(c(3:15), plot_HI, data = all3_m, data_sub=all3_sub_m)

pdf("Dataset_shift_minDiss.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI, nrow=3, ncol=5))

dev.off()


plots_HI <- lapply(c(3:15), plot_HI, data = all3_s, data_sub=all3_sub_s)

pdf("Dataset_shift_sumDiss.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI, nrow=3, ncol=5))

dev.off()


plots_HI <- lapply(c(3:15), plot_HI, data = all3_r, data_sub=all3_sub_r)

pdf("Dataset_shift_random.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI, nrow=3, ncol=5))

dev.off()


plots_HI_z <- lapply(c(3:15), plot_HI_zoom, data = all3_m, data_sub=all3_sub_m)

pdf("Dataset_shift_minDiss_zoom.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI_z, nrow=3, ncol=5))

dev.off()


plots_HI_z <- lapply(c(3:15), plot_HI_zoom, data = all3_s, data_sub=all3_sub_s)

pdf("Dataset_shift_sumDiss_zoom.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI_z, nrow=3, ncol=5))

dev.off()

plots_HI_z <- lapply(c(3:15), plot_HI_zoom, data = all3_r, data_sub=all3_sub_r)

pdf("Dataset_shift_random_zoom.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI_z, nrow=3, ncol=5))

dev.off()


# same on sets obtained after scaling


plots_HI <- lapply(c(3:15), plot_HI, data = all3_m_sc, data_sub=all3_sub_m_sc)

pdf("Dataset_shift_minDiss_sc.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI, nrow=3, ncol=5))

dev.off()


plots_HI <- lapply(c(3:15), plot_HI, data = all3_s_sc, data_sub=all3_sub_s_sc)

pdf("Dataset_shift_sumDiss_sc.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI, nrow=3, ncol=5))

dev.off()


plots_HI_z <- lapply(c(3:15), plot_HI_zoom, data = all3_m_sc, data_sub=all3_sub_m_sc)

pdf("Dataset_shift_minDiss_zoom_sc.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI_z, nrow=3, ncol=5))

dev.off()


plots_HI_z <- lapply(c(3:15), plot_HI_zoom, data = all3_s_sc, data_sub=all3_sub_s_sc)

pdf("Dataset_shift_sumDiss_zoom_sc.pdf", 20, 12)

do.call(grid.arrange, c(plots_HI_z, nrow=3, ncol=5))

dev.off()


# Only very few differences (of <= 0.03 units) between results using test set sampled based on scaled data or not
# Will use the test set sampled ased on scaled data

# When comparing the distribtions and HI obtained with minDiss or sumDiss, we basically obtain the same
# In a few cases it is slightly different, but only by 0.01
# In most cases, the differences are because HI_train-test is higher with sumDiss (by 0.04 at most, mostly by 0.01), or HI_NA-test is higher with minDiss (by 0.01)
# Because it is more important to have a high HI_NA-test than a high HI_train-test, and because the minDiss test distribution visually look more representative of the non-assessed, we will use the minDiss test set for the ML

# The difference between minDiss and random is also weak. It seems that in some cases one of the test sets covers extreme parts of the predictor
# distribution better than the other, but there is no consistency as per which test set does it. 
# Check the balance of both test sets in conservation status:
summary(test_set_m_sc$IUCN_TS_sum)
#LC    NE nonLC 
#33     0    42
summary(test_set_r$IUCN_TS_sum)
#LC    NE nonLC 
#36     0    39
# The random one is better balanced (although both are well balanced)
# but the m would give more opportunities to study the performance at classifying nonLC, while still giving a lot opportunities
# to emasure performance at classifying LC
# Check the balance of the remaining training sets in both cases:
summary(train_set_m_sc$IUCN_TS_sum)
#LC    NE nonLC 
#97     0   128 
summary(train_set_r$IUCN_TS_sum)
#LC    NE nonLC 
#94     0   131
# They are both almost identical so the choice of the test set won't affect the balance of the training set regarding the 
# conservation status (it will be unbalanced anyway - we will deal with this during the ML)

# Check the balance of the train and test sets in terms of continents and genera

traits <-read.table("Continental_distribution.txt", 
                    header=T, row.names=NULL, sep="\t", stringsAsFactors = F, na.strings = F)

cons_m <- all3_m_sc
cons_r <- all3_r


cons_m2 <- cons_m[,c(1, 21)] 
cons_r2 <- cons_r[,c(1, 21)]

names(cons_m2)[1] <- "species"
names(cons_r2)[1] <- "species"

df_m <- full_join(traits, cons_m2, by="species")
df_r <- full_join(traits, cons_r2, by="species")

for (x in 1:length(df_m$species)){
  if(is.na(df_m$Subset[x])) {
    df_m$Subset[x] <- "NO_ML"
  } else if (df_m$Subset[x] == "NOT_ASSESSED") {
    df_m$Subset[x] <- "PREDICTED"
  }
}

for (x in 1:length(df_r$species)){
  if(is.na(df_r$Subset[x])) {
    df_r$Subset[x] <- "NO_ML"
  } else if (df_r$Subset[x] == "NOT_ASSESSED") {
    df_r$Subset[x] <- "PREDICTED"
  }
}


# remove the species that will not be included in the ML
# because we care about the difference between what we use to train and test the model and what we will predict
# create a column genus, which will be our approximation for phylogenetic position

df_onlyML_r <- subset(df_r, df_r$Subset != "NO_ML")
df_onlyML_m <- subset(df_m, df_m$Subset != "NO_ML")

df_onlyML_r$genus <- gsub(" .*", "", df_onlyML_r$species)
df_onlyML_m$genus <- gsub(" .*", "", df_onlyML_m$species)

df_onlyML_r_test <- subset(df_onlyML_r, df_onlyML_r$Subset == "TEST")
df_onlyML_r_train <- subset(df_onlyML_r, df_onlyML_r$Subset == "TRAIN")
df_onlyML_m_test <- subset(df_onlyML_m, df_onlyML_m$Subset == "TEST")
df_onlyML_m_train <- subset(df_onlyML_m, df_onlyML_m$Subset == "TRAIN")


summary(as.factor(df_onlyML_m_test$Continent))
#AfricaWAsia     Americas EAsiaPacific 
#17           25           33 
summary(as.factor(df_onlyML_r_test$Continent))
#AfricaWAsia     Americas EAsiaPacific 
#20           26           29 

summary(as.factor(df_onlyML_m_test$genus)) # see plots below
summary(as.factor(df_onlyML_r_test$genus)) # see plots below

length(levels(as.factor(df_onlyML_m_test$genus)))
# 28
length(levels(as.factor(df_onlyML_r_test$genus)))
# 37

# --> the minDiss test set has less genera, both are similarly balanced in terms of continents

summary(as.factor(df_onlyML_m_train$Continent))
#AfricaWAsia     Americas EAsiaPacific 
#57           68          100 
summary(as.factor(df_onlyML_r_train$Continent))
#AfricaWAsia     Americas EAsiaPacific 
#54           67          104

summary(as.factor(df_onlyML_m_train$genus)) # see plots below
summary(as.factor(df_onlyML_r_train$genus)) # see plots below

length(levels(as.factor(df_onlyML_m_train$genus)))
# 77
length(levels(as.factor(df_onlyML_r_train$genus)))
# 71

# --> the minDiss train set has more genera, both are similarly balanced in terms of continents
# it is more important to have as many genera as possible in the training than in the test set
# so we will use the minDiss test (and train) sets


# continent bias
pCm <- ggplot(data=df_onlyML_m, aes(x=Continent, fill=Subset)) + 
  geom_bar(position = "dodge")+ 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Continental representation, minDiss test set")

pCr <- ggplot(data=df_onlyML_r, aes(x=Continent, fill=Subset)) + 
  geom_bar(position = "dodge")+ 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Continental representation, random test set")

# taxonomy bias
pTm <- ggplot(data=df_onlyML_m, aes(x=genus, fill=Subset)) + 
  geom_bar(position = "dodge") + 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Taxonomic representation, minDIss test set")

pTr <- ggplot(data=df_onlyML_r, aes(x=genus, fill=Subset)) + 
  geom_bar(position = "dodge") + 
  scale_fill_grey() +
  theme_classic() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25))+
  labs(title="Taxonomic representation, random test set")


# export

pdf("Geographic_and_taxonomic_biases_random_vs_minDiss_test_sets.pdf", 24, 12) 
grid.arrange(pCm, pTm, pCr, pTr, nrow = 2)
dev.off()

setEPS()
postscript("Geographic_and_taxonomic_biases_random_vs_minDiss_test_sets.eps", height = 12, width = 24)
grid.arrange(pC, pT, nrow = 1)
dev.off()

pdf("Taxonomic_biases_random_vs_minDiss_test_sets.pdf", 24, 12) 
grid.arrange(pTm, pTr, nrow = 2)
dev.off()


# --> similar balances.



### Decide what predictors to keep

# The HI are globally high (>0.81 for NA-train and >0.95 for NA-test), but we can see that some predictors have parts of their distribution in the non-assessed that is
# under-represented in the train and test datasets. This is especially true for:
# HFI (low values)
# HPD (low values)
# Tseas (low values)
# Pseas (low values)
# Tflex (high values)
# Pflex (high values)
# as visible on the Dataset_shift_minDiss_zoom_sc.pdf plots

# Make a list of these predictors to create a dataset without them later
removeShift <- c("HFI", "HPD", "Tseas", "Pseas", "Tflex", "Pflex")
write(removeShift, "predictors-to-remove_shift.txt")

# import all data and list of predictors to remove
all3_m_sc <- read.table("Data_split_minDiss_scaled.txt", stringsAsFactors = F, header=T) # data are not scaled, but the TEST split has been done after using scaled data
removeShift <- read.file("predictors-to-remove_shift.txt", header = F)[,1]

# make a dataset with everything but in a nice order, with no-predictor variables at the beginning
# keep other iucn category columns just in case want to predict other classes one day
all4 <- all3_m_sc[,1:2]
all4 <- cbind(all4, all3_m_sc[,17:21], all3_m_sc[,3:15])

# Generate an additional datasets without shifted predictors
  
all4_S <- all4[, !(names(all4) %in% removeShift)]

# Generate copies of these data where predictors are rescaled between 0 and 1

all4_sc <- all4
all4_S_sc <- all4_S

for (i in 8:length(names(all4_sc))) {
  all4_sc[,i] <- rescale(all4[,i])
}

for (i in 8:length(names(all4_S_sc))) {
  all4_S_sc[,i] <- rescale(all4_S[,i])
}


# export all datasets

write.table(all4, "Data_for_ML.txt", row.names = FALSE, sep="\t")
write.table(all4_S, "Data_for_ML_Shifted_removed.txt", row.names = FALSE, sep="\t")
write.table(all4_sc, "Data_for_ML_scaled.txt", row.names = FALSE, sep="\t")
write.table(all4_S_sc, "Data_for_ML_Shifted_removed_scaled.txt", row.names = FALSE, sep="\t")

