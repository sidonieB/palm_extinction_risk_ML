### Compare models


library(reshape2)
library(pROC)
library(ggplot2)
library(gridExtra)
library(caret)
library(gplots)
library(plyr)
library(dplyr)
library(scales)
library(stringr)
library(tidyr)
library(ggfortify)
library(DMwR)
library(randomForest)
library(ggpubr)


# Import results with final options saved
trained_rf <- readRDS(file = "trained_models_rf.rds")
test_rf <- readRDS(file = "test_results_rf.rds")
trained_nnD <- readRDS(file = "trained_models_nnD.rds")
test_nnD <- readRDS(file = "test_results_nnD.rds")
trained_nnDS <- readRDS(file = "trained_models_nnDS.rds")
test_nnDS <- readRDS(file = "test_results_nnDS.rds")
trained_rf_mt <- readRDS(file = "trained_models_rf_mt.rds")
test_rf_mt <- readRDS(file = "test_results_rf_mt.rds")

trained_nn <- do.call(c, list(trained_nnD, trained_nnDS))
test_nn <- do.call(c, list(test_nnD, test_nnDS))

# keep only the Kappa models for rf
trained_rf2 <- vector(mode = "list", length = length(trained_nn))
y = 0
for (x in 1:length(trained_rf)) {
  if (length(grep("Kappa", names(trained_rf)[x])) > 0) {
    y = y+ 1
    trained_rf2[[y]] <- trained_rf[[x]]
    names(trained_rf2)[y] <- names(trained_rf)[x]
  }
}


test_rf2 <- vector(mode = "list", length = length(test_nn))
y = 0
for (x in 1:length(test_rf)) {
  if (length(grep("Kappa", names(test_rf)[x])) > 0) {
    y = y+ 1
    test_rf2[[y]] <- test_rf[[x]]
    names(test_rf2)[y] <- names(test_rf)[x]
  }
}

trained_rf <- trained_rf2
test_rf <- test_rf2

names(trained_rf) <- gsub("^M_", "RF___", names(trained_rf))
names(trained_nn) <- gsub("^M_", "NN___", names(trained_nn))
names(trained_rf_mt) <- gsub("^M_", "RFt___", names(trained_rf_mt))

names(test_rf) <- gsub("^M_", "RF___", names(test_rf))
names(test_nn) <- gsub("^M_", "NN___", names(test_nn))
names(test_rf_mt) <- gsub("^M_", "RFt___", names(test_rf_mt))




# Make a table of all interesting indicators obtained from the test data

test_df <- as.data.frame(c(names(trained_rf), names(trained_nn), names(trained_rf_mt)))
names(test_df)[1] <- "model"
test_df$model <- gsub("_Kappa.*|_ROC.*|_DPM.*", "", test_df$model)
test_df$Method <- rep(NA, length(test_df$model))
test_df$Model <- rep(NA, length(test_df$model))
for (x in 1:length(test_df$model)) {
test_df$Method[x] <- strsplit(test_df$model[x], "___")[[1]][1]
test_df$Model[x] <- strsplit(test_df$model[x], "___")[[1]][2]
}
test_df$Sample <- rep(NA, length(test_df$model))
test_df$AUCROC <- rep(0, length(test_df$model))
test_df$Sensitivity <- rep(0, length(test_df$model))
test_df$Specificity <- rep(0, length(test_df$model))
test_df$DPM <- rep(0, length(test_df$model))
test_df$Accuracy <- rep(0, length(test_df$model))
test_df$Balanced_accuracy <- rep(0, length(test_df$model))
test_df$Kappa <- rep(0, length(test_df$model))
test_df$Precision <- rep(0, length(test_df$model))
test_df$NegativePrecision <- rep(0, length(test_df$model))
test_df$Recall <- rep(0, length(test_df$model))
test_df$F1score <- rep(0, length(test_df$model))


# Fill the table - Careful!!! Relevel the test observations to get nonLC as positive class

all_models <- do.call(c, list(trained_rf, trained_nn, trained_rf_mt))
all_tests <- do.call(c, list(test_rf, test_nn, test_rf_mt))

for (m in 1:length(all_models)) {
  mname <- names(all_models)[m]
  mname2 <- gsub("_Kappa.*|_ROC.*|_DPM.*", "", mname)
  sname <- paste(mname2, "Test", sep="_")
  test_df$Sample[which(test_df$model == mname2)] <- sname
  
  roc_m <- roc(all_tests[[m]]$IUCN_TS_sum, predict(all_models[[m]], all_tests[[m]][,8:(length(names(all_tests[[m]]))-3)], type = "prob")[,1], 
               levels = rev(levels(as.factor(as.character(all_tests[[m]]$IUCN_TS_sum)))))
  
  test_df$AUCROC[which(test_df$model == mname2)] <- roc_m$auc[1]
  
  cm_m <- confusionMatrix(data=as.factor(all_tests[[m]]$pred_cl), reference=relevel(as.factor(as.character(all_tests[[m]]$IUCN_TS_sum)), "nonLC"), mode = "everything")
  
  test_df$Accuracy[which(test_df$model == mname2)] <- cm_m$overall[[1]]
  test_df$Kappa[which(test_df$model == mname2)] <- cm_m$overall[[2]]
  test_df$Balanced_accuracy[which(test_df$model == mname2)] <- cm_m$byClass[[11]]
  test_df$Sensitivity[which(test_df$model == mname2)] <- cm_m$byClass[[1]]
  test_df$Specificity[which(test_df$model == mname2)] <- cm_m$byClass[[2]]
  Sens <- cm_m$byClass[[1]]
  Spec <- cm_m$byClass[[2]]
  test_df$Precision[which(test_df$model == mname2)] <- cm_m$byClass[[5]]
  test_df$NegativePrecision[which(test_df$model == mname2)] <- cm_m$byClass[[4]]
  test_df$Recall[which(test_df$model == mname2)] <- cm_m$byClass[[6]]
  test_df$F1score[which(test_df$model == mname2)] <- cm_m$byClass[[7]]
  coords <- matrix(c(1, 1, Spec, Sens), ncol = 2, byrow = TRUE)
  colnames(coords) <- c("Spec", "Sens")
  rownames(coords) <- c("Best", "Current")
  test_df$DPM[which(test_df$model == mname2)] <- dist(coords)[1]
}

# export table
write.table(test_df, "Performance_AllModels.txt", sep="\t", row.names = F)



# Plot results for each test species (heatmap)

# join predictions of all models as P(nonLC) or class (for mt), add one column at the beginning for the "Observed"
# transform the obs, or *_mt classes as probabilities (0 and 1 with 1 = nonLC)

# build the dataset by starting to make a column for the real assessment value of each test species
test_preds <- as.data.frame(test_rf[[1]][,1])
names(test_preds)[1] <- "Species"
test_preds$Observed <- as.character(test_rf[[1]]$IUCN_TS_sum)
for (z in 1:length(test_preds$Observed)) {
  if (test_preds$Observed[z] == "nonLC") {
    test_preds$Observed[z] <- 1
  } else {
    test_preds$Observed[z] <- 0
  }
}

# fill it with the test results
for (x in 1:length(all_tests)){
  test_preds$newC <- all_tests[[x]]$nonLC
  if (length(grep("RFt___", names(all_tests)[x])) > 0) {
    print("yes")
  test_preds$newC <- as.character(all_tests[[x]]$pred_cl)
  for (z in 1:length(test_preds$newC)) {
    if (test_preds$newC[z] == "nonLC") {
      test_preds$newC[z] <- 1
    } else {
      test_preds$newC[z] <- 0
    }
  }
  }
  names(test_preds)[length(test_preds)] <- gsub("_Kappa.*|_ROC.*|_DPM.*", "", names(all_tests)[x])
}


tail(test_preds)


# do heatmap

# change df to matrix after ordering it conveniently and using Species as row names instead of a column
test_preds2 <- test_preds[order(test_preds$Observed, test_preds$Species, decreasing = F),]
row.names(test_preds2) <- test_preds2$Species
test_preds2 <- test_preds2[,2:length(names(test_preds2))]
test_preds3 <- test_preds2
test_preds3[1:length(names(test_preds3))] <- sapply(test_preds2[1:length(names(test_preds2))],as.numeric)

test_preds2_m  <- data.matrix(test_preds3)

head(test_preds2_m)
tail(test_preds2_m)

# heatmap
# custom legend with colors for bins 0.2
breaks = c(0, 0.25, 0.5, 0.75, 1)
hm.colors = c('#0571b0','#92c5de', '#f4a582','#ca0020')


pdf("Performance_test_heatmap.pdf", 42, 60) 
heatmap.2(test_preds2_m, dendrogram="none", Rowv=NA, 
          Colv=NA, breaks=breaks, col = hm.colors, 
          scale="none", trace="none", margins=c(25,15), keysize = 1, key.xlab = "P(non-LC)", key.title = "",
          density.info="none", cexRow = 4,
          cexCol = 4, offsetRow=-235, offsetCol = -310)
dev.off()



setEPS()
postscript("Performance_test_heatmap.eps", height = 60, width = 42)
heatmap.2(test_preds2_m, dendrogram="none", Rowv=NA, 
          Colv=NA, breaks=breaks, col = hm.colors, 
          scale="none", trace="none", margins=c(25,15), keysize = 1, key.xlab = "P(non-LC)", key.title = "",
          density.info="none", cexRow = 4,
          cexCol = 4, offsetRow=-235, offsetCol = -310)

dev.off()


# Score species depending on the prediction success 
# (add 1 for each model that classified it well - 0.5 is not a good classification)
test_preds2$ClassificationScore <- rep(NA, length(test_preds2$Observed))
for (x in 1:length(test_preds2$Observed)) {
  y = 0
    if (test_preds2$Observed[x] == 0) {
      for (i in 2:(length(names(test_preds2))-1)) {
        if (test_preds2[x,i] < 0.5) {
          y = y+1
        }
      }
    } else {
      for (i in 2:(length(names(test_preds2))-1)) {
        if (test_preds2[x,i] > 0.5) {
          y = y+1
        }        
      }
    }
  test_preds2$ClassificationScore[x] <- y
}

# Add a column to have it as a proportion of models which got it right, i.e. score*100/48
test_preds2$ClassificationScoreProp <- round(test_preds2$ClassificationScore*100/48)





# Explore links between classification score, number of occurrences, assessment method and assessment category
# to see if consistently wrongly predicted species have things in commmon

# get original predictor data (not scaled) and add columns to make clear distinction between Threat Search vs IUCN, and recode categories for consistency
# do it only for the test set

preds <- read.table("Data_for_ML.txt", stringsAsFactors = F, header=T)
preds_sub <- subset(preds, preds$Subset == "TEST")
preds_sub$Assessment_scheme <- rep(NA, length(preds_sub$accepted_name)) 
preds_sub$Subcategory_homogenized <- rep(NA, length(preds_sub$accepted_name))
for (x in 1:length(preds_sub$accepted_name)) {
  if (preds_sub$IUCN_TS_detail[x] == "Near Threatened") {
    preds_sub$Assessment_scheme[x] <- "TS"
    preds_sub$Subcategory_homogenized[x] <- "NT"
  } else if (preds_sub$IUCN_TS_detail[x] == "Not Threatened") {
    preds_sub$Assessment_scheme[x] <- "TS"
    preds_sub$Subcategory_homogenized[x] <- "LC"    
  } else if (preds_sub$IUCN_TS_detail[x] == "Possibly Threatened") {
    preds_sub$Assessment_scheme[x] <- "TS"
    preds_sub$Subcategory_homogenized[x] <- "PT"    
  } else if (preds_sub$IUCN_TS_detail[x] == "Threatened,Data Deficient") {
    preds_sub$Assessment_scheme[x] <- "TS"
    preds_sub$Subcategory_homogenized[x] <- "T"    
  } else if (preds_sub$IUCN_TS_detail[x] == "Threatened") {
    preds_sub$Assessment_scheme[x] <- "TS"
    preds_sub$Subcategory_homogenized[x] <- "T"    
  } else {
    preds_sub$Assessment_scheme[x] <- "RL"
    preds_sub$Subcategory_homogenized[x] <- preds_sub$IUCN_TS_detail[x]    
  }
}

# join this with the classification scores
score <- data.frame(score = test_preds2$ClassificationScore)
score$accepted_name <- rownames(test_preds2)
score$scoreProp <- test_preds2$ClassificationScoreProp
preds2 <- full_join(preds_sub, score, by="accepted_name")
preds2_2 <- preds2
names(preds2_2)[1] <- "Species"
preds6 <- full_join(preds2_2, test_preds, by="Species")


# Number of occurrences

p_occ <- ggplot(data=preds2, aes(x = Nocc, y = scoreProp, color=IUCN_TS_sum)) +
geom_point() +
scale_color_manual(values=c('#0571b0','#ca0020'))+ 
#geom_vline(xintercept = 60, linetype="dashed", color = "grey40") +
theme_bw()+
theme(legend.position="none")


# same but add correlation test and distinguishes tests for threatened or non-threatened species
# also try after removing outliers

p_occ_2 <- ggscatter(data=preds2, x = "Nocc", y = "scoreProp", add = "reg.line") +
  #geom_point() +
  #scale_color_manual(values=c('#0571b0','#ca0020'))+ 
  #geom_smooth(method=lm, se=FALSE, fullrange=F)+
  stat_cor(method = "spearman", label.x = 200, label.y = 80)+
  #geom_vline(xintercept = 100, linetype="dashed", color = "grey40") +
  theme_bw()+
  theme(legend.position="none")

p_occ_2_sep <- ggscatter(data=preds2, x = "Nocc", y = "scoreProp", color="IUCN_TS_sum", add = "reg.line") +
  #geom_point() +
  scale_color_manual(values=c('#0571b0','#ca0020'))+ 
  #geom_smooth(method=lm, se=FALSE, fullrange=F)+
  stat_cor(method = "spearman", aes(color = IUCN_TS_sum), label.x = 200, label.y = c(80, 75))+
  #geom_vline(xintercept = 100, linetype="dashed", color = "grey40") +
  theme_bw()+
  theme(legend.position="none")


# remove 5 species with >= 100 occs
preds2_inf100 <- preds2[preds2$Nocc < 100,]

  p_occ_3 <- ggscatter(data=preds2_inf100, x = "Nocc", y = "scoreProp", add = "reg.line") +
    #geom_point() +
    #scale_color_manual(values=c('#0571b0','#ca0020'))+ 
    #geom_smooth(method=lm, se=FALSE, fullrange=F)+
    stat_cor(method = "spearman", label.x = 60, label.y = 80)+
  #geom_vline(xintercept = 100, linetype="dashed", color = "grey40") +
  theme_bw()+
    theme(legend.position="none")
  
  p_occ_3_sep <- ggscatter(data=preds2_inf100, x = "Nocc", y = "scoreProp", color="IUCN_TS_sum", add = "reg.line") +
    #geom_point() +
    scale_color_manual(values=c('#0571b0','#ca0020'))+ 
    #geom_smooth(method=lm, se=FALSE, fullrange=F)+
    stat_cor(method = "spearman", aes(color = IUCN_TS_sum), label.x = 50, label.y = c(80, 75))+
    #geom_vline(xintercept = 100, linetype="dashed", color = "grey40") +
    theme_bw()+
    theme(legend.position="none")

  
  
# same for the best-case and worst case models only
  preds7 <- preds6
  preds7$RFt___DS_1_down <- as.numeric(preds7$RFt___DS_1_down)
  
  p_occ_sep_BS <- ggscatter(data=preds7, x = "Nocc", y = "RFt___DS_1_down", color="IUCN_TS_sum", add = "reg.line") +
    #geom_point() +
    scale_color_manual(values=c('#0571b0','#ca0020'))+ 
    #geom_smooth(method=lm, se=FALSE, fullrange=F)+
    stat_cor(method = "spearman", aes(color = IUCN_TS_sum), label.x = 200, label.y = c(1.5, 1.4))+
    #geom_vline(xintercept = 100, linetype="dashed", color = "grey40") +
    theme_bw()+
    theme(legend.position="none")  

  p_occ_sep_WS <- ggscatter(data=preds7, x = "Nocc", y = "NN___DS_1_smote", color="IUCN_TS_sum", add = "reg.line") +
    #geom_point() +
    scale_color_manual(values=c('#0571b0','#ca0020'))+ 
    #geom_smooth(method=lm, se=FALSE, fullrange=F)+
    stat_cor(method = "spearman", aes(color = IUCN_TS_sum), label.x = 200, label.y = c(0.6, 0.65))+
    geom_hline(yintercept = 0.5, linetype="dashed", color = "grey40") +
    theme_bw()+
    theme(legend.position="none")
  
  # remove 5 species with >= 100 occs
  preds7_inf100 <- preds7[preds7$Nocc < 100,]
  
  p_occ_sep_BS2 <- ggscatter(data=preds7_inf100, x = "Nocc", y = "RFt___DS_1_down", color="IUCN_TS_sum", add = "reg.line") +
    #geom_point() +
    scale_color_manual(values=c('#0571b0','#ca0020'))+ 
    #geom_smooth(method=lm, se=FALSE, fullrange=F)+
    stat_cor(method = "spearman", aes(color = IUCN_TS_sum), label.x = 25, label.y = c(1.5, 1.4))+
    #geom_vline(xintercept = 100, linetype="dashed", color = "grey40") +
    theme_bw()+
    theme(legend.position="none")  
  
  p_occ_sep_WS2 <- ggscatter(data=preds7_inf100, x = "Nocc", y = "NN___DS_1_smote", color="IUCN_TS_sum", add = "reg.line") +
    #geom_point() +
    scale_color_manual(values=c('#0571b0','#ca0020'))+ 
    #geom_smooth(method=lm, se=FALSE, fullrange=F)+
    stat_cor(method = "spearman", aes(color = IUCN_TS_sum), label.x = 25, label.y = c(0.65, 0.7))+
    geom_hline(yintercept = 0.5, linetype="dashed", color = "grey40") +
    theme_bw()+
    theme(legend.position="none")
  

 
# Assessment scheme and subcategory
  
preds2$Assessment_scheme <- as.factor(preds2$Assessment_scheme)
preds2$Subcategory_homogenized <- as.factor(preds2$Subcategory_homogenized)
preds7$Assessment_scheme <- as.factor(preds7$Assessment_scheme)
preds7$Subcategory_homogenized <- as.factor(preds7$Subcategory_homogenized)

#write.table(preds2, "Test_species_scores.txt", row.names = F, sep="\t")
write.table(preds6, "Test_species_scores_allModels.txt", row.names = F, sep="\t")

p_TS <- ggplot(data=preds2, aes(x=Assessment_scheme, y=scoreProp, fill=IUCN_TS_sum)) + 
  #geom_boxplot(outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1, binwidth = 1)+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#0571b0','#ca0020'))+
  theme_bw()+
  theme(legend.position="none")

p_TS_BS <- ggplot(data=preds7, aes(x=Assessment_scheme, y=RFt___DS_1_down, fill=IUCN_TS_sum)) + 
  #geom_boxplot(outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=0.03, binwidth = 0.5)+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#0571b0','#ca0020'))+
  theme_bw()+
  theme(legend.position="none")

p_TS_WS <- ggplot(data=preds7, aes(x=Assessment_scheme, y=NN___DS_1_smote, fill=IUCN_TS_sum)) + 
  #geom_boxplot(outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=0.1, binwidth = 0.1)+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#0571b0','#ca0020'))+
  theme_bw()+
  geom_hline(yintercept = 0.5, linetype="dashed", color = "grey40") +
  theme(legend.position="none")


p_Cat <- ggplot(data=preds2, aes(x=Subcategory_homogenized, y=scoreProp, fill=IUCN_TS_sum)) + 
  #geom_boxplot(outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=1, binwidth = 1)+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#0571b0','#ca0020'))+
  theme_bw()+
  theme(legend.position="none")

p_Cat_BS <- ggplot(data=preds7, aes(x=Subcategory_homogenized, y=RFt___DS_1_down, fill=IUCN_TS_sum)) + 
  #geom_boxplot(outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=0.03, binwidth = 0.5)+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#0571b0','#ca0020'))+
  theme_bw()+
  theme(legend.position="none")

p_Cat_WS <- ggplot(data=preds7, aes(x=Subcategory_homogenized, y=NN___DS_1_smote, fill=IUCN_TS_sum)) + 
  #geom_boxplot(outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), dotsize=0.1, binwidth = 0.1)+
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#0571b0','#ca0020'))+
  theme_bw()+
  geom_hline(yintercept = 0.5, linetype="dashed", color = "grey40") +
  theme(legend.position="none")


pdf("Performance_Occ_TS_Cat.pdf", 15, 15) 
grid.arrange(p_TS, p_TS_BS, p_TS_WS, 
             p_Cat, p_Cat_BS, p_Cat_WS, 
             p_occ_3_sep, p_occ_sep_BS2, p_occ_sep_WS2, nrow=3)
dev.off()

setEPS()
postscript("Performance_Occ_TS_Cat.eps", height = 15, width = 15)
grid.arrange(p_TS, p_TS_BS, p_TS_WS, 
             p_Cat, p_Cat_BS, p_Cat_WS, 
             p_occ_3_sep, p_occ_sep_BS2, p_occ_sep_WS2, nrow=3)
dev.off()




