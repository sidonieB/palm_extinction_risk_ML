library(reshape2)
library(pROC)
library(ggplot2)
library(gridExtra)
library(caret)
library(gplots)
library(plyr)
library(DMwR)
library(dplyr)

### Predict the status of NE species using each trained model

# import trained models

trained_rf <- readRDS(file = "trained_models_rf.rds")
trained_nnD <- readRDS(file = "trained_models_nnD.rds")
trained_nnDS <- readRDS(file = "trained_models_nnDS.rds")
trained_rf_mt <- readRDS(file = "trained_models_rf_mt.rds")

trained_nn <- do.call(c, list(trained_nnD, trained_nnDS))

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

trained_rf <- trained_rf2

names(trained_rf) <- gsub("^M_", "RF___", names(trained_rf))
names(trained_nn) <- gsub("^M_", "NN___", names(trained_nn))
names(trained_rf_mt) <- gsub("^M_", "RFt___", names(trained_rf_mt))

# concatenate all models

trained_all <- vector(mode = "list")
trained_all <- append(trained_all, trained_rf)
trained_all <- append(trained_all, trained_rf_mt)
trained_all <- append(trained_all, trained_nn)

# export concatenated models for shap analysis later (separate script)
saveRDS(trained_all, file = "trained_models_all48.rds")
names(trained_all)


# import data for predictions and to combine

df <- read.table("Data_for_ML_scaled.txt", stringsAsFactors = F, sep="\t", header=T)
unused <- read.table("Unused_assessed.txt", stringsAsFactors = F, sep="\t", header=T)

# format data

pred_set <- subset(df, df$Subset == "NOT_ASSESSED")
train_set <- subset(df, df$Subset == "TRAIN")
test_set <- subset(df, df$Subset == "TEST")
unused2 <- unused[,c(1:2, 17:21, 3:15)]

# imput the NAs of the prediction set (in the same way we did for the test set before)
pred_set_noNA <- pred_set
pred_set_knn <- knnImputation(pred_set[,8:length(names(pred_set))], k=5, meth="median")
pred_set_noNA[,8:length(names(pred_set_noNA))] <- pred_set_knn

# predict the pred set without NAs with each model and append predictions to the original pred_set (that still has the NAs - not important)

for (m in 1:length(trained_all)) {
  trained_m <- trained_all[[m]]
  mname <- names(trained_all)[m]
  pred_m_cl <- predict(trained_m, pred_set_noNA)
  pred_set$pred_cl <- pred_m_cl
  names(pred_set)[length(names(pred_set))] <- paste(mname, "pred_cl", sep="_")
  pred_m_pr <- predict(trained_m, pred_set_noNA, type="prob")
  pred_set <- cbind(pred_set, pred_m_pr)
  names(pred_set)[(length(names(pred_set))-1):length(names(pred_set))] <- c(paste(mname, names(pred_m_pr)[1], sep="_"), paste(mname, names(pred_m_pr)[2], sep="_"))
}

# join to the other sets
df2 <- rbind(train_set, unused2)
df2_2 <- rbind(df2, test_set)
df3 <- bind_rows(df2_2, pred_set)

names(df3) <- gsub("_Kappa_rf", "", names(df3))
names(df3) <- gsub("_Kappa_nn", "", names(df3))
names(df3) <- gsub("_DPM_rf_mt", "", names(df3))

# export to be able to start from there if needed
write.table(df3, "Predictions_all_data_all_models.txt", sep="\t", row.names = F)

# create new dataset with only preds from chosen models and relevant iucn columns, i.e. without the predictors
#
# highest sensitivity
names(df3)
keep <- c(names(df3)[1], names(df3)[4], names(df3)[6:7], names(df3)[156:158]) # change 21:23 for desired model (chose highest precision)
df4 <- df3[,(names(df3) %in% keep)]
names(df4)

# rename last columns for easier handling
names(df4)[(length(names(df4))-2):length(names(df4))] <- c("Predicted_class", "Predicted_PnonLC", "Predicted_PLC")
names(df4)

# export
write.table(df4, "Predictions_selected_model_HighSe.txt", sep="\t", row.names = F)


# highest precision
names(df3)
keep <- c(names(df3)[1], names(df3)[4], names(df3)[6:7], names(df3)[114:116]) # change 21:23 for desired model (chose highest precision)
df4 <- df3[,(names(df3) %in% keep)]
names(df4)

# rename last columns for easier handling
names(df4)[(length(names(df4))-2):length(names(df4))] <- c("Predicted_class", "Predicted_PnonLC", "Predicted_PLC")
names(df4)

# export
write.table(df4, "Predictions_selected_model_HighP.txt", sep="\t", row.names = F)

