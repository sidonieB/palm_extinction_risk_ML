############################################
# Machine learning
############################################
# Before being predicted, the test set will have its NAs imputted using 5-NN median 
# otherwise may be done too roughly during the prediction: I could not make sure of how it would be done.
# (we will do the same with the pred set when we do the final predictions, but not here - we don't use the pred set here)
# The training set will be 5-NN imputted during the ML process using the preProcess function, to not bias the cv process (if input before, it creates data leakage between the in-bag and oob samples).

# levels of the conservation status variable are "releveled" so that the model sets nonLC as the positive class

# nnet is not deep learning: only one layer. Likely to be enough according to various sources
# number of neurones per layer (size parameter) could be set around the mean between number of input neurons (= number of predictors) and number of output neurons (= 1 if class or 2 if softmax)
# number of input is dep on preprocess. Usually number of variables after preProcess in rf was around 8, max 13. So could try size = c(3, 5, 7, 9, 11) 
# weight decay: too small=underfit, too large=overfit; try decay = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 8, 10, 12) (see how long it takes)
# number of iterations: try 100 (default) and see how it goes: 
# if model reached the max number of iterations, the "latter"convergence" output parameter will be 1 and the model may not have reached convergence, so we should increase the number of iterations

# new trial, with more size and decay values and more iterations only for Kappa, only for D data


library(lattice)
library(ggplot2)
library(caret)
library(gdata)
library(grid)
library(plyr)
library(DMwR)
library(ranger)
library(dplyr)
library(e1071)
library(randomForest)
library(nnet)
library(RANN)

## Import data and ensure numeric values are indeed set as numeric

D <- read.table("../../Data_for_ML_scaled.txt", header=T, sep="\t")
#DS <- read.table("../../Data_for_ML_Shifted_removed_scaled.txt", header=T, sep="\t")

D$Nocc <- as.numeric(D$Nocc)
D$NsubPop <- as.numeric(D$NsubPop)
D$Nloc <- as.numeric(D$Nloc)
D$AOO <- as.numeric(D$AOO)
D$ECO <- as.numeric(D$ECO)
D$TDWG3 <- as.numeric(D$TDWG3)

# DS$Nocc <- as.numeric(DS$Nocc)
# DS$NsubPop <- as.numeric(DS$NsubPop)
# DS$Nloc <- as.numeric(DS$Nloc)
# DS$AOO <- as.numeric(DS$AOO)
# DS$ECO <- as.numeric(DS$ECO)
# DS$TDWG3 <- as.numeric(DS$TDWG3)


## make lists/vectors of options (change values inside lists if want to try different things)
df_v <- list(D) #, DS)
names(df_v) <- c("D") #, "DS")
cor_cutoff_v <- c(0.75, 1)
sampling_option_v <- c("none", "smote", "up", "down")
metric_option_v <- c("Kappa")

## make a list in which to store all trained models and a list for the test results
model_num <- length(names(df_v)) * length(cor_cutoff_v) * length(sampling_option_v) * length(metric_option_v)
trained_models_nn <- vector(mode = "list", length = model_num)
test_results_nn <- vector(mode = "list", length = model_num)

## initialise the loop
v <- 0

## RUN FUNCTIONS THAT ARE AT THE END OF THE SCRIPT BEFORE PROCEEDING!!!!

## loop through these

for (w in 1:length(df_v)) {
  for (x in 1:length(cor_cutoff_v)) {
    for (y in 1:length(sampling_option_v)) {
      for (z in 1:length(metric_option_v)) {
        
        # train the model for the selected options using the function defined below
        # the function returns the trained model (trained_m), the model name (m), and the results when testing it on the test set (tested_m)
        iucn_model_nn(df=df_v[[w]], df_name=names(df_v)[w],cor_cutoff=cor_cutoff_v[x], sampling_option=sampling_option_v[y], metric_option=metric_option_v[z])
        
        v <- v + 1
        
        # append the trained model in the list and give it the model name
        trained_models_nn[[v]] <- trained_m
        names(trained_models_nn)[v] <- m_nn
        
        # same for the test results
        test_results_nn[[v]] <- test_set_wPreds
        names(test_results_nn)[v] <- m_nn
        
        # rename the trained model and the test results in case needed later
        assign(m_nn,trained_m)
        assign(paste(m_nn, "results", sep="_"),test_set_wPreds)
        
      }
    }
  }
}



### Save model objects

saveRDS(trained_models_nn, file = "trained_models_nn.rds")
saveRDS(test_results_nn, file = "test_results_nn.rds")

# check params selected, if at upper or lower bound, expand/shift range of params values
# check also convergence. If 0: ok, converged before the 100 iterations
params_df <- as.data.frame(names(trained_models_nn))
names(params_df)[1] <- "name"
params_df$size <- rep(0, length(params_df$name))
params_df$decay <- rep(0, length(params_df$name))
params_df$convergence <- rep(0, length(params_df$name))
for (i in 1:length(trained_models_nn)) {
  params_df$size[i] <- trained_models_nn[[i]]$finalModel$n[2]
  params_df$decay[i] <- trained_models_nn[[i]]$finalModel$decay
  params_df$convergence[i] <- trained_models_nn[[i]]$finalModel$convergence
}

write.table(params_df, "parameters_and_convergenceD.txt", sep = "\t")

# --> size at bounds and convergence not reached for some models: run them again longer and with more parameters
# other models are fine with the size and decay elected by the Kappa method: 
# we will run them again as ROC with the parameters fixed as found by Kappa so that we get ROC, Sens and Spec for them





### define function to train the model depending on options and to optimize size and decay

iucn_model_nn <- function(df, df_name, cor_cutoff, sampling_option, metric_option) {

## set option and model names

# set option names

cor_name <- as.character(cor_cutoff)
samp_name <- sampling_option
metric_name <- metric_option

# set model name accordingly, ensure it is set outside of the function too

m_nn <<- paste("M", df_name, cor_name, samp_name, metric_name, "nn", sep = "_")
print(paste("Running", m_nn))

## Reset sampling option to NULL if none because that's what is needed downstream
if (sampling_option == "none") {
  sampling_option <- NULL
}


## Split data and imput NAs in the test data

# split between training, test and to predict sets
pred_set <- subset(df, df$Subset == "NOT_ASSESSED")
train_set <- subset(df, df$Subset == "TRAIN")
test_set <- subset(df, df$Subset == "TEST")

# export subsets in case needed later (they will be the same for models using same df, and therefore overwrite each-other)
write.table(pred_set, paste(df_name, "toPredict_set.txt", sep="_"), row.names = FALSE, sep="\t")
write.table(train_set,paste(df_name, "training_set.txt", sep="_"), row.names = FALSE, sep="\t")
write.table(test_set, paste(df_name, "test_set.txt", sep="_"), row.names = FALSE, sep="\t")

# imput the NAs of the test set
test_set_knn <- knnImputation(test_set[,8:length(names(test_set))], k=5, meth="median")
test_set[,8:length(names(test_set))] <- test_set_knn

#seperate predictors and response variables for TRAINING and TESTING
# Careful that 8 still corresponds to the first predictor 
# (should be fine if did not modify previous scripts to shuffle/remove columns)
train_set_x <- train_set[,8:length(names(train_set))] # predictors
train_set_y <- relevel(as.factor(as.character(train_set$IUCN_TS_sum)), "nonLC") # response: conservation status, tweaked to have only the factor levels present in the subset otherwise error later
test_set_x <- test_set[,8:length(names(test_set))] # predictors
test_set_y <- relevel(as.factor(as.character(test_set$IUCN_TS_sum)), "nonLC") # response: conservation status


## Set the training parameters
# read the options selected
# always try multiple size and decay
# always do 5-NN NA inputtation (preliminary trials show similar performance than bagging and much quicker than the latter)

# depending on the metric, the summaryFunction option will be different
if (metric_option == "ROC") {
  summaryFunction_option <- twoClassSummary
} else if (metric_option == "Kappa") {
  summaryFunction_option <- defaultSummary
} else {
  print("ERROR, wrong metric option")
}


# set seeds so that reproducible (very important to do it here AND once before running train)
# we are doing 10 repeats of 10-Fold CV. We will fit a model that evaluates as many as 5*10 combinations of size and decay parameter values
# We also add 20 to that value to have seeds set in the same way than when we also tune the threshold in other runs
# and we need a value for the last model so:
repeats <- 10
folds <- 10
mtry_num <- 101
length_seeds <- (repeats * folds) + 1
length_seeds_iter <- mtry_num
my_seeds <- vector(mode = "list", length = length_seeds)
set.seed(1)
for(i in 1:(length_seeds - 1)) my_seeds[[i]] <- sample.int(1000, length_seeds_iter)
set.seed(1)
my_seeds[[length_seeds]] <- sample.int(1000, 1)


# all options are provided for easier understanding, reproductibility and to facilitate later changes.
m_control <- trainControl(    
  preProcOptions = list(cutoff = cor_cutoff), # set correlation threshold for predictors removal
  method = "repeatedcv", number = 10, repeats = 10, # 10-fold CV, repeat 10 times
  search = "grid", # tuning on the grid of parameter values provided below, instead of on random parameter values
  verboseIter = TRUE, # print training progress
  returnData = TRUE, # saves the data
  returnResamp = "final", # return performance on the resampled for the optimal tuning parameters - change to all if want to plot performance depending on parameter values
  savePredictions = "final",# saves the predictions for the held out samples, but only for the optimal tuning parameters - can change to "all" or "none"
  classProbs = TRUE, # compute class probability in addition to the predicted class
  summaryFunction = summaryFunction_option, # summarizes performance accross resamples. Use this because we predict two classes. It will give the area under the ROC, sensitivity and specificity - see defaultSummary()
                                     # could probably leave it to defaultSummary, in which case I think it defaults to the "metric" defined in train below, but I am not sure.
  selectionFunction = "best", # choses the parameter values associated with the highest model performance to make the final model. Performance is measured using Accuracy or Kappa (better if unbalanced classes) or something else as defined by the "metric" option in train below.
  sampling=sampling_option, # additional sampling conducted after the resampling for the cv. Here we can deal with class unbalance.
  index = NULL, # in case do not want folds of the cv to be random
  indexOut = NULL, # in case do not want folds of the cv to be random
  indexFinal = NULL,  # in case do not want to fit the final model on all samples
  timingSamps = 0, # change if want to measure how long does the prediction take (see help)
  predictionBounds = rep(FALSE,2), # do not constrain prediction outcome (see help if want to constrain it)
  seeds = my_seeds, # set seeds internally, important!
  #adaptive = list(min = 5, alpha = 0.05, method ="gls", complete = TRUE), # only relevant if use an adaptive method
  trim = FALSE, # set to TRUE if wants to make space by removing some model attributes (may prevent some downstream analyses)
  allowParallel = TRUE) # to prallelize, if available


# Set the preprocessing method for below (m_pp_method), depending on the correlation cutoff
# The cor function does not accept a cutoff of 1 (nor a cutoff of 0.95 apparently but not sure if the case for all data so leave it at 1
# change 1 to something lower if get an error)
if (cor_cutoff < 1) {
m_pp_method <- c("corr", "knnImpute")
} else {
m_pp_method <- c("knnImpute")
}

# Set a grid to try multiple values of size and decay
myGrid <- expand.grid(size = c(1,2,3,5,7, 9, 11, 13), decay = c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 8, 10))
# myGrid <- expand.grid(size = c(1,2,3,5,7, 9, 11, 13), decay = c(0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 8, 10))
# myGrid <- expand.grid(size = 3, decay = 0.1)

capture.output(myGrid, file=paste(m_nn, "myGrid.txt", sep="_"))

## Train the model and assign the result outside of the function
set.seed(1)
trained_m <<- train(x = train_set_x, y = train_set_y, 
                  metric = metric_option, # assessment of model performance during tuning will be done with this metric
                  method = "nnet",
                  preProcess= m_pp_method, # preprocessing method defined above
                  trControl = m_control, # training parameters defined above
                  tuneGrid = myGrid,  # tuning parameters defined above
                  maxit = 1000, # number of iterations
                  importance = T) # assess variable imporatance


## Test the best model (the one trained with the best mtry) and assign results out of the function

tested_m_cl <- predict(trained_m, test_set)
test_set$pred_cl <- tested_m_cl
tested_m_pr <- predict(trained_m, test_set, type="prob")
test_set_wPreds <<- cbind(test_set, tested_m_pr)


## Export outputs

# model training and test details
capture.output(trained_m, file=paste(m_nn, "training_stdout.txt", sep="_")) 
write.table(trained_m$results, paste(m_nn, "training_results.txt", sep="_"), row.names = FALSE)
write.table(trained_m$pred, paste(m_nn, "training_cv-preds.txt", sep="_"), row.names = FALSE)
write.table(test_set_wPreds, paste(m_nn, "test_results.txt", sep="_"), row.names = FALSE)

# Export also the confusion matrix on the resamples and on the test
# change the positive class to be nLC and get precision in addition to sensitivity and specificity
# specificity will be fraction of total negative (LC) samples that are well predicted
# sensitivity will be fraction of total positive (nLC) samples that are well predicted = recall
# precision will be fraction of samples predicted as the positive class (nLC) that are indeed of that class (nLC)

cm2 <- confusionMatrix(data = trained_m$pred$pred, reference = trained_m$pred$ob, mode = "everything")
capture.output(cm2, file=paste(m_nn, "training_ConfMat.txt", sep="_"))

cm3 <- confusionMatrix(data=as.factor(test_set_wPreds$pred_cl), reference=relevel(as.factor(as.character(test_set_wPreds$IUCN_TS_sum)), "nonLC"), mode = "everything")
capture.output(cm3, file=paste(m_nn, "testing_ConfMat.txt", sep="_"))


}







