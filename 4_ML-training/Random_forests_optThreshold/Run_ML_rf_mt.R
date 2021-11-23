############################################
# Machine learning
############################################
# Before being predicted, the test set will have its NAs imputted using 5-NN median 
# otherwise may be done too roughly during the prediction: I could not make sure of how it would be done.
# (we will do the same with the pred set when we do the final predictions, but not here - we don't use the pred set here)
# The training set will be 5-NN imputted during the ML process using the preProcess function, to not bias the cv process (if input before, it creates data leakage between the in-bag and oob samples).

# levels of the conservation status variable are "releveled" so that the model sets nonLC as the positive class

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

## Import data and ensure numeric values are indeed set as numeric

D <- read.table("../Data_for_ML_scaled.txt", header=T, sep="\t")
DS <- read.table("../Data_for_ML_Shifted_removed_scaled.txt", header=T, sep="\t")

D$Nocc <- as.numeric(D$Nocc)
D$NsubPop <- as.numeric(D$NsubPop)
D$Nloc <- as.numeric(D$Nloc)
D$AOO <- as.numeric(D$AOO)
D$ECO <- as.numeric(D$ECO)
D$TDWG3 <- as.numeric(D$TDWG3)

DS$Nocc <- as.numeric(DS$Nocc)
DS$NsubPop <- as.numeric(DS$NsubPop)
DS$Nloc <- as.numeric(DS$Nloc)
DS$AOO <- as.numeric(DS$AOO)
DS$ECO <- as.numeric(DS$ECO)
DS$TDWG3 <- as.numeric(DS$TDWG3)



## make lists/vectors of options (change values inside lists if want to try different things)
df_v <- list(D, DS)
names(df_v) <- c("D", "DS")
cor_cutoff_v <- c(0.75, 1)
sampling_option_v <- c("none", "smote", "up", "down")
#metric_option_v <- c("ROC", "Kappa") # the metric is now the distance to perfect model (DPM), set directly in the new function

## make a list in which to store all trained models and a list for the test results
model_num <- length(names(df_v)) * length(cor_cutoff_v) * length(sampling_option_v) #* length(metric_option_v)
trained_models_rf_mt <- vector(mode = "list", length = model_num)
test_results_rf_mt <- vector(mode = "list", length = model_num)

## initialise the loop
v <- 0

## RUN FUNCTIONS THAT ARE AT THE END OF THE SCRIPT BEFORE PROCEEDING!!!!

## loop through these

for (w in 1:length(df_v)) {
  for (x in 1:length(cor_cutoff_v)) {
    for (y in 1:length(sampling_option_v)) {
     # for (z in 1:length(metric_option_v)) {
        
        # train the model for the selected options using the function defined below
        # the function returns the trained model (trained_m), the model name (m), and the results when testing it on the test set (tested_m)
        iucn_model_rf_mt(df=df_v[[w]], df_name=names(df_v)[w],cor_cutoff=cor_cutoff_v[x], sampling_option=sampling_option_v[y]) #, metric_option=metric_option_v[z])
        
        v <- v + 1
        
        # append the trained model in the list and give it the model name
        trained_models_rf_mt[[v]] <- trained_m
        names(trained_models_rf_mt)[v] <- m_rf_mt
        
        # same for the test results
        test_results_rf_mt[[v]] <- test_set_wPreds
        names(test_results_rf_mt)[v] <- m_rf_mt
        
        # rename the trained model and the test results in case needed later
        assign(m_rf_mt,trained_m)
        assign(paste(m_rf_mt, "results", sep="_"),test_set_wPreds)
        
      }
    }
  }
#}


### Save model objects

saveRDS(trained_models_rf_mt, file = "trained_models_rf_mt.rds")
saveRDS(test_results_rf_mt, file = "test_results_rf_mt.rds")










### define function to train the model depending on options

iucn_model_rf_mt <- function(df, df_name, cor_cutoff, sampling_option) {
  
  ## set option and model names
  
  # set option names
  
  cor_name <- as.character(cor_cutoff)
  samp_name <- sampling_option
  
  # set model name accordingly, ensure it is set outside of the function too
  
  m_rf_mt <<- paste("M", df_name, cor_name, samp_name, "DPM_rf_mt", sep = "_")
  print(paste("Running", m_rf_mt))
  
  ## Reset sampling option to NULL if none because that's what is needed downstream
  if (sampling_option == "none") {
    sampling_option <- NULL
  }
  

  ## Split data
  
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
  
  
  # set seeds so that reproducible (very important to do it here AND once before running train)
  # we are doing 10 repeats of 10-Fold CV. We will fit a model that evaluates as much as 13 values of mtry 
  # (well in fact 12 because number of independant variables) and 20 values of threshold
  # and we need a value for the last model so:
  repeats <- 10
  folds <- 10
  mtry_num <- 33
  length_seeds <- (repeats * folds) + 1
  length_seeds_iter <- mtry_num
  my_seeds <- vector(mode = "list", length = length_seeds)
  set.seed(1)
  for(i in 1:(length_seeds - 1)) my_seeds[[i]] <- sample.int(1000, length_seeds_iter)
  set.seed(1)
  my_seeds[[length_seeds]] <- sample.int(1000, 1)
  
  ## Set the training parameters
  # read the options selected
  # always preProcess with range with rangeBounds set to (0,1) (so that knn inputtation not biased)
  # mtry set to optimal mtry
  # always do k-NN NA inputtation (preliminary trials show similar performance than bagging and much quicker than the latter)
  
  m_control <- trainControl(    
    preProcOptions = list(cutoff = cor_cutoff),
    method = "repeatedcv", number = 10, repeats = 10,
    classProbs = TRUE,
    savePredictions = "final",
    summaryFunction = fourStats, # custom
    seeds = my_seeds,
    sampling=sampling_option)
  
  

  # Assess the number of variables that will be left after the preprocessing 
  # They will be stored in m_pp$method$remove
  # Also use this opportunity to set the preprocessing method for below (m_pp_method), depending on the correlation cutoff
  # The cor function does not accept a cutoff of 1 (nor a cutoff of 0.95 apparently but not sure if the case for all data so leave it at 1
  # change 1 to something lower if get an error)
  if (cor_cutoff < 1) {
    m_pp_method <- c("corr", "knnImpute")
  } else {
    m_pp_method <- c("knnImpute")
  }
  
  m_pp <- preProcess(train_set_x,
                     method = m_pp_method,
                     cutoff = cor_cutoff)
  
  # Set a grid to try multiple values of mtry 
  # Controls how many predictors are exposed to the splitting search routine at each split
  # More increases correlation among trees, ie bad, but also increases individual tree performance, ie good --> tradeoff
  
  # Set the number of predictors to try to be equal to the number of predictors that will be left in the df after correlated predictors have been removed minus one
  if (cor_cutoff < 1) {
    opt_mtry <- 1:(length(names(train_set_x)) - length(m_pp$method$remove) - 1)
  } else {
    opt_mtry <- 1:(length(names(train_set_x)) - 1) 
  }
  
  
 
  ## train the model and find optimal mtry+threshold (the one that minimises DPM) 
  # could not find a way to pass mtry here or above in trainControl so the function using mtry has to be sourced for each new opt_mtry
  # mtry_fixed is a misleading name kept because used elsewhere: in fact the mtry is optimized during the training, together with the threshold
  
  thresh_code$grid <- function(x, y, len = NULL, search = "grid", mtry_fixed=opt_mtry) {
    if(search == "grid") {
      grid <- expand.grid(mtry = mtry_fixed, threshold = seq(.01, .99, length = len))
    } else {
      grid <- expand.grid(mtry = mtry_fixed, threshold = runif(runif, min = 0, max = 1))
    }
    grid
  }
  
  
  
  set.seed(1)
  trained_m <<- train(x = train_set_x, y = train_set_y, 
                     metric = "DPM", # custom: the distance to the perfect model (see function defined below)
                     method = thresh_code, # custom (defined below)
                     maximize = FALSE, # custom (will minimize the new metric DPM)
                     tuneLength = 20, # custom (trying 20 values of threshold between 0 and 1)
                     preProcess= m_pp_method,
                     trControl = m_control,
                     num.trees = 1000,
                     importance = T)
  
  

    
  ## Test the best model (the one trained with the best mtry+threshold) and assign results out of the function
  
  tested_m_cl <- predict(trained_m, test_set)
  test_set$pred_cl <- tested_m_cl
  tested_m_pr <- predict(trained_m, test_set, type="prob")
  test_set_wPreds <<- cbind(test_set, tested_m_pr)
  
  
  ## Export outputs
  
  # model training and test details
  capture.output(trained_m, file=paste(m_rf_mt, "training_stdout.txt", sep="_")) 
  write.table(trained_m$results, paste(m_rf_mt, "training_results.txt", sep="_"), row.names = FALSE)
  write.table(trained_m$pred, paste(m_rf_mt, "training_cv-preds.txt", sep="_"), row.names = FALSE)
  write.table(test_set_wPreds, paste(m_rf_mt, "test_results.txt", sep="_"), row.names = FALSE)
  
  # Export also the confusion matrix on the resamples and on the test
  # change the positive class to be nLC and get precision in addition to sensitivity and specificity
  # specificity will be fraction of total negative (LC) samples that are well predicted
  # sensitivity will be fraction of total positive (nLC) samples that are well predicted = recall
  # precision will be fraction of samples predicted as the positive class (nLC) that are indeed of that class (nLC)
  
   cm2 <- confusionMatrix(data = trained_m$pred$pred, reference = trained_m$pred$ob, mode = "everything")
  capture.output(cm2, file=paste(m_rf_mt, "training_ConfMat.txt", sep="_"))
  
  
  cm3 <- confusionMatrix(data=as.factor(test_set_wPreds$pred_cl), reference=relevel(as.factor(as.character(test_set_wPreds$IUCN_TS_sum)), "nonLC"), mode = "everything")
  capture.output(cm3, file=paste(m_rf_mt, "testing_ConfMat.txt", sep="_"))
  
  
}






### functions needed inside iucn_model_rf_mt
# Copied from https://topepo.github.io/caret/using-your-own-model-in-train.html#illustrative-example-5-optimizing-probability-thresholds-for-class-imbalances
# Adapted to tune mtry too
# I changed a bit the comments


# Initialise the new method
# Copy the model code from the original random forest method:
thresh_code <- getModelInfo("rf", regex = FALSE)[[1]]
thresh_code$type <- c("Classification")

# Add the threshold as another tuning parameter to the existing one
thresh_code$parameters <- data.frame(parameter = c("mtry", "threshold"),
                                     class = c("numeric", "numeric"),
                                     label = c("#Randomly Selected Predictors",
                                               "Probability Cutoff"))

# Set new tuning grid options for the method
# if search set to grid (the default), set the grid to cover all combinations of the mtry passed by opt_mtry, which will be estimated from the number of variables 
# (--> need to load the function inside iucn_model_t_rf_mt, not here)
# and the number of thresholds ("len"), which will be passed via tuneLength in train
# if search is set to something else, do same for the mtry but something I don't understand for the thresholds
# (we will keep search set to grid anyway)


thresh_code$grid <- function(x, y, len = NULL, search = "grid", mtry_fixed=opt_mtry) {
  if(search == "grid") {
    grid <- expand.grid(mtry = mtry_fixed, threshold = seq(.01, .99, length = len))
  } else {
    grid <- expand.grid(mtry = mtry_fixed, threshold = runif(runif, min = 0, max = 1))
  }
  grid
}

# Based on grid, build a list of submodels, one submodel per mtry
# for each submodel, all threshold values are listed

thresh_code$loop = function(grid) {   
  library(plyr)
  loop <- ddply(grid, c("mtry"), function(x) c(threshold = max(x$threshold)))
  submodels <- vector(mode = "list", length = nrow(loop))
  for(i in seq(along = loop$threshold)) {
    index <- which(grid$mtry == loop$mtry[i])
    cuts <- grid[index, "threshold"] 
    submodels[[i]] <- data.frame(threshold = cuts[cuts != loop$threshold[i]])
  }    
  list(loop = loop, submodels = submodels)
}

# Fit the model independently of the threshold parameter
# implement a check because would not work when more than 2 classes

thresh_code$fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
  if(length(levels(y)) != 2)
    stop("This works only for 2-class problems")
  randomForest(x, y, mtry = param$mtry, ...)
}


# Set the predict method so that get a single probability prediction and use the chosen threshold to get the predicted class, 
# instead of the default 0.5

thresh_code$predict = function(modelFit, newdata, submodels = NULL) {
  class1Prob <- predict(modelFit, newdata, type = "prob")[, modelFit$obsLevels[1]]
  out <- ifelse(class1Prob >= modelFit$tuneValue$threshold,modelFit$obsLevels[1],modelFit$obsLevels[2])

    if(!is.null(submodels)) {
    tmp2 <- out
    out <- vector(mode = "list", length = length(submodels$threshold))
    out[[1]] <- tmp2
    
    for(i in seq(along = submodels$threshold)) {
      out[[i+1]] <- ifelse(class1Prob >= submodels$threshold[[i]],
                           modelFit$obsLevels[1], 
                           modelFit$obsLevels[2])
    }
  } 
  out  
}



# The probabilities are always the same but we have to create mulitple versions of the method$prob to evaluate the data 
# across thresholds
thresh_code$prob = function(modelFit, newdata, submodels = NULL) {
  out <- as.data.frame(predict(modelFit, newdata, type = "prob"))
  if(!is.null(submodels)) {
    probs <- out
    out <- vector(mode = "list", length = length(submodels$threshold)+1)
    out <- lapply(out, function(x) probs)
  } 
  out 
}



## design new summary function that calculates the distance "DPM" (distance to perfect model) = distance between the model and the best possible model
# The latter would have both Sens and Spe = 1
# so DPM represents how different are sens and spec for the current candidate value of the probability threshold

fourStats <- function (data, lev = levels(data$obs), model = NULL) {
  
  out <- c(twoClassSummary(data, lev = levels(data$obs), model = NULL)) # gets the model summary as normal, i.e. retrieves ROC, sens and spe
  coords <- matrix(c(1, 1, out["Spec"], out["Sens"]), ncol = 2, byrow = TRUE) # build a matrix with our sens and spe and the best model sens and spe (ie 1 and 1)
  colnames(coords) <- c("Spec", "Sens")
  rownames(coords) <- c("Best", "Current")
  c(out, DPM = dist(coords)[1]) # calculate the distance between the values of the matrix
}













