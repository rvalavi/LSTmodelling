# k-fold cross-validation
kf <- 10
set.seed(2732)
spfolds$foldID <- dismo::kfold(nrow(train_data), kf)

# function to calculate RMSE
RMSE <- function(observed, predicted){
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}


# Random Forest -----------------------------------------------------------
library(randomForest)

# create a table for the predictions
predictions <- mutate(resp, LST = train_data$LST)
models_rf <- list()
# cross-validation for RF
for(k in 1:kf){
  trainSet <- which(spfolds$foldID != k)
  testSet <- which(spfolds$foldID == k)
  mod_rf <- randomForest(LST ~ ., 
                         data = train_data[trainSet,], 
                         ntree = 500, 
                         importance = TRUE)
  predictions[testSet,"rf_cv"] <- predict(mod_rf, newdata = train_data[testSet,])
  models_rf[[k]] <- mod_rf
  print(k)
}


# BRT ---------------------------------------------------------------------
# Boosted regression trees
library(gbm)

for(k in 1:kf){
  trainSet <- which(spfolds$foldID != k)
  testSet <- which(spfolds$foldID == k)
  mod_brt <- gbm::gbm(formula = LST ~ .,
                      distribution = "gaussian",
                      data = train_data[trainSet, ],
                      n.trees = 10000,
                      shrinkage = 0.001,
                      bag.fraction = 0.75,
                      interaction.depth = 5,
                      train.fraction = 1.0,
                      n.minobsinnode = 2, 
                      cv.folds = 10,
                      n.cores = 8)
  bestTree <- gbm::gbm.perf(mod_brt)
  predictions[testSet,"brt_cv"] <- predict(mod_brt, train_data[testSet,], n.trees = bestTree)
  print(k)
}


# GAM ---------------------------------------------------------------------
# Genralised additive model
library(mgcv)

myform <- LST ~ s(NDVI) + s(DEM) + s(Slope) + s(Solar) + s(Road) + LandUse 
models_gam <- list()
for(k in 1:kf){
  trainSet <- which(spfolds$foldID != k)
  testSet <- which(spfolds$foldID == k)
  mod_gam <- mgcv::gam(as.formula(myform),
                       data = train_data[trainSet, ],
                       gaussian(link = "identity"),
                       select = FALSE)
  summary(mod_gam)
  models_gam[[k]] <- mod_gam
  plot(mod_gam, pages=1, rug=TRUE, shade=TRUE)
  predictions[testSet,"gam_cv"] <- predict(mod_gam, train_data[testSet,])
  print(k)
}


# SVM ---------------------------------------------------------------------
# Support vector machine
library(e1071)

models_svm <- list()
for(k in 1:kf){
  trainSet <- which(spfolds$foldID != k)
  testSet <- which(spfolds$foldID == k)
  mod_svm <- e1071::svm(LST ~ .,
                        data = train_data[trainSet, ],
                        kernel = "radial",
                        cost = 2,
                        # gamma = .1,
                        scale = TRUE)
  models_svm[[k]] <- mod_svm
  predictions[testSet,"svm_cv"] <- predict(mod_svm, train_data[testSet,])
  print(k)
}


# Full models for spatial prediction --------------------------------------
fulmod_rf <- randomForest(LST ~ ., 
                          data = train_data, 
                          ntree = 500, 
                          importance = TRUE)

fulmod_brt <- gbm::gbm(formula = LST ~ .,
                       distribution = "gaussian",
                       data = train_data,
                       n.trees = 10000,
                       shrinkage = 0.001,
                       bag.fraction = 0.75,
                       interaction.depth = 6,
                       train.fraction = 1.0,
                       n.minobsinnode = 2, 
                       cv.folds = 10,
                       n.cores = 8)

myform <- LST ~ s(NDVI) + s(DEM) + s(Slope) + s(Solar) + s(Road) + LandUse 
fulmod_gam <- mgcv::gam(as.formula(myform),
                        data = train_data,
                        gaussian(link = "identity"),
                        select = FALSE)

fulmod_svm <- e1071::svm(LST ~ .,
                         data = train_data,
                         kernel = "radial",
                         cost = 2,
                         # gamma = .1,
                         scale = TRUE)
