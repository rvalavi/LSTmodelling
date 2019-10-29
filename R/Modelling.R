# k-fold cross-validation
library(dismo)

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


# Create the R^2 string for printing
rf_tp <- bquote(R^2 == .(round(summary(lm(LST ~ rf_cv, predictions))$r.squared, 3)))
rf_tp2 <- bquote(RMSE == .(round(RMSE(predictions$LST, predictions$rf_cv), 3)))
# Change the mypanel function to use the lm1 object
rf_panel <- function(x,y,...){
  panel.hexbinplot(x, y, ...)
  panel.text(30,44, labels = rf_tp)
  panel.text(30.7,42.66, labels = rf_tp2)
}
# hexplot for the correlation of obs vs. pred
hex_rf <- hexbinplot(LST ~ rf_cv,
                     data = predictions,
                     aspect = 1,
                     xbins = 20,
                     type = "r",
                     cex.labels = .9,
                     cex.title = 1.0,
                     xaxis.cex = 1,
                     yaxis.cex = 1,
                     xaxis.fontface = 1,
                     yaxis.fontface = 1,
                     xlab.cex = 1.5,
                     ylab.cex = 1.5,
                     maxcnt = 17,
                     panel = rf_panel,
                     colramp = function(n) rev(terrain.colors(16)),
                     xlab = "Predicted LST",
                     ylab = "Landsat LST",
                     main = "RF"
)
hex_rf

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


# Create the R^2 string for printing
brt_tp <- bquote(R^2 == .(round(summary(lm(LST ~ brt_cv, predictions))$r.squared, 3)))
brt_tp2 <- bquote(RMSE == .(round(RMSE(predictions$LST, predictions$brt_cv), 3)))
# Change the mypanel function to use the lm1 object
brt_panel <- function(x,y,...){
  panel.hexbinplot(x, y, ...)
  panel.text(30,44, labels = brt_tp)
  panel.text(30.7,42.66, labels = brt_tp2)
}
# hexplot for the correlation of obs vs. pred
hex_brt <- hexbin::hexbinplot(LST ~ brt_cv,
                              data = predictions,
                              aspect = 1,
                              xbins = 20,
                              type = "r",
                              cex.labels = .9,
                              cex.title = 1.0,
                              xaxis.cex = 1,
                              yaxis.cex = 1,
                              xaxis.fontface = 1,
                              yaxis.fontface = 1,
                              xlab.cex = 1.5,
                              ylab.cex = 1.5,
                              maxcnt = 17,
                              panel = brt_panel,
                              colramp = function(n) rev(terrain.colors(16)),
                              xlab = "Predicted LST",
                              ylab = "Landsat LST",
                              main = "BRT"
)
hex_brt

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


# Create the R^2 string for printing
gam_tp <- bquote(R^2 == .(round(summary(lm(LST ~ gam_cv, predictions))$r.squared, 3)))
gam_tp2 <- bquote(RMSE == .(round(RMSE(predictions$LST, predictions$gam_cv), 3)))
# Change the mypanel function to use the lm1 object
gam_panel <- function(x,y,...){
  panel.hexbinplot(x, y, ...)
  panel.text(30,44, labels = gam_tp)
  panel.text(30.7,42.66, labels = gam_tp2)
}
# hexplot for the correlation of obs vs. pred
hex_gam <- hexbinplot(LST ~ gam_cv,
                      data = predictions,
                      aspect = 1,
                      xbins = 20,
                      type = "r",
                      cex.labels = .9,
                      cex.title = 1.0,
                      xaxis.cex = 1,
                      yaxis.cex = 1,
                      xaxis.fontface = 1,
                      yaxis.fontface = 1,
                      xlab.cex = 1.5,
                      ylab.cex = 1.5,
                      maxcnt = 17,
                      panel = gam_panel,
                      colramp = function(n) rev(terrain.colors(16)),
                      xlab = "Predicted LST",
                      ylab = "Landsat LST",
                      main = "GAM"
)
hex_gam

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


# Create the R^2 string for printing
svm_tp <- bquote(R^2 == .(round(summary(lm(LST ~ svm_cv, predictions))$r.squared, 3)))
svm_tp2 <- bquote(RMSE == .(round(RMSE(predictions$LST, predictions$svm_cv), 3)))
# Change the mypanel function to use the lm1 object
svm_panel <- function(x,y,...){
  panel.hexbinplot(x, y, ...)
  panel.text(30,44, labels = svm_tp)
  panel.text(30.7,42.66, labels = svm_tp2)
}
# hexplot for the correlation of obs vs. pred
hex_svm <- hexbinplot(LST ~ svm_cv,
                      data = predictions,
                      aspect = 1,
                      xbins = 20,
                      type = "r",
                      cex.labels = 0.9,
                      cex.title = 1.0,
                      xaxis.cex = 1,
                      yaxis.cex = 1,
                      xaxis.fontface = 1,
                      yaxis.fontface = 1,
                      xlab.cex = 1.5,
                      ylab.cex = 1.5,
                      maxcnt = 17,
                      panel = svm_panel,
                      colramp = function(n) rev(terrain.colors(16)),
                      xlab = "Predicted LST",
                      ylab = "Landsat LST",
                      main = "SVM"
)
hex_svm


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
