# create response curves in ggplot
ggResponse <- function(models, covariates, colPlot=3, responseName="Prediction", ...){
  require(raster)
  require(blockCV)
  require(dplyr)
  require(reshape)
  n <- 0
  categoricals <- c()
  if(is(covariates, "Raster")){
    nlayer <- raster::nlayers(covariates)
    meanVars <- matrix(nrow=100, ncol=nlayer)
    meanVars <- as.data.frame(meanVars)
    names(meanVars) <- names(covariates)
    ranges <- predictions <- meanVars
    categoricals <- names(covariates)[which(is.factor(covariates))]
    if(anyNA(maxValue(covariates))){
      naminmax <- which(is.na(maxValue(covariates)))
      for(l in naminmax){
        covariates[[l]] <- raster::setMinMax(covariates[[l]])
      }
    }
    # calculate the means and ranges for non-categorical vars
    for(i in 1:nlayer){
      if(!is.factor(covariates[[i]])){
        ranges[,i] <- seq(minValue(covariates[[i]]), maxValue(covariates[[i]]), length.out = 100)
        meanVars[,i] <- rep(mean(values(covariates[[i]]), na.rm=TRUE), 100)
      }
    }
  } else if(is(covariates, "data.frame")){
    nlayer <- ncol(covariates)
    meanVars <- matrix(nrow=100, ncol=nlayer)
    meanVars <- as.data.frame(meanVars)
    names(meanVars) <- names(covariates)
    ranges <- predictions <- meanVars
    for(b in 1:nlayer){
      if(is.factor(covariates[,b])){
        n <- n + 1
        categoricals[n] <- names(covariates)[b]
      }
    }
    # calculate the means and ranges for non-categorical vars
    for(i in 1:nlayer){
      if(!is.factor(covariates[[i]])){
        ranges[,i] <- seq(min(covariates[,i]), max(covariates[,i]), length.out = 100)
        meanVars[,i] <- rep(mean(covariates[,i]), 100)
      }
    }
  } else{
    stop("covariates should be a raster layer or data.frame object contining variables used in the model")
  }
  # calculate the means and ranges for categorical vars
  if(length(categoricals) > 0){
    cats <- which(names(covariates) %in% categoricals) # categorical vars
    if(is(covariates, "data.frame")){
      for(ct in cats){
        commCats <- names(which(table(covariates[,ct]) == max(table(covariates[,ct]))))
        level <- unlist(levels(covariates[,ct]))
        ranges[,ct] <- c(level, sample(level, 100 - length(level), replace = T))
        meanVars[,ct] <- rep(commCats, 100)
        ranges[,ct] <- as.factor(ranges[,ct])
        meanVars[,ct] <- as.factor(meanVars[,ct])
      }
    } else{
      for(ct in cats){
        commCats <- names(which(table(values(covariates[[ct]])) == max(table(values(covariates[[ct]])))))
        level <- unlist(levels(covariates[[ct]]))
        ranges[,ct] <- c(level, sample(level, 100 - length(level), replace = T))
        meanVars[,ct] <- rep(commCats, 100)
        ranges[,ct] <- as.factor(ranges[,ct])
        meanVars[,ct] <- as.factor(meanVars[,ct])
      }
    }
  }
  # predict with the model
  for(j in 1:nlayer){
    mydf <- cbind(ranges[,j], meanVars[,-j])
    names(mydf)[1] <- colnames(meanVars)[j]
    for(c in categoricals){
      if(is(covariates, "Raster")){
        levels(mydf[, c]) <- as.character(unlist(levels(covariates[[c]])))
      } else{
        levels(mydf[, c]) <- levels(covariates[,c])
      }
    }
    predictions[,j] <- predict(models, mydf, ...)
  }
  # change the cats to numeric for melting
  if(length(categoricals) > 0){
    for(ct in cats){
      ranges[,ct] <- as.numeric(as.character(ranges[,ct]))
      predictions[,ct] <- as.numeric(as.character(predictions[,ct]))
    }
  }
  val <- reshape::melt(ranges)
  # nrow(val);head(val, 10)
  prd <- reshape::melt(predictions)
  # nrow(prd); head(prd, 10)
  finaltable <- dplyr::bind_cols(val, prd)
  yMin <- min(finaltable$value1)
  yMax <- max(finaltable$value1)
  # create the plots
  pp <- list()
  for(k in 1:nlayer){
    down <- k*100-99
    up <- k*100
    pp[[k]] <- ggplot(data=finaltable[down:up,], aes(x=value, y=value1)) + geom_line() + 
      xlab(stringr::str_to_sentence(names(covariates)[k])) + scale_y_continuous(name=responseName, limits = c(yMin, yMax)) +
      theme_bw()
  }
  # create plot for categorical variables
  if(length(categoricals) > 0){
    for(ct in cats){
      lutable <- finaltable[which(finaltable$variable == names(covariates)[ct]),]
      lutable$value <- as.factor(lutable$value)
      catText <- "ggplot(data=lutable, aes(x=value, y=value1)) + geom_point(size=0.1) + 
      xlab(stringr::str_to_sentence(names(covariates)[ct])) + scale_y_continuous(name=responseName, limits = c(yMin, yMax)) + theme_bw() +"
      for(i in 1:length(level)){
        n <- n + 1 
        if(i < length(level)){
          tmp <- sprintf("geom_segment(aes(x=(%d - 0.5),xend=(%d + 0.5),y=lutable[%d,'value1'],yend=lutable[%d,'value1']), size=1.2) +", i,i,i,i)
        } else{
          tmp <- sprintf("geom_segment(aes(x=(%d - 0.5),xend=(%d + 0.5),y=lutable[%d,'value1'],yend=lutable[%d,'value1']), size=1.2)", i,i,i,i)
        }
        catText <- paste(catText, tmp)
      }
      pp[[ct]] <- eval(parse(text = catText))
    }
  }
  # create final plot
  responsePlot <- "blockCV:::multiplot("
  for(i in 1:length(pp)){
    if(i < length(pp)){
      tmp <- sprintf("pp[[%d]], ", i)
    } else{
      tmp <- sprintf("pp[[%d]], cols = %d)", i, colPlot)
    }
    responsePlot <- paste0(responsePlot, tmp)
  }
  # plot final responses
  print(eval(parse(text = responsePlot)))
  return(pp)
}
