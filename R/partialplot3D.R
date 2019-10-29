# create partial plot 3D
response3D <- function(models, covariates, vars = c(1, 2), theta = 130, 
                       phi = 40, xlab = NULL, ylab = NULL, zlim = NULL, 
                       legend = FALSE, label = FALSE, col = "gray", 
                       responseName = "Prediction", ...){
  require(raster)
  require(plot3D)
  # require(dplyr)
  n <- 0
  categoricals <- c()
  if(is(covariates, "Raster")){
    nlayer <- raster::nlayers(covariates)
    meanVars <- matrix(nrow=30, ncol=nlayer)
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
        ranges[,i] <- seq(minValue(covariates[[i]]), maxValue(covariates[[i]]), length.out = 30)
        meanVars[,i] <- rep(mean(values(covariates[[i]]), na.rm=TRUE), 30)
      }
    }
  } else if(is(covariates, "data.frame")){
    nlayer <- ncol(covariates)
    meanVars <- matrix(nrow=30, ncol=nlayer)
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
        ranges[,i] <- seq(min(covariates[,i]), max(covariates[,i]), length.out = 30)
        meanVars[,i] <- rep(mean(covariates[,i]), 30)
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
        ranges[,ct] <- c(level, sample(level, 30 - length(level), replace = T))
        meanVars[,ct] <- rep(commCats, 30)
        ranges[,ct] <- as.factor(ranges[,ct])
        meanVars[,ct] <- as.factor(meanVars[,ct])
      }
    } else{
      for(ct in cats){
        commCats <- names(which(table(values(covariates[[ct]])) == max(table(values(covariates[[ct]])))))
        level <- unlist(levels(covariates[[ct]]))
        ranges[,ct] <- c(level, sample(level, 30 - length(level), replace = T))
        meanVars[,ct] <- rep(commCats, 30)
        ranges[,ct] <- as.factor(ranges[,ct])
        meanVars[,ct] <- as.factor(meanVars[,ct])
      }
    }
  }
  # change the cats to numeric for building stack
  if(length(categoricals) > 0){
    for(i in 1:length(categoricals)){
      ct <- categoricals[i]
      ranges[,ct] <- as.numeric(as.character(ranges[,ct]))
      meanVars[,ct] <- as.numeric(as.character(meanVars[,ct]))
    }
  }
  plot_res <- 30
  plot_col <- 30
  plot_row <- 30
  x <- matrix(rep(ranges[,vars[1]], plot_res), nrow = plot_row, ncol = plot_col, byrow = TRUE)
  y <- matrix(rep(ranges[,vars[2]], plot_res), nrow = plot_row, ncol = plot_col, byrow = FALSE)
  mat13D <- raster::raster(matrix(rep(ranges[,vars[1]], plot_res), nrow = plot_row, ncol = plot_col, byrow = TRUE))
  mat23D <- raster::raster(matrix(rep(ranges[,vars[2]], plot_res), nrow = plot_row, ncol = plot_col, byrow = FALSE))
  # condition on categorical covariates
  if(names(covariates)[vars[1]] %in% categoricals){
    plot_col <- plot_res <- length(unique(covariates[,vars[1]]))
    x <- matrix(as.numeric(levels(covariates[,vars[1]])), nrow = plot_row, ncol = plot_col, byrow = TRUE)
    y <- matrix(rep(ranges[,vars[2]], plot_res), nrow = plot_row, ncol = plot_col, byrow = FALSE)
    yc <- as.numeric(levels(covariates[,vars[1]]))
    xc <- ranges[,vars[2]]
    mat23D <- raster::raster(matrix(rep(ranges[,vars[2]], plot_res), nrow = plot_row, ncol = plot_col, byrow = FALSE))
    mat13D <- raster::raster(x)
    direct <- "x"
  }
  if(names(covariates)[vars[2]] %in% categoricals){
    plot_row <- plot_res <- length(unique(covariates[,vars[2]]))
    x <- matrix(rep(ranges[,vars[1]], plot_res), nrow = plot_row, ncol = plot_col, byrow = TRUE)
    y <- matrix(as.numeric(levels(covariates[,vars[2]])), nrow = plot_row, ncol = plot_col, byrow = FALSE)
    xc <- as.numeric(levels(covariates[,vars[2]]))
    yc <- ranges[,vars[1]]
    mat13D <- raster::raster(matrix(rep(ranges[,vars[1]], plot_res), nrow = plot_row, ncol = plot_col, byrow = TRUE))
    mat23D <- raster::raster(y)
    direct <- "y"
  }
  # vars <- vars
  stack3D <- setNames(raster::stack(mat13D, mat23D), names(covariates)[vars])
  # plot(stack3D)
  rest_of_layers <- names(covariates)[-vars]
  for(m in rest_of_layers){
    print(m)
    stack3D <- raster::stack(stack3D,
                             setNames(raster::raster(matrix(meanVars[1,m], nrow = plot_row, ncol = plot_col)), m)
    )
  }
  # changing back to categoricals
  if(length(categoricals) > 0){
    for(i in 1:length(categoricals)){
      ct <- categoricals[i]
      stack3D[[ct]] <- raster::ratify(stack3D[[ct]])
    }
  }
  prediction <- raster::predict(stack3D, models, ...)
  if(is.null(zlim)){
    zlim <- c(raster::minValue(prediction), raster::maxValue(prediction))
  }
  if(label == TRUE){
    lab <- "detailed"
  } else{
    lab <- "simple"
  }
  if(is.null(xlab)){
    xlab <- names(covariates)[vars[1]]
  }
  if(is.null(ylab)){
    ylab <- names(covariates)[vars[2]]
  }
  if(any(names(covariates)[vars] %in% categoricals)){
    if(is.null(xlab)){
      xlab <- names(covariates)[vars[2]]
    }
    if(is.null(ylab)){
      ylab <- names(covariates)[vars[1]]
    }
    plot3D::hist3D(z = as.matrix(prediction),
                   y = yc, 
                   x = xc, 
                   border = "black",
                   theta = theta, 
                   phi = phi,
                   shade = 0.5, 
                   space = 0.1,
                   col = col,
                   box = TRUE, 
                   facets = TRUE,
                   zlim = zlim,
                   inttype = 1,
                   bty = "b", 
                   ticktype = lab,
                   colkey = legend,
                   ltheta = 180, 
                   lphi = 160,
                   xlab = xlab,
                   ylab = ylab,
                   zlab = responseName) 
    # add ribbon
    plot3D::ribbon3D(z = as.matrix(prediction),
                     y = yc,
                     x = xc,
                     along = direct,
                     theta = theta, 
                     phi = phi,
                     shade = 0.5, 
                     curtain = FALSE,
                     space = 0.1,
                     col = col,
                     box = TRUE, 
                     border = "black",
                     add = TRUE,
                     facets = TRUE,
                     zlim = zlim,
                     inttype = 1,
                     bty = "b", 
                     ticktype = lab,
                     colkey = legend,
                     ltheta = 180, 
                     lphi = 160,
                     xlab = xlab,
                     ylab = ylab,
                     zlab = responseName)
  } else{
    plot3D::surf3D(x = x, y = y, z = as.matrix(prediction), 
                   theta = theta, 
                   phi = phi,
                   shade = 0.5, 
                   box = TRUE, 
                   border = "black", 
                   zlim = zlim,
                   inttype = 1,
                   bty = "b", 
                   curtain = TRUE,
                   # colvar = colvar,
                   col = col, 
                   ticktype = lab,
                   colkey = legend,
                   ltheta = 180, 
                   lphi = 160,
                   xlab = xlab,
                   ylab = ylab,
                   zlab = responseName)
  }
}
