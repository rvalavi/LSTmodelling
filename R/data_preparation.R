# data preparation
library(tidyverse)
library(raster)
library(sf)
library(viridis)

# loading raster data
covs <- raster::brick("data/covariates.grd")
plot(covs)

# generate a random sample
set.seed(761)
rp <- dismo::randomPoints(covs[[1]], 490) %>% # get some samples form water bodies
  rbind(xyFromCell(covs$LandUse, sample(which(values(covs$LandUse) == 2), 10)))
# convert to sf object
resp <- st_as_sf(as.data.frame(rp), coords = c("x", "y"), crs = st_crs(32639))
plot(covs$LST)
points(rp)

# extract the predictor variables for the sample sites
train_data <- raster::extract(covs, rp) %>%
  as.data.frame()

# convert lu to factor
train_data$LandUse <- as.factor(train_data$LandUse)
# merge the factor levels
train_data$LandUse <- fct_recode(train_data$LandUse, "1" = "4", "4" = "5", "5" = "6")
levels(train_data$LandUse)

anyNA(train_data) # check for any NA values


# plot the data
ggplot(aes(x = Road, y = NDVI, col = LST), data = train_data) +
  geom_point(size = 2, alpha = .8) +
  scale_color_viridis(option = "A", direction = 1) +
  theme_bw() +
  theme(axis.text = element_text(size = 11)) +
  xlab("Distance to road")
