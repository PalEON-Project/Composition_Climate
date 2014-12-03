#  Setup some data and libraries:
library(raster)
library(ggplot2)
library(mgcv)
library(reshape2)
library(plyr)

base.rast <- raster('data/albers/density//Abies_balsamea_density_alb.tif')
base.rast <- setValues(base.rast, NA)
