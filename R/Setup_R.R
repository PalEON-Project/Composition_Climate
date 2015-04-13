#  Setup some data and libraries:
library(raster)
library(ggplot2)
library(mgcv)
library(reshape2)
library(plyr)
library(ncdf4)
library(gridExtra)

base.rast <- raster('data/albers/density//Abies_balsamea_density_alb.tif')
base.rast <- setValues(base.rast, 1:ncell(base.rast))

pls.rast <- raster('data/input//dens_v0.9-4.tif')
unit.raster <- setValues(pls.rast, 1:ncell(pls.rast))
