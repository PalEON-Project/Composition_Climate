#  Setup some data and libraries:
library(colorspace)
library(raster)
library(ggplot2)
library(mgcv)
library(reshape2)
library(plyr)
library(gridExtra)
library(rgdal)
library(RColorBrewer)
library(analogue)
library(maptools)

base.rast <- raster('data/albers/density//Abies_balsamea_density_alb.tif')
base.rast <- setValues(base.rast, 1:ncell(base.rast))

pls.rast <- raster('data/input//dens_v0.9-4.tif')
unit.raster <- setValues(pls.rast, 1:ncell(pls.rast))

land_use <- raster('../../../Maps/nlcd_2011_landcover_2011_edition_2014_10_10//nlcd_cropped_repro.tif')

#  We need to change the raster to a data.frame:
land_use <- data.frame(cell = extract(unit.raster, xyFromCell(land_use, 1:ncell(land_use))),
                       xyFromCell(land_use, 1:ncell(land_use)),
                       value = getValues(land_use))
reassign <- read.csv('data/input/landuse_reassignment.csv', stringsAsFactors=FALSE)
land_use$class <- factor(reassign[match(land_use$value, reassign[,1]),2])

model.proj <- '+init=epsg:3175'

usa <- readOGR('data/input/shapes/usa/us.shp', 'us')
canada <- readOGR('data/input/shapes/canada/PROVINCE.SHP', 'PROVINCE')
