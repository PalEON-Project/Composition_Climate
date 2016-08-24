#  Setup some data and libraries:
#install.packages(c('colorspace', 'raster', 'ggplot2', 'mgcv', 'reshape2', 'plyr', 
#                   'gridExtra', 'rgdal', 'RColorBrewer', 'analogue', 'maptools',
#                   'captioner'))

library(colorspace)
library(captioner)
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

fig_nums <- captioner(prefix = "Figure")
tab_nums <- captioner(prefix = "Table")

base_rast <- raster("data/albers/density//Abies_balsamea_density_alb.tif")
base_rast <- setValues(base_rast, 1:ncell(base_rast))

pls_rast <- raster("data/input//dens_v0.9-4.tif")
unit_raster <- setValues(pls_rast, 1:ncell(pls_rast))

land_use <- raster("data/input/nlcd_cropped_repro.tif")

#  We need to change the raster to a data.frame:
land_use <- data.frame(cell = extract(unit_raster,
                                      xyFromCell(land_use, 1:ncell(land_use))),
                       xyFromCell(land_use, 1:ncell(land_use)),
                       value = getValues(land_use))

reassign <- read.csv("data/input/landuse_reassignment.csv",
                     stringsAsFactors = FALSE)
land_use$class <- factor(reassign[match(land_use$value, reassign[,1]),2])

model_proj <- "+init=epsg:3175"

usa <- readOGR("data/input/shapes/usa/us.shp", "us")
canada <- readOGR("data/input/shapes/canada/PROVINCE.SHP", "PROVINCE")