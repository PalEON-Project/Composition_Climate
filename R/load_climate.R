#  Load the climate data:
path.root <- getwd()

if(regexpr('goring', path.root)>0){
  path.gis <- 'C:/Users/goring/Documents/GriddedClimate/CIEE_data'
}
if(regexpr('Goring', path.root)>0){
  path.gis = 'C:/Users/Simon Goring/Documents/GriddedClimate/CIEE_data'
}

## climate normal layers:
#naming convention: biox_norm1 = 1900-1975; biox_norm2 = 1985-2010
#  
then.files = list.files(path.gis, pattern = "_norm1.asc", 
                        recursive=TRUE, full.names = TRUE)
now.files = list.files(path.gis, pattern = "_norm2.asc", 
                       recursive = TRUE, full.names = TRUE)

load.stack <- function(filenames){
  #  Loads a raster stack based on the filenames provided, then crops to
  #  the extent provided in 'crop.ext'
  loaded.stack <- stack(filenames)
  proj4string(loaded.stack) <- '+proj=longlat +ellps=WGS84'
  
  return(loaded.stack)
}  

#  For each year, by variable:
## year-specific climate layers
bioclim.files <- list.files(path.gis, pattern=".asc$", recursive = TRUE)
bioclim.files <- bioclim.files[!(regexpr('norm', bioclim.files)>0 | regexpr('elevstd', bioclim.files) > 0)]

#  Get just the bioclim variable names, not the full path.
bio.vals <- unique(substr(bioclim.files, 
                          regexpr('/', bioclim.files) + 1,
                          regexpr('_', bioclim.files)-1))

right.clim <- function(x)regexpr(x, bioclim.files) > 0

#  Load each of the bioclim variables.  Takes some time.  This could
#  be sped up a bit if we parallelize.
#  The result is a list of RasterBricks for each bioclim value, 
#  representing each year (in sequence) from 1901 - 2010.

if(!'then.rast.rds' %in% list.files(path='data/output')){
  then.rast <- lapply(bio.vals, function(x){
    load.stack(paste0(path.gis, '/', bioclim.files[right.clim(x)][1:30]))
  })
  
  then.rast <- lapply(then.rast, 
                      function(x)resample(projectRaster(x, crs=CRS('+init=epsg:3175')), 
                                          base.rast))

  saveRDS(then.rast, file = 'data/output/then.rast.rds')
} else {
  then.rast <- readRDS('data/output/then.rast.rds')
}

if(!'now.rast.rds' %in% list.files(path='data/output')){
  now.rast <- lapply(bio.vals, function(x){
        load.stack(paste0(path.gis, '/', bioclim.files[right.clim(x)][76:101]))
        })

  now.rast <- lapply(now.rast, function(x)resample(projectRaster(x, crs=CRS('+init=epsg:3175')), 
                                             base.rast))
  saveRDS(now.rast, file = 'data/output/now.rast.rds')
} else {
  now.rast <- readRDS('data/output/now.rast.RData')
}
