#  This file generates the gridded climate data used for the analysis.
#  We take gridded PRISM data, at a 4 x 4km resolution, obtained from PRISM:
#  http://www.prism.oregonstate.edu/
#
#  From there we have two datasets, a modern 'normal' dataset used for the FIA data
#  and a historical data, a normal we generate for 1895 - 1924.  We extract the data and
#  then crop, reproject, and resample so that it meets the same resolution as the base
#  PalEON dataset.

library(raster)
library(snowfall)

sfInit(parallel = TRUE, cpus = 3)

paleon.grid <- raster('data/input/dens_v0.9-4.tif')

regrid <- function(i, class){
  
  variables <- c('ppt', 'tmax', 'tmin', 'tmean')
  
  if(class == 'norm'){
    all.files <- list.files(path = paste0('../../GriddedClimate/prism/PRISM_',variables[i], 
                     '_30yr_normal_4kmM2_all_bil/'), full.names=TRUE)
    all.files <- all.files[regexpr('\\d{2}_bil.bil$', all.files)>0]
    
    norm <- crop(stack(all.files), extent(c(-99, -64, 35, 51)))
    
    
    if(i == 1){
      norm <- calc(norm, mean) * 12 # Annual Precip
    }
    
    if(i == 2){
      # Maximum annual minimum:
      norm <- setValues(norm[[1]], apply(getValues(norm), 1, max))
    }
    
    if(i == 3){
      # Minimum annual maximum:
      norm <- setValues(norm[[1]], apply(getValues(norm), 1, min))
    }
    
    if(i == 4){  
      norm <- calc(norm, mean)      # Mean Annual Temperature
    }
    
    
  } else {
    setfiles <- list.files(paste0('../../GriddedClimate/prism/PRISM_',variables[i], 
                       '_stable_4kmM2_189501_198012_bil'),full.names=TRUE)
    
    # This gives all monthly values, we're going to stack them and then
    #  average all of them:
    
    annual <- regexpr('_(190[0-9]|189[0-9]|191[0-9]|192[0-4])\\d{2}_bil.bil$', setfiles)
    norm <- stack(setfiles[annual>0])
    norm <- crop(norm, extent(c(-99, -64, 35, 51)))
    
    if(i == 1){
      norm <- calc(norm, mean) * 12 # Annual Precip
    }
    
    if(i == 2){
      # Maximum annual minimum:
      years <- rep(1:30, each = 12)
      norm <- do.call(stack, lapply(1:30, function(x){
                                sub <- subset(norm, (1:360)[years == x])
                                setValues(norm[[1]], apply(getValues(sub), 1, max))
                              }))
      norm <- calc(norm, mean)
    }
    
    if(i == 3){
      # Minimum annual maximum:
      years <- rep(1:30, each = 12)
      
      norm <- do.call(stack, lapply(1:30, function(x){
        sub <- subset(norm, (1:360)[years == x])
        setValues(norm[[1]], apply(getValues(sub), 1, min))
      }))
      norm <- calc(norm, mean)
    }

    if(i == 4){  
      norm <- calc(norm, mean)      # Mean Annual Temperature
    }

  }

  norm.out <- resample(projectRaster(norm, paleon.grid), paleon.grid)
  writeRaster(norm.out, 
              paste0('data/input/climate/', variables[i], '_', class, '.tif'),
              overwrite = TRUE)
}

#  Very slow, but only needs to be done once:
lapply(c('norm'), function(x) lapply(1:4, function(y) regrid(y, x)))
