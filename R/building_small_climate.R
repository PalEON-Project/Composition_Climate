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


#  We need to get a plot of each value through time for the domain:
variables <- c('ppt', 'tmax', 'tmin', 'tmean')

clim <- list()

for(i in 1:4){
  
  setfiles <- list.files(paste0('../../GriddedClimate/prism/PRISM_',variables[i], 
                                '_stable_4kmM2_189501_198012_bil'),full.names=TRUE)
  setfiles <- setfiles[regexpr('.bil$', setfiles)>0]
  
  # This gives all monthly values, we're going to stack them and then
  #  average all of them:
  
  clim[[i]] <- data.frame(param = variables[i],
                          val =   rep(NA, length(setfiles)),
                          year =  c(rep(c(1:12), 119), 1:10),
                          month = c(rep(1895:2013, each=12), rep(2014, 10)))
  
  for(j in 1:length(setfiles)){
    
    norm <- raster(setfiles[j])
    norm <- projectRaster(crop(norm, extent(c(-99, -64, 35, 51))), paleon.grid)
    norm[is.na(paleon.grid)] <- NA
    clim[[i]]$val[j] <- mean(getValues(norm), na.rm=TRUE)
  }    
  
}

clim.year <- data.frame(ppt = as.numeric(dcast(clim[[1]], param ~ month, fun.aggregate = sum, value.var = 'val', na.rm=TRUE)[,-1]),
                        tmx = as.numeric(dcast(clim[[2]], param ~ month, fun.aggregate = max, value.var = 'val', na.rm=TRUE)[,-1]),
                        tmi = as.numeric(dcast(clim[[3]], param ~ month, fun.aggregate = min, value.var = 'val', na.rm=TRUE)[,-1]),
                        tmn = as.numeric(dcast(clim[[4]], param ~ month, fun.aggregate = mean, value.var = 'val', na.rm=TRUE)[,-1]))

clim.yout <- data.frame(melt(clim.year), year = 1895:2014)
write.csv(clim.yout, 'data//input/annualclim.csv')
