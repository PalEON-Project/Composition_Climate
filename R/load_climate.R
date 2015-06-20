#  Load the climate data:
#  The large prism datasets live on my hard drive, but are accessible through:
#

#  To crop and reproject we use the R file, included hereL
#  'R/building_small_climate.R'

load_climate <- function(){
  clim.files <- list.files('data/input/climate', full.names = TRUE)

  now.rast <- stack(clim.files[regexpr('norm', clim.files)>0])
  then.rast <- stack(clim.files[regexpr('early', clim.files)>0])

  bio.vals <- c('ppt', 'tmax', 'tmean', 'tmin')
  return(list(now = now.rast, then = then.rast, names = bio.vals))
}

climate <- load_climate()
