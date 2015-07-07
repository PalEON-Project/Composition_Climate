climate_ts <- function(unit.raster, pls_data, now_rast, then_rast){
    climate <- read.csv('data/input/annualclim.csv')[,-1]
  
  #  load the PDSI data as a function so the memory use is contained in the function.
  #  Load the PDSI data for the region:
  #  download from - http://www.ncdc.noaa.gov/paleo/pdsi.html
  #  We shortcut this if the pdsi_df file exists:
  if('pdsi_df.rds' %in% list.files('data/output/')){
    pdsi_df <- readRDS('data/output/pdsi_df.rds')  
  } else{
    pdsi_out <- function(unit.raster, pls_data){
      NCDF <- stack('data/input/NADAv2-2008.nc')
      ncdf.resamp <- resample(projectRaster(NCDF, crs = '+init=epsg:3175'), y = unit.raster)
      pdsi <- getValues(ncdf.resamp)[pls.data$cell,]
      pdsi_df <- data.frame(variable = 'pdsi', 
                            value = colMeans(pdsi, na.rm=TRUE)[1:206],
                            year = 2006:1801)
      pdsi_df
    }
    
    pdsi_df <- pdsi_out(unit.raster, pls.data)
    saveRDS(pdsi_df, 'data/output/pdsi_df.rds')
  }
  
  climate <- rbind(climate, pdsi_df)
  
  # This removes any raster cells without PLS data to get a regional climate
  #  picture for the plotting:
  clim_pattern <- function(rast.in, pls.rast, year){
  
    rast.mask <- rast.in
    rast.mask[is.na(pls.rast)] <- NA
  
    data.frame(variable = c('P[ann]', 'T[max]', 'T[mean]', 'T[min]'),
               value = colMeans(getValues(rast.mask), na.rm=TRUE),
               year = year)
  }
  
  modern <- clim_pattern(now_rast, pls.rast, 1995)
  past   <- clim_pattern(then_rast, pls.rast, 1910)
  
  climate$variable <- factor(climate$variable, 
                             levels = c('tmx', 'tmn', 'tmi', 'ppt', 'pdsi'),
                             labels = c('T[max]', 'T[mean]', 'T[min]', 'P[ann]', 'PDSI'))
  
  clim_plot <- ggplot(data = climate, aes(x = year, y = value)) + geom_line(size = 2) +
    geom_point(data = modern, aes(x = year, y = value), size = 4, color = 'red') +
    geom_point(data = past, aes(x = year, y = value), size = 4, color = 'red') +
    facet_grid(variable~., scales='free_y', labeller = label_parsed) +
    ylab('') +
    xlab('Year (CE)') +
    theme_bw() +
    theme(axis.text = element_text(family = 'serif', size = 12),
          axis.title = element_text(family = 'serif', face = 'bold', size = 14),
          strip.text = element_text(family = 'serif', face='bold', size = 14),
          legend.position = "none")
    
  ggsave(clim_plot, filename = 'figures/clim_plot.tiff', width = 8, height = 6, dpi = 150)
  list(modern = modern, past = past, climate = climate, plot = clim_plot)
}