climate_ts <- function(unit.raster, pls_data, now_rast, then_rast){
    climate_ann <- read.csv('data/input/annualclim.csv')[,-1]
  
    #  This is the code used to change t_mean to t_diff:
#     climate_ann$value[climate_ann$variable == 'tmn'] <- 
#       climate_ann$value[climate_ann$variable == 'tmx'] -
#       climate_ann$value[climate_ann$variable == 'tmi'] 
#     climate_ann$variable <- as.character(climate_ann$variable)
#     climate_ann$variable[climate_ann$variable == 'tmn'] <- 'tdf'
#     write.csv(climate_ann, 'data/input/annualclim.csv')
  #  load the PDSI data as a function so the memory use is contained in the function.
  #  Load the PDSI data for the region:
  #  download from - http://www.ncdc.noaa.gov/paleo/pdsi.html
  #  We shortcut this if the pdsi_df file exists:
  if ('pdsi_df.rds' %in% list.files('data/output/')) {
    pdsi_df <- readRDS('data/output/pdsi_df.rds')  
  } else {
    pdsi_out <- function(unit.raster, pls_data) {
      NCDF <- stack('data/input/NADAv2-2008.nc')
      
      ncdf.resamp <- raster::resample(raster::projectRaster(NCDF, 
                                                            crs = '+init=epsg:3175'), 
                                      y = unit.raster)
      
      pdsi <- raster::getValues(ncdf.resamp)[pls_data$cell,]
      
      pdsi_df <- data.frame(variable = 'pdsi', 
                            value    = colMeans(pdsi, na.rm = TRUE)[1:206],
                            year     = 2006:1801)
      pdsi_df
    }
    
    pdsi_df <- pdsi_out(unit.raster, pls_data)
    saveRDS(pdsi_df, 'data/output/pdsi_df.rds')
  }
  
  climate_ann <- rbind(climate_ann, pdsi_df)
  
  # This removes any raster cells without PLS data to get a regional climate
  #  picture for the plotting:
  clim_pattern <- function(rast_in, pls_rast, year){
  
    rast_mask <- rast_in
    rast_mask[!(1:ncell(rast_in)) %in% pls_data$cell] <- NA
  
    data.frame(variable = c('P[ann]', 'T[max]', 'T[diff]', 'T[min]'),
               value    = colMeans(getValues(rast_mask), na.rm = TRUE),
               year     = year)
  }
  
  modern <- clim_pattern(now_rast, pls_rast, 1995)
  past   <- clim_pattern(then_rast, pls_rast, 1910)
  
  climate_ann$variable <- factor(climate_ann$variable, 
                             levels = c('tmx', 'tdf', 'tmi', 'ppt', 'pdsi'),
                             labels = c('T[max]', 'T[diff]', 'T[min]', 'P[ann]', 'PDSI'))
  
  clim_plot <- ggplot(data = climate_ann, aes(x = year, y = value)) + geom_line(size = 1.3) +
    geom_vline(data = modern, aes(xintercept = year), size = 4, color = 'red') +
    geom_vline(data = past,   aes(xintercept = year), size = 4, color = 'red') +
    facet_grid(variable~.,    scales = 'free_y', labeller = label_parsed) +
    ylab('') +
    xlab('Year (CE)') +
    coord_cartesian(c(1801, 2015)) +
    theme_bw() +
    theme(axis.text = element_text(family = 'serif', size = 12),
          axis.title = element_text(family = 'serif', face = 'bold', size = 14),
          strip.text = element_text(family = 'serif', face = 'bold', size = 16),
          legend.position = "none")
    
  ggsave(clim_plot, filename = 'figures/clim_plot.tiff', width = 8, height = 6, dpi = 150)
  list(modern = modern, past = past, climate = climate_ann, plot = clim_plot)
}