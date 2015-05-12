climate <- read.csv('data/input/annualclim.csv')[,-1]

#load the PDSI data:
#  Load the PDSI data for the region:

NCDF <- stack('data/NADApdsi04.nc')
ncdf.resamp <- resample(projectRaster(NCDF, crs = '+init=epsg:3175'), y = unit.raster)
pdsi <- getValues(ncdf.resamp)[pls.data$cell,]
pdsi_df <- data.frame(variable = 'pdsi', 
                      value = colMeans(pdsi, na.rm=TRUE)[1:205],
                      year = 2005:1801)

climate <- rbind(climate, pdsi_df)

now.mask <- now.rast
now.mask[is.na(pls.rast)] <- NA

then.mask <- then.rast
then.mask[is.na(pls.rast)] <- NA

modern <- colMeans(getValues(now.mask), na.rm=TRUE)
past   <- colMeans(getValues(then.mask), na.rm=TRUE)

modern <- data.frame(variable = c('ppt', 'tmx', 'tmn', 'tmi'),
                     value = colMeans(getValues(now.mask), na.rm=TRUE),
                     year = 1990)

past <- data.frame(variable = c('ppt', 'tmx', 'tmn', 'tmi'),
                     value = colMeans(getValues(then.mask), na.rm=TRUE),
                     year = 1905)

clim_plot <- ggplot(data = climate, aes(x = year, y = value)) + geom_line(size = 2) +
  geom_point(data = modern, aes(x = year, y = value), size = 4, color = 'red') +
  geom_point(data = past, aes(x = year, y = value), size = 4, color = 'red') +
  facet_wrap(~variable, ncol = 1, scales='free_y') +
  theme(axis.text = element_text(family = 'serif', size = 12),
        axis.title = element_text(family = 'serif', face = 'bold', size = 14),
        strip.text = element_text(family = 'serif', face='bold', size = 14),
        legend.position = "none") +
  theme_bw()
  
ggsave(clim_plot, filename = 'figures/clim_plot.tiff', width = 8, height = 6, dpi = 150)
