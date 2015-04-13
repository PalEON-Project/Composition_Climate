source('R/load_fia.R')

pls.data <- read.csv('data/plss_density_alb_v0.9-4.csv', row.names=1)
colnames(pls.data) <- tolower(colnames(pls.data))

colnames(agg.dens) <- tolower(gsub(' ', '.', colnames(agg.dens), fixed=TRUE))

colnames(agg.dens) <- gsub('/', '.', colnames(agg.dens), fixed=TRUE)

agg.dens$poplar <- agg.dens$poplar.tulip.poplar

best.taxa <- c('beech', 'ironwood', 'hemlock', 'fir', 'cedar.juniper', 
              'basswood', 'spruce', 'elm', 'ash', 'pine', 'oak', 
              'tamarack', 'maple', 'poplar', 'birch')

agg.dens <- agg.dens[,c('cell', best.taxa)]
pls.data <- pls.data[,c('cell', best.taxa)]

agg.dens[is.na(agg.dens)] <- 0

cell.xy <- data.frame(cell = 1:ncell(unit.raster),
                     xyFromCell(unit.raster, 1:ncell(unit.raster)))
