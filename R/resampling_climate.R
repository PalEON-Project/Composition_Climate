#  This code takes a climate variable and the FIA & PLS data
#  and returns a table with the climate variable for each taxon
#  within each cell in each of the four eras (PLS/PLS, FIA/FIA 
#  and the cross values).

resample_climate <- function(climate, agg_dens, pls_data){
  agg_dens[,-1] <- tran(agg_dens[,-1], method = 'proportion')
  pls_data[,-1] <- tran(pls_data[,-1], method = 'proportion')
  
  then_rast <- resample(climate$then, unit.raster)
  now_rast <- resample(climate$now, unit.raster)
  
  #For each taxon sample:
  #  1.  The FIA data for one of the now_rast years for each of the climate variables.
  #  2.  One of the then_rast years & one of the western layers
  #  3.  Get vectors of values
  
  # These are the climate ranges:
  ranges <- data.frame(min = c(400,  15, 0, -30),
                       max = c(2200, 35, 16, 2))
  
  # To do the analysis we want to get the climate layers for the past:
  get_dens_era <- function(taxon, clim.var){
    
    #  Function resamples proportion data from the composition model for the region
    #  and then draws from the composition model and the climate data.
    #  The function outputs a table with six columns:
    #  data      - the proportion composition within the cell
    #  class     - the era of observation
    #  cell      - the literal cell # in the raster
    #  plss.clim - the climate extracted from the 1900 - 1930s era
    #  fia.clim  - the climate extracted from the 1970 - 2000s era
    #
    #  Each table has PLSS/pres (taxon present in cell at PLSS)
    #                 FIA/pres  (taxon present in cell at PLSS)
    #                 FIA/np    (taxon present in cell at PLSS but not in FIA)
    #taxon.vals <- ncvar_get(western, taxon, c(1,1,1), c(-1, -1, -1))
  
    #  bind all the variables together:  
    values <- data.frame(cell  = c(rep(pls_data$cell, 2), rep(agg_dens$cell, 2)),
                         x     = c(rep(xyFromCell(unit.raster, pls_data$cell)[,'x'],2),
                                   rep(xyFromCell(unit.raster, agg_dens$cell)[,'x'],2)),
                         y     = c(rep(xyFromCell(unit.raster, pls_data$cell)[,'y'],2),
                                   rep(xyFromCell(unit.raster, agg_dens$cell)[,'y'],2)),
                         data  = c(rep(as.numeric(pls_data[,taxon]), 2),
                                   rep(agg_dens[,taxon],2)),
                         base  = c(rep('PLSS', nrow(pls_data)), 
                                   rep('PLSS', nrow(pls_data)), 
                                   rep('FIA', nrow(agg_dens)),
                                   rep('FIA', nrow(agg_dens))),
                         c.ref = c(rep('PLSS', nrow(pls_data)), 
                                   rep('FIA', nrow(pls_data)), 
                                   rep('FIA', nrow(agg_dens)),
                                   rep('PLSS', nrow(agg_dens))),
                         stringsAsFactors = FALSE)
    
  #  The problem here is that both 
    
    plss.plss.clim <- getValues(then_rast)[subset(values, base %in% 'PLSS' & c.ref %in% 'PLSS')$cell, clim.var]
    plss.fia.clim  <- getValues(now_rast)[subset(values, base %in% 'PLSS' & c.ref %in% 'FIA')$cell, clim.var]
    
    fia.fia.clim  <- getValues(now_rast)[subset(values, base %in% 'FIA' & c.ref %in% 'FIA')$cell, clim.var]
    fia.plss.clim   <- getValues(then_rast)[subset(values, base %in% 'FIA' & c.ref %in% 'PLSS')$cell ,clim.var]
    
    values$clim <- c(plss.plss.clim, plss.fia.clim,
                     fia.fia.clim, fia.plss.clim)
    
    values
    
  }
  
  #  This gives us samples from the data for both PLSS and FIA: 
  data.tables <-   llply(colnames(pls_data)[-1],
                          function(x){
                            lapply(1:4, function(y){
                              do.call(rbind.data.frame, 
                                lapply(1, function(z){
                                  aa <- try(get_dens_era(x, y))
                                  if(class(aa) == 'try-error'){
                                    aa <- data.frame(data=NA, base=NA, c.ref = NA, cell=NA, pres=NA, climate = NA)
                                  }
                                  aa}))
                            })
})
  
  data.out <- lapply(1:15, function(x)lapply(1:4, function(y){
    aa <- data.tables[[x]][[y]]
    aa$taxon <- colnames(pls_data)[-1][x]
    aa$climate <- climate$names[y]
    aa}))
  
  vegclim.table <- do.call(rbind.data.frame, lapply(data.out, function(x) do.call(rbind.data.frame, x)))
  vegclim.table <- vegclim.table[!is.na(vegclim.table$data),]
  
  vegclim.table$comb <- factor(with(vegclim.table, paste0(base, c.ref)), 
                               levels = c('PLSSPLSS','FIAFIA', 'PLSSFIA', 'FIAPLSS'),
                               labels = c('PLS Era', 'FIA Era', 'NA1', 'NA2'))
  
  vegclim.table$taxon <- factor(vegclim.table$taxon, 
                                levels = c('tamarack', 'pine', 'spruce', 'fir', 'hemlock', 'cedar.juniper',
                                           'poplar', 'maple', 'birch', 'beech', 'ironwood', 'basswood', 'ash', 'elm', 'oak'),
                                labels = c('Tamarack', 'Pine', 'Spruce', 'Fir', 'Hemlock', 'Cedar',
                                           'Poplar', 'Maple', 'Birch', 'Beech', 'Ironwood', 'Basswood', 'Ash', 'Elm', 'Oak'))
  vegclim.table
}