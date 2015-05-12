
agg.dens[,-1] <- agg.dens[,-1] / rowSums(agg.dens[,-1])
pls.data[,-1] <- pls.data[,-1]/rowSums(pls.data[,-1])

then.rast <- resample(then.rast, unit.raster)
now.rast <- resample(now.rast, unit.raster)

#For each taxon sample:
#  1.  The FIA data for one of the now.rast years for each of the climate variables.
#  2.  One of the then.rast years & one of the western layers
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
  values <- data.frame(cell  = c(rep(pls.data$cell, 2), rep(agg.dens$cell, 2)),
                       x     = c(rep(xyFromCell(unit.raster, pls.data$cell)[,'x'],2),
                                 rep(xyFromCell(unit.raster, agg.dens$cell)[,'x'],2)),
                       y     = c(rep(xyFromCell(unit.raster, pls.data$cell)[,'y'],2),
                                 rep(xyFromCell(unit.raster, agg.dens$cell)[,'y'],2)),
                       data  = c(rep(as.numeric(pls.data[,taxon]), 2),
                                 rep(agg.dens[,taxon],2)),
                       base  = c(rep('PLSS', nrow(pls.data)), 
                                 rep('PLSS', nrow(pls.data)), 
                                 rep('FIA', nrow(agg.dens)),
                                 rep('FIA', nrow(agg.dens))),
                       c.ref = c(rep('PLSS', nrow(pls.data)), 
                                 rep('FIA', nrow(pls.data)), 
                                 rep('FIA', nrow(agg.dens)),
                                 rep('PLSS', nrow(agg.dens))),
                       stringsAsFactors = FALSE)
  
#  The problem here is that both 
  
  plss.plss.clim <- getValues(then.rast)[subset(values, base %in% 'PLSS' & c.ref %in% 'PLSS')$cell, clim.var]
  plss.fia.clim  <- getValues(now.rast)[subset(values, base %in% 'PLSS' & c.ref %in% 'FIA')$cell, clim.var]
  
  fia.fia.clim  <- getValues(now.rast)[subset(values, base %in% 'FIA' & c.ref %in% 'FIA')$cell, clim.var]
  fia.plss.clim   <- getValues(then.rast)[subset(values, base %in% 'FIA' & c.ref %in% 'PLSS')$cell ,clim.var]
  
  values$clim <- c(plss.plss.clim, plss.fia.clim,
                   fia.fia.clim, fia.plss.clim)
  
  values
  
}

#  This gives us samples from the data for both PLSS and FIA: 
data.tables <-   llply(colnames(pls.data)[-1],
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
  aa$taxon <- colnames(pls.data)[-1][x]
  aa$climate <- bio.vals[y]
  aa}))

vegclim.table <- do.call(rbind.data.frame, lapply(data.out, function(x) do.call(rbind.data.frame, x)))
vegclim.table <- vegclim.table[!is.na(vegclim.table$data),]
