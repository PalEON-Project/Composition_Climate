#  To run the hellinger distances:

hellinger <- function(clim1, clim2){
  
  center <- mean(c(clim1, clim2), na.rm = TRUE)
  ranger <- diff(range(c(clim1, clim2),na.rm=TRUE))
  
  clim1 <- (clim1 - center) / ranger
  clim2 <- (clim2 - center) / ranger
  
  clim.range <- range(c(clim1, clim2), na.rm=TRUE)
  
  kernel1 <- scale(density(clim1, na.rm=TRUE, 
                     from = clim.range[1], to = clim.range[2])$y)
  kernel2 <- scale(density(clim2, na.rm=TRUE, 
                     from = clim.range[1], to = clim.range[2])$y)
  
  #  Calculate hellinger distance (between two kernels)
  1/sqrt(2) * sqrt(sum(sqrt(kernel1 - min(kernel1)) - sqrt(kernel2-min(kernel2)))^2)
  
}

all.hellinger <- ldply(1:length(best.taxa), function(x){
      ldply(1:4, function(y){
        data.subset <- vegclim.table[vegclim.table$taxon %in% best.taxa[x] & 
                                       vegclim.table$climate %in% bio.vals[y] & 
                                       vegclim.table$data > 0,]
        
        full_shift <- hellinger(subset(data.subset, base == 'PLSS' & c.ref == 'PLSS')$clim,
                               subset(data.subset, base == 'FIA' & c.ref == 'FIA')$clim)
        
        clim_shift <- hellinger(subset(data.subset, base == 'PLSS' & c.ref == 'PLSS')$clim,
                               subset(data.subset, base == 'PLSS' & c.ref == 'FIA')$clim)
        
        lu_shift   <- hellinger(subset(data.subset, base == 'FIA' & c.ref == 'FIA')$clim,
                               subset(data.subset, base == 'PLSS' & c.ref == 'FIA')$clim)
        
        data.frame(taxon = best.taxa[x],
                   climate = bio.vals[y],
                   full.change = full_shift,
                   fixedv_c = c(clim_shift, lu_shift),
                   type     = c("climate", "land use"))
      })
    })
