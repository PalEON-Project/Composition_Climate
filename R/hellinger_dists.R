#  To run the hellinger distances:

hellinger <- function(clim1, clim2){
  
  center <- mean(c(clim1, clim2), na.rm = TRUE)
  ranger <- diff(range(c(clim1, clim2),na.rm=TRUE))
  
  clim1 <- (clim1 - center) / ranger
  clim2 <- (clim2 - center) / ranger
  
  clim.range <- range(c(clim1, clim2), na.rm=TRUE)
  
  kernel1 <- density(clim1, na.rm=TRUE, 
                     from = clim.range[1], to = clim.range[2])$y
  kernel2 <- density(clim2, na.rm=TRUE, 
                     from = clim.range[1], to = clim.range[2])$y
  
  #  Calculate hellinger distance (between two kernels)
  1/sqrt(2) * sqrt(sum(sqrt(kernel1) - sqrt(kernel2))^2)
  
}

all.hellinger <- ldply(1:length(best.taxa), function(x){
      ldply(1:4, function(y){
        data.subset <- vegclim.table[vegclim.table$taxon %in% best.taxa[x] & 
                                       vegclim.table$climate %in% bio.vals[y] & 
                                       vegclim.table$data > 0,]
        
        p_shift_c <- hellinger(subset(data.subset, base == 'PLSS' & c.ref == 'PLSS')$clim,
                               subset(data.subset, base == 'PLSS' & c.ref == 'FIA')$clim)
        f_shift_c <- hellinger(subset(data.subset, base == 'FIA' & c.ref == 'PLSS')$clim,
                               subset(data.subset, base == 'FIA' & c.ref == 'FIA')$clim)
        p_shift_v <- hellinger(subset(data.subset, base == 'PLSS' & c.ref == 'PLSS')$clim,
                               subset(data.subset, base == 'FIA' & c.ref == 'PLSS')$clim)
        f_shift_v <- hellinger(subset(data.subset, base == 'FIA' & c.ref == 'FIA')$clim,
                               subset(data.subset, base == 'PLSS' & c.ref == 'FIA')$clim)
        
        data.frame(taxon = best.taxa[x],
                   climate = bio.vals[y],
                   fixedv_c = c(p_shift_c,f_shift_c),
                   shift=c('PLS', 'FIA'),
                   fixedc_v = c(p_shift_v, f_shift_v))
      })
    })
