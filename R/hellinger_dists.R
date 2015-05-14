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
                   fixedv_c = c(clim_shift - lu_shift))
      })
    })

levels(all.hellinger$climate) <- c('Annual Precipitation', 
                                   'Maximum Temperature',
                                   'Mean Temperature',
                                   'Minimum Temperature')

levels(all.hellinger$taxon) <- c('Beech', 'Ironwood', 'Hemlock', 'Fir', 'Cedar',
                                 'Basswood', 'Spruce', 'Elm', 'Ash', 'Pine',
                                 'Oak', 'Tamarack', 'Maple', 'Poplar', 'Birch')

hell.plot <- ggplot(all.hellinger, aes(x = full.change, y = fixedv_c)) + 
  geom_text(aes(label = taxon), family = 'serif', size = 4, fontface = 'italic') +
  facet_wrap(~climate) +
  theme_bw() +
  scale_x_sqrt() +
  #scale_y_sqrt() +
  xlab('Total Change - PLS - FIA') +
  ylab('Climate Change - Land Use Change') +
  theme(axis.text = element_text(family = 'serif', size = 12),
        axis.title = element_text(family = 'serif', face = 'bold', size = 14),
        strip.text = element_text(family = 'serif', face='bold', size = 14),
        legend.position = "none") +
  coord_cartesian(xlim=c(0,150), ylim=c(-85, 85)) +
  geom_abline(intercept=0, slope=0) +
  annotate(geom = 'text', x = 7, y = 70, 
           label='Climate Dominates', family='serif',
           fontface = 'bold', size = 5) +
  annotate(geom = 'text', x = 7, y = -40, 
           label='Land Use Dominates', family='serif',
           fontface = 'bold', size = 5)

ggsave(plot = hell.plot, filename = 'figures/hellingerplot.tiff', 
       width = 8, height = 8, dpi = 150)
