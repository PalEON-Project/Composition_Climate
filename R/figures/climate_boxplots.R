# plot resampling:

ppt.plot <- ggplot(subset(vegclim.table, data > 0.05 & climate == 'ppt'), 
                   aes(factor(paste0(base, ', ', c.ref)), clim)) + 
  geom_boxplot(alpha=0.2) +
  xlab('') + ylab('') +
  facet_grid(taxon ~ climate, scales='free') + theme_bw()
tma.plot <- ggplot(subset(vegclim.table, data > 0.05 & climate == 'tmax'), 
                   aes(factor(paste0(base, ', ', c.ref)), clim)) + 
  geom_boxplot(alpha=0.2) +
  xlab('') + ylab('') +            
  facet_grid(taxon ~ climate, scales='free') + theme_bw()
tmi.plot <- ggplot(subset(vegclim.table, data > 0.05 & climate == 'tmin'), 
                   aes(factor(paste0(base, ', ', c.ref)), clim)) + 
  geom_boxplot(alpha=0.2) +
  xlab('') + ylab('') +
  facet_grid(taxon ~ climate, scales='free') + theme_bw()
tme.plot <- ggplot(subset(vegclim.table, data > 0.05 & climate == 'tmean'), 
                   aes(factor(paste0(base, ', ', c.ref)), clim)) + 
  geom_boxplot(alpha=0.2) +
  xlab('') + ylab('') +
  facet_grid(taxon ~ climate, scales='free') + theme_bw()

grid.arrange(ppt.plot, tma.plot, tmi.plot, tme.plot, nrow=1, sub='')
