# Plotting the total climate-shift for the domain.
clim_boxplots <- function(clim_table = vegclim_table){
  
  clim_table$climate <- factor(clim_table$climate, 
                             levels = c('tmax', 'tmean', 'tmin', 'ppt'),
                             labels = c('T[max]', 'T[mean]', 'T[min]', 'P[ann]'))
  
  boxplot_clim <- ggplot(data = subset(na.omit(clim_table), 
                                       comb %in% c('FIA Era', 'PLS Era') & data>0 &
                                         !taxon %in% 'Ironwood')) +
    geom_boxplot(aes(y = clim, x = taxon, fill = comb)) +
    scale_fill_grey() +
    facet_grid(climate~., scales='free_y', labeller = label_parsed) +
    theme_bw() +
    theme(axis.text = element_text(family = 'serif', face='bold', size = 12),
          axis.title = element_blank(),
          strip.text = element_text(family = 'serif', face='bold', size = 14),
          legend.title = element_text(family = 'serif', face='bold', size = 14),
          legend.text = element_text(family = 'serif', size = 12),
          axis.text.x=element_text(angle=45, hjust = 1))
  
  ggsave(plot = boxplot_clim, filename = 'figures/boxplot_clim.tiff', width = 9, height = 6, dpi = 150)
  boxplot_clim
}