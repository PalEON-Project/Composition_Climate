# Plotting the total climate-shift for the domain.
clim_boxplots <- function(clim_table = vegclim_table){
  
  clim_table$taxon <- mapvalues(clim_table$taxon, from = levels(clim_table$taxon),
                                to = c("Larix", "Pinus", "Picea", "Abies", "Tsuga",
                                       "Thuja/Juniperus", "Populus", "Acer", "Betula", "Fagus", 
                                       "Ostrya/Carpinus", "Tilia", "Fraxinus", "Ulmus", "Quercus"))
  
  clim_table$climate <- factor(clim_table$climate, 
                             levels = c('tmax', 'tmean', 'tmin', 'ppt'),
                             labels = c('T[max]', 'T[mean]', 'T[min]', 'P[ann]'))
  
  boxplot_clim <- ggplot(data = subset(na.omit(clim_table), 
                                       comb %in% c('FIA Era', 'PLS Era') & data > 0 &
                                         !taxon %in% 'Ironwood')) +
    geom_boxplot(aes(y = clim, x = taxon, fill = comb), lwd = 0.2, width = 1, outlier.size = 0.5) +
    scale_fill_grey(name = 'Climate Era') +
    facet_grid(climate~., scales = 'free_y', labeller = label_parsed) +
    theme_bw() +
    theme(axis.text = element_text(family = 'serif', face = 'bold', size = 12),
          axis.title = element_blank(),
          strip.text = element_text(family = 'serif', face = 'bold', size = 14),
          legend.title = element_text(family = 'serif', face = 'bold', size = 14),
          legend.text = element_text(family = 'serif', size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold.italic'),
          plot.margin = grid::unit(c(10,10,10,10), "mm"))
  
  ggsave(plot = boxplot_clim, filename = 'figures/boxplot_clim.pdf', width = 6.7, height = 4.5)
  ggsave(plot = boxplot_clim, filename = 'Final_Figures/boxplot_clim.pdf', width = 6.7, height = 4.5)
  boxplot_clim
}