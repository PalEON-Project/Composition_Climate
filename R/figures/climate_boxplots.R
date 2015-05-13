# Plotting the total climate-shift for the domain.

boxplot_clim <- ggplot(data = subset(na.omit(vegclim.table), 
                                     comb %in% c('FIA Era', 'PLS Era') & data>0)) +
  geom_boxplot(aes(y = clim, x = taxon, fill = comb)) +
  scale_fill_grey() +
  facet_wrap(~climate, ncol = 1, scales='free_y') +
  theme_bw() +
  theme(axis.text = element_text(family = 'serif', face='bold', size = 12),
        axis.title = element_blank(),
        strip.text = element_text(family = 'serif', face='bold', size = 14),
        legend.title = element_text(family = 'serif', face='bold', size = 14),
        legend.text = element_text(family = 'serif', size = 12))

ggsave(plot = boxplot_clim, filename = 'figure/boxplot_clim.tiff', width = 8, height = 6, dpi = 150)