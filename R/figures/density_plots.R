#  Forget the backwards plot.
density.plots <- ggplot(subset(na.omit(vegclim.table), data > 0 & 
                                 !(taxon %in% 'ironwood') & 
                                 (base == 'PLSS' | (base == 'FIA' & c.ref == 'FIA')))) + 
    geom_density(aes(x = clim, y = ..scaled.., size = base, fill = c.ref), alpha = 0.2) +
  scale_size_manual(values = c(0.5, 1)) +
    facet_grid(taxon ~ climate, scales = 'free') +
  theme_bw()  +
    theme(axis.text = element_text(family = 'serif', size = 12),
          axis.title = element_text(family = 'serif', face = 'bold', size = 14),
          strip.text = element_text(family = 'serif', face='bold', size = 14))

ggsave(density.plots, filename = 'figures/densityplots.tiff', height = 8, width = 8, dpi = 150)
