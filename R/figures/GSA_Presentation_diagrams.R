#  GSA Presentation Examples:

plot.densities <- function(subset.boolean){

  full.shift <- vegclim.table[subset.boolean,]
  
  full.shift$climate <- factor(full.shift$climate, levels = c('tmax', 'ppt'),
                             labels = c('T[max]', 'P[ann]'))
  
  full.shift$base <- factor(full.shift$base, levels = c('PLSS', 'FIA'))
  full.shift$c.ref <- factor(full.shift$c.ref, levels = c('PLSS', 'FIA'))
  
  subset.dplot <- ggplot(data = full.shift) +
    geom_density(aes(x = clim, y = ..scaled.., fill = base, 
                     linetype = c.ref, size = c.ref),
                 alpha = 0.2) +
    scale_size_manual(values = c(1, 1.5)) +
    scale_fill_manual(values = c('#F8766D', '#00BFC4')) +
    scale_linetype_manual(values = c(1, 2), drop = FALSE) +
    facet_grid(taxon~climate, scale = 'free_x', labeller = label_parsed) +
    theme_bw() +
    xlab('Climate Variable') +
    ylab('Kernel Density') +
    theme(axis.text = element_text(family = 'serif', size = 12),
          axis.title = element_text(family = 'serif', face = 'bold', size = 14),
          strip.text = element_text(family = 'serif', face='bold', size = 14),
          legend.position = 'none')
  
  subset.dplot
}

full.shift <- with(vegclim.table, taxon %in% c('Spruce', 'Hemlock', 'Elm') & 
                    data > 0 & climate  %in% c('ppt', 'tmax') &
                    ((base == 'PLSS' & c.ref == 'PLSS') | (base == 'FIA' & c.ref == 'FIA')))

plss.shift <- with(vegclim.table, taxon %in% c('Spruce', 'Hemlock', 'Elm') & 
                     data > 0 & climate %in% c('ppt', 'tmax') &
                     (base == 'PLSS'))

land.use   <- with(vegclim.table, taxon %in% c('Spruce', 'Hemlock', 'Elm') & 
                     data > 0 & climate %in% c('ppt', 'tmax') &
                     ((base == 'PLSS' & c.ref == 'FIA') | (base == 'FIA' & c.ref == 'FIA')))

three.levels <-  with(vegclim.table, taxon %in% c('Spruce', 'Hemlock', 'Elm') & 
                       data > 0 & climate %in% c('ppt', 'tmax') &
                       ((base == 'PLSS') | (base == 'FIA' & c.ref == 'FIA')))


ggsave(plot = plot.densities(full.shift), 
       filename='figures/full_shift.tiff', width = 8, height = 6, dpi = 150)

ggsave(plot = plot.densities(plss.shift), 
       filename='figures/plss_shift.tiff', width = 8, height = 6, dpi = 150)

ggsave(plot = plot.densities(land.use), 
       filename='figures/land_use.tiff', width = 8, height = 6, dpi = 150)

ggsave(plot = plot.densities(three.levels), 
       filename='figures/hellingerplot.tiff', width = 8, height = 6, dpi = 150)
