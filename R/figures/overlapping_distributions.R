
overlapping_clim <- ggplot(vegclim.table, aes(x = clim)) + 
  geom_density(aes(fill = c.ref), alpha = 0.5) + 
  facet_wrap(~climate, scale='free', nrow=1) +
  xlab('Climate Parameter') + ylab('Kernel Density') +
  theme(axis.text.x = element_text(family = 'serif', size = 12),
        axis.text.y = element_blank(),
        axis.title = element_text(family = 'serif', face = 'bold', size = 14),
        strip.text = element_text(family = 'serif', face='bold', size = 14),
        legend.title = element_text(family = 'serif', face='bold', size = 14)) +
  labs(fill = 'Era') +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

ggsave(overlapping_clim, 
       filename = 'figures/overlappingdensity.tiff', 
       dpi = 150, width = 8, height = 4)