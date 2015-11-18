#  This is a hellinger plot of the four curves for Oak in vegclim_table:

library(png)
library(grid)

ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

model_grid <- readPNG('figures/model_grid.png')

oak <- subset(vegclim_table, taxon == 'Tamarack' & data > 0)

oak_curves <- ggplot(oak) +
  geom_density(adjust=1, size = 1.2,
                         aes(x = clim, color = base, linetype = c.ref)) +
            facet_wrap(~climate, scale = 'free') +
            theme_bw() +
            scale_color_manual(name = 'Vegetation', values = rev(ggplotColours(2))) +
            scale_linetype_manual(name = 'Climate', values = c(2,1)) +
            theme(axis.text = element_text(family = 'serif', face='bold', size = 12),
                  axis.title = element_blank(),
                  strip.text = element_text(family = 'serif', face='bold', size = 14),
                  legend.title = element_text(family = 'serif', face='bold', size = 14),
                  legend.text = element_text(family = 'serif', size = 12),
                  legend.position = 'none')

ggsave(plot=grid.arrange(rasterGrob(model_grid), oak_curves, ncol=2,
                         widths = c(3,4)), 
       file = 'figures/figureFour_panels.tiff',
       width=6, height = 3, dpi = 300)
