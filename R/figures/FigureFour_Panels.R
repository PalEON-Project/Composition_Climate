#  This is a hellinger plot of the four curves for Oak in vegclim_table:

library(png)
library(grid)

ggplotColours <- function(n=6, h=c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

model_grid <- readPNG('figures/model_grid.png')

  fourpanel <- function(clim_table, taxon) {
    veg <- clim_table %>% filter(taxon == taxon & data > 0 & climate == "ppt")
    
    veg$base <- factor(veg$base, levels = c("PLSS", "FIA"))
    
    veg$c.ref = factor(veg$c.ref, levels = c("PLSS", "FIA"))
    
    curves <- ggplot(veg) +
      geom_density(adjust = 1, size = 1.1,
                             aes(x = clim, color = base, linetype = c.ref)) +
      scale_color_manual(name = 'Vegetation', values = (ggplotColours(2))) +
      scale_linetype_manual(name = 'Climate', values = c(1,3)) +
      facet_grid(c.ref ~ base) +
      theme_bw() +
      scale_x_continuous(breaks = c(500, 700, 900)) +
      theme(axis.text = element_text(family = 'serif', face = 'bold', size = 12),
            axis.title = element_blank(),
            strip.text = element_text(family = 'serif', face = 'bold', size = 14),
            legend.title = element_text(family = 'serif', face = 'bold', size = 14),
            legend.text = element_text(family = 'serif', size = 12),
            legend.position = 'none')
    
    ggsave(plot = grid.arrange(rasterGrob(model_grid), curves, ncol = 2,
                             widths = c(3,4)), 
           file = 'figures/figureFour_panels.tiff',
           width = 6, height = 3, dpi = 300)
    
    curves
  }
  
  aa <- fourpanel(vegclim_table, "Oak")
