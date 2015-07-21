#  Trying to generate the MDS ordination for the climate data.
library(vegan)

hellinger.mat <- function(x){
  x$prop <- x$fixedv_c / x$full.change
  new.climate <-  dcast(data = x, 
                        formula = taxon ~ climate, 
                        fun.aggregate = sum, 
                        value.var = 'prop')
  
  rownames(new.climate) <- new.climate[,1]
  new.climate[,-1]
}

multi.var.hellinger <- hellinger.mat(subset(all.hellinger, !taxon %in% 'Ironwood'))

#  We drop beech because the change in precip is so big it blows everything up.
multi.sqrt <- apply(multi.var.hellinger[-1,], 2, function(x) sqrt(x - min(x)+0.1))

ordination <- metaMDS(multi.sqrt, k = 2)

ordi_plot <- data.frame(taxon = rownames(ordination$points),
                        ordination$points)
ordi_spec <- data.frame(clim = rownames(ordination$species),
                        ordination$species)

ordi.plot <- ggplot(ordi_plot, aes(x = MDS1, y = MDS2, label = taxon)) +
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 2) +
  geom_vline(xintercept=0, alpha = 0.5, linetype = 2) +
  geom_text(family = 'serif') +
  geom_text(data = ordi_spec, aes(x = MDS1, y = MDS2, label = clim), fontface = 'italic') +
  theme_bw() +
  theme(axis.text = element_text(family = 'serif', size = 12),
    axis.title = element_text(family = 'serif', face = 'bold', size = 14),
    strip.text = element_text(family = 'serif', face='bold', size = 14),
    legend.position = "none") +
  coord_cartesian(xlim=c(-0.55,.55), ylim=c(-0.35, .35))

ggsave(plot = ordi.plot, filename = 'figures/ordi_plot.tiff', dpi = 150, width = 8, height = 8)
