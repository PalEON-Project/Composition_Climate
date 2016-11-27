#  To run the hellinger distances:

hellinger <- function(clim1, clim2){
  
  center <- mean(c(clim1, clim2), na.rm = TRUE)
  ranger <- diff(range(c(clim1, clim2),na.rm = TRUE))
  
  clim1 <- (clim1 - center) / ranger
  clim2 <- (clim2 - center) / ranger
  
  clim.range <- range(c(clim1, clim2), na.rm = TRUE)
  
  kernel1 <- scale(density(clim1, na.rm = TRUE, 
                     from = clim.range[1], 
                     to = clim.range[2], bw = "SJ")$y)
  kernel2 <- scale(density(clim2, na.rm = TRUE, 
                     from = clim.range[1], 
                     to = clim.range[2], bw = "SJ")$y)
  
  #  Calculate hellinger distance (between two kernels)
  1/sqrt(2) * sqrt(sum(sqrt(kernel1 - min(kernel1)) - sqrt(kernel2 - min(kernel2))) ^ 2)
  
}

hellinger_plot <- function(){
  
  vegclim_table$taxon <- mapvalues(vegclim_table$taxon, from = levels(vegclim_table$taxon),
                                to = c("Larix", "Pinus", "Picea", "Abies", "Tsuga",
                                       "Thuja/Juniperus", "Populus", "Acer", "Betula", "Fagus", 
                                       "Ostrya/Carpinus", "Tilia", "Fraxinus", "Ulmus", "Quercus"))
  
  best.taxa <- as.character(unique(vegclim_table$taxon))
  
  all.hellinger <- ldply(1:length(best.taxa), function(x){
        ldply(1:4, function(y){
          
          tree     <- best.taxa[x]
          clim_var <- climate$name[y]
          
          data.subset <- vegclim_table %>% filter(taxon == tree & climate %in% clim_var & data > 0)
          
          VHCH <- data.subset %>% filter(base == 'PLSS' & c.ref == 'PLSS') %>% select(clim) %>% unlist
          VHCM <- data.subset %>% filter(base == 'PLSS' & c.ref == 'FIA') %>% select(clim) %>% unlist
          VMCM <- data.subset %>% filter(base == 'FIA' & c.ref == 'FIA') %>% select(clim) %>% unlist
          
          full_shift <- hellinger(VHCH, VMCM)
          
          clim_shift <- hellinger(VHCH, VHCM)
          
          lu_shift   <- hellinger(VMCM, VHCM)
          
          data.frame(taxon = best.taxa[x],
                     climate = climate$name[y],
                     full.change = full_shift,
                     fixedv_c = c(clim_shift - lu_shift))
        })
      })
  
  levels(all.hellinger$climate) <- c('Annual Precipitation', 
                                     'Maximum Temperature',
                                     'Mean Temperature',
                                     'Minimum Temperature')

  hell_plot <- ggplot(data = subset(all.hellinger, !taxon %in% 'Ostrya/Carpinus'), 
                      aes(x = full.change, y = fixedv_c)) + 
    geom_text(aes(label = taxon), 
              family = 'serif', size = 4, fontface = 'italic') +
    facet_wrap(~climate) +
    theme_bw() +
    scale_x_log10(breaks = c(2, 10, 25, 50, 125), minor_breaks = NULL) +
    scale_y_continuous(minor_breaks = NULL) +
    xlab(expression("Total Change"*" - "*PLS*" - "*FIA*" - "*d[tot])) +
    ylab(expression('Climate Change - Land Use Change - '*d[c]*' - '*d[v])) +
    theme(axis.text = element_text(family = 'serif', size = 12),
          axis.title = element_text(family = 'serif', face = 'bold', size = 12),
          strip.text = element_text(family = 'serif', face = 'bold', size = 12),
          legend.position = "none",
          plot.margin = grid::unit(c(10,10,10,10), "mm")) +
    coord_cartesian(xlim = c(1,125), ylim = c(-65, 85)) +
    geom_abline(intercept = 0, slope = 0, alpha = 0.5) +
    annotate(geom = 'text', x = 7, y = 70, 
             label = 'Climate Dominates', family = 'serif',
             fontface = 'bold', size = 3) +
    annotate(geom = 'text', x = 7, y = -40, 
             label = 'Land Use Dominates', family = 'serif',
             fontface = 'bold', size = 3)
  
  ggsave(plot = hell_plot, filename = 'figures/hellingerplot.png', 
         width = 6.7, height = 6, dpi = 300)
  ggsave(plot = hell_plot, filename = 'Final_Figures/hellingerplot.pdf', 
         width = 6.7, height = 6, dpi = 300)
  
  list(hell_plot, all.hellinger)
}