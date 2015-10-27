#  Plotting for taxa that shows the change between historical and modern distributions.

get_loss_gain <- function(unit.raster, best.taxa, pls_data, agg_dens){
  loss.gain <- data.frame(cell = 1:ncell(unit.raster))
  
  agg_dens[is.na(agg_dens)] <- 0
  
  for(i in 1:length(best.taxa)){
    loss.gain[,(i+1)] <- NA
    colnames(loss.gain)[i+1] <- best.taxa[i]
    
    sub.pls <- na.omit(pls_data[pls_data[,best.taxa[i]] > 0 & !is.na(pls_data$cell),])
    sub.fia <- na.omit(agg_dens[agg_dens[,best.taxa[i]] > 0 & !is.na(agg_dens$cell),])
    
    #  PLS is there, but FIA is not.
    loss.gain[sub.pls$cell[!sub.pls$cell %in% sub.fia$cell],best.taxa[i]] <- 'Loss'
    #  PLS and FIA are there:
    loss.gain[sub.pls$cell[sub.pls$cell %in% sub.fia$cell],best.taxa[i]] <- 'Presence'
    # FIA is there, PLS is not:
    loss.gain[sub.fia$cell[!sub.fia$cell %in% sub.pls$cell],best.taxa[i]] <- 'Gain'
    
  }
  
  loss.melt <- na.omit(melt(loss.gain, id = 'cell'))
  
  loss.melt <- data.frame(xyFromCell(unit.raster, loss.melt$cell),
                          loss.melt)
  loss.melt
}

loss_plot <- function(){
  
  if(model.proj == '+init=epsg:4326'){
    ext <- c(-98, -83, 42, 50)
  }
  
  if(model.proj == '+init=epsg:3175'){
    ext <- c(-100000, 1050000, 600000, 1600000)
  }
  
  loss_melt <- get_loss_gain(unit.raster, best.taxa, pls_data, agg_dens)

  loss_melt$rgb <- brewer.pal(n = 3, 'Dark2')[as.numeric(factor(loss_melt$value, 
                                                                levels = c('Presence', 'Loss', 'Gain')))]
  
  loss_melt$rgb[loss_melt$value == 'Gain'] <- "#FF0000"
  loss_melt$rgb[loss_melt$value == 'Presence'] <- colorspace::hex(colorspace::RGB(0.1, 0.1, 0.1))
  loss_melt$rgb[loss_melt$value == 'Loss'] <- colorspace::hex(colorspace::RGB(2/255,119/255,160/255))
  
  loss_melt$variable <- factor(loss_melt$variable,
                               levels = c('tamarack', 'pine', 'spruce', 'fir', 'hemlock', 'cedar.juniper',
                                          'poplar', 'maple', 'birch', 'beech', 'ironwood', 'basswood', 'ash', 'elm', 'oak'),
                               labels = c('Tamarack', 'Pine', 'Spruce', 'Fir', 'Hemlock', 'Cedar',
                                          'Poplar', 'Maple', 'Birch', 'Beech', 'Ironwood', 'Basswood', 'Ash', 'Elm', 'Oak'))
  
  loss_plot <- natural$base_map +
      geom_tile(data = natural$rast_table, aes(x = x, y = y, fill = rgb)) +
      geom_tile(data = subset(loss_melt, variable %in% 'Pine') , aes(x = x, y = y, fill = rgb)) +
      #facet_wrap(~variable, ncol = 5) + theme_bw() +
      scale_fill_identity() +
  #    geom_path(data = umw.domain, aes(x = long, y = lat, group = group)) +
  #    geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            strip.text = element_text(family = 'serif', face='bold', size = 14),
            legend.title = element_text(family = 'serif', face='bold', size = 14),
            legend.text = element_text(family = 'serif', size = 12)) +
    coord_equal(xlim = ext[1:2], ylim=ext[3:4])
  
  ggsave(plot = loss_plot, filename = 'figures/loss_plots.tiff', dpi = 150, width = 6, height = 4)
  
  short_plot <- natural$base_map +
    geom_tile(data = natural$rast_table, aes(x = x, y = y, fill = rgb)) +
    geom_tile(data = subset(loss_melt, variable %in% 'Oak') , aes(x = x, y = y, fill = rgb)) +
    scale_fill_identity() +
    facet_wrap(~variable, ncol = 5) + theme_bw() +
    #geom_path(data = umw.domain, aes(x = long, y = lat, group = group), color = 'black') +
    #geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(family = 'serif', face='bold', size = 14),
          legend.title = element_text(family = 'serif', face='bold', size = 14),
          legend.text = element_text(family = 'serif', size = 12)) +
    coord_equal(xlim = ext[1:2], ylim=ext[3:4])
  
  ggsave(plot = short_plot, filename = 'figures/short_plots.tiff', dpi = 150, width = 6, height = 5)
  
  
  list(loss_melt, loss_plot)
}

