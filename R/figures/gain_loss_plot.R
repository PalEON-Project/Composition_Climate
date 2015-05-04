loss.gain <- data.frame(cell = 1:ncell(unit.raster))

for(i in 1:length(best.taxa)){
  loss.gain[,(i+1)] <- NA
  colnames(loss.gain)[i+1] <- best.taxa[i]
  
  sub.pls <- na.omit(pls.data[pls.data[,best.taxa[i]] > 0 & !is.na(pls.data$cell),])
  sub.fia <- na.omit(agg.dens[agg.dens[,best.taxa[i]] > 0 & !is.na(agg.dens$cell),])
  
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

ggplot(loss.melt, aes(x = x, y = y, fill = factor(value))) + geom_tile() +
  facet_wrap(~variable) + theme_bw() +
  scale_fill_brewer(type='qual', name = 'Change State') +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(family = 'serif', face='bold', size = 14),
        legend.title = element_text(family = 'serif', face='bold', size = 14),
        legend.text = element_text(family = 'serif', size = 12))

