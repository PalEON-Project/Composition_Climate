## ----setup-knitr, echo = FALSE, message = FALSE--------------------------
knitr::opts_chunk$set(
  comment = " ",
  error = TRUE,
  eval = TRUE,
  cache = FALSE,
  tidy = TRUE

)

library("pander")
library("plyr")
panderOptions('table.style', 'rmarkdown')

## ---- echo = FALSE, warning = FALSE, message=FALSE, results='hide'-------
# Load climate data and generate the pre-sett and modern bounds:
#  Pre-settlement climate: 1900 - 1930
#  Modern climate: 1960 - 1990

source('R/Setup_R.R')
source('R/load_climate.R') # gives `climate` with elements `then`, `now` and `names`
source('R/load_veg.R')

best.taxa <- c('beech', 'ironwood', 'hemlock', 'fir', 'cedar.juniper', 
              'basswood', 'spruce', 'elm', 'ash', 'pine', 'oak', 
              'tamarack', 'maple', 'poplar', 'birch')

pls_data <- get_pls('0.9-4')
agg_dens <- get_fia()[,colnames(pls_data)]

land_use <- land_use[land_use$cell %in% c(pls_data$cell, agg_dens$cell),]

source('R/natural_earth.R')

natural <- natural_earth()
  

## ----figure1, echo = TRUE, warning = FALSE, message=FALSE----------------

land_use$class[land_use$class == 'Water'] <- NA

land_use_plot <- natural[[2]] + 
  geom_tile(data = land_use, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(values = c('#ff7f00', '#e41a1c', '#33a02c'))


source('R/figures/climate_ts.R')

climate_values <- climate_ts(unit.raster, pls_data, climate$now, climate$then)

clim_change <- climate_values$modern$value - climate_values$past$value


## ----figureTime, fig.width=6, fig.height=5-------------------------------
grid.arrange(climate_values$plot, land_use_plot, nrow=1, widths = c(3, 4))

clim_change[3]

## ----loss_gain_fig, warning=FALSE, message=FALSE, echo=FALSE-------------

source('R/figures/gain_loss_plot.R')

#  This takes a bit of time.
loss_melt <- loss_plot()

loss_table <- as.data.frame.matrix(table(loss_melt[[1]]$variable, loss_melt[[1]]$value))
loss_table$Gain <- loss_table$Gain / (loss_table$Loss + loss_table$Presence)
loss_table[,2:3] <- tran(loss_table[,2:3], method = 'proportion')

loss_melt[[2]]

## ---- results='as.is'----------------------------------------------------
pandoc.table(loss_table, justify = "left", style = 'simple')

## ----densities, warning=FALSE, message=FALSE, echo=FALSE-----------------

source('R/resampling_climate.R')
vegclim_table <- resample_climate(climate, agg_dens, pls_data)

source('R/figures/climate_boxplots.R')
boxplot_clim <- clim_boxplots(clim_table = vegclim_table)
boxplot_clim


## ----climate_shifts, warning=FALSE, message=FALSE, echo=FALSE------------

source('R/hellinger_dists.R')

#  Returns a list with the ggplot (object 1) and all.hellinger.
hell_plot <- hellinger_plot()

hell_plot[[1]]


