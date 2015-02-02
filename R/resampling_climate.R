landuse <- list.files('../../../Dropbox/REU Work/Land Use/',pattern='.tif$', full.names=TRUE)

landuse <- landuse[c(1, 4, 6, 7, 8, 9, 10, 16)]
landuse.r <- (stack(landuse))
landuse.r[landuse.r>100] <- NA
landuse <- sum(landuse.r)>0

agg.dens[is.na(agg.dens)] <- 0
agg.dens[,-1] <- agg.dens[,-1] / rowSums(agg.dens[,-1])

all.taxa <- names(western$var)

all.taxa <- all.taxa[!all.taxa%in% c('Other hardwood', 'Atlantic White Cedar',
  'Buckeye', 'Ironwood', 'Walnut', 'Black gum/sweet gum')]

unit.raster <- setValues(base.rast, 1:ncell(base.rast))
proj4string(landuse) <- '+proj=longlat +ellps=WGS84'

landuse <- projectRaster(landuse, unit.raster) > 0.9

#For each taxon sample:
#  1.  The FIA data for one of the now.rast years for each of the climate variables.
#  2.  One of the then.rast years & one of the western layers
#  3.  Get vectors of values

hellinger <- function(kernel1, kernel2){
  #  Calculate hellinger distance (between two kernels)
  1/sqrt(2) * sqrt(sum(sqrt(kernel1) - sqrt(kernel2))^2)
}

ranges <- data.frame(min = c(100, 0, 220, 155, -383),
                     max = c(879, 403, 558, 405, 27))

# To do the analysis we want to get the climate layers for the past:
get_dens_era <- function(taxon, clim.var){
  
  #  Function resamples proportion data from the composition model for the region
  #  and then draws from the composition model and the climate data.
  #  The function outputs a table with six columns:
  #  data      - the proportion composition within the cell
  #  class     - the era of observation
  #  cell      - the literal cell # in the raster
  #  plss.clim - the climate extracted from the 1900 - 1930s era
  #  fia.clim  - the climate extracted from the 1970 - 2000s era
  #
  #  Each table has PLSS/pres (taxon present in cell at PLSS)
  #                 FIA/pres  (taxon present in cell at PLSS)
  #                 FIA/np    (taxon present in cell at PLSS but not in FIA)

  layer_ncdf <- sample(500, 1)                      # there are 500 layers in the ncdf
  taxon.vals <- ncvar_get(western, taxon, c(1,1,1), c(-1, -1, -1))

  #  bind all the 
  
  values <- data.frame(x     = c(rep(western.grid$x, 2),
                                 rep(xyFromCell(base.rast, agg.dens$cell)[,'x'],2),
                                 xyFromCell(landuse, agg.dens$cell)[,'x']),
                       y     = c(rep(western.grid$y, 2),
                                 rep(xyFromCell(base.rast, agg.dens$cell)[,'y'],2),
                                 xyFromCell(landuse, agg.dens$cell)[,'y']),
                       data  = c(rep(as.numeric(taxon.vals[,,layer_ncdf]), 2),
                                 rep(agg.dens[,taxon],2),
                                 landuse[agg.dens$cell]),
                       base  = c(rep('PLSS', length(western.grid$x)), 
                                 rep('PLSS', length(western.grid$x)), 
                                 rep('FIA', nrow(agg.dens)),
                                 rep('FIA', nrow(agg.dens)),
                                 rep('Land Use', nrow(agg.dens))),
                       c.ref = c(rep('PLSS', length(western.grid$x)), 
                                 rep('FIA', length(western.grid$x)), 
                                 rep('FIA', nrow(agg.dens)),
                                 rep('PLSS', nrow(agg.dens)),
                                 rep('FIA', nrow(agg.dens))),
                       stringsAsFactors = FALSE)
  
  coordinates(values) <- ~ x + y
  proj4string(values) <- CRS('+init=epsg:3175')
  
  values$cell <- extract(unit.raster, values)
  values$pres <- values$data > 0.1
  
  values <- values[values$cell %in% subset(values@data, base %in% 'PLSS')$cell,]
  
  mod_clim <- sample(nlayers(now.rast[[1]]), 1, replace = TRUE)
  pls_clim <- sample(nlayers(then.rast[[1]]), 1, replace = TRUE)
  
  plss.plss.clim <- do.call(cbind.data.frame, 
                            lapply(then.rast, function(x)extract(x[[pls_clim]], 
                                                                 values[values$base %in% 'PLSS' & values$c.ref %in% 'PLSS',])))[,clim.var]
  plss.fia.clim  <- do.call(cbind.data.frame, 
                            lapply(now.rast,  function(x)extract(x[[mod_clim]], values[values$base %in% 'PLSS' & values$c.ref %in% 'FIA',])))[,clim.var]

  fia.fia.clim  <- do.call(cbind.data.frame, lapply(now.rast,  function(x)extract(x[[mod_clim]], values[values$base %in% 'FIA' & values$c.ref %in% 'FIA',])))[,clim.var]
  fia.plss.clim <- do.call(cbind.data.frame, lapply(then.rast, function(x)extract(x[[pls_clim]], values[values$base %in% 'FIA' & values$c.ref %in% 'PLSS',])))[,clim.var]
  
  land.clim <- do.call(cbind.data.frame, lapply(now.rast, function(x)extract(x[[mod_clim]], values[values$base %in% 'Land Use',])))[,clim.var]
  
  values$clim <- c(plss.plss.clim, plss.fia.clim,
                   fia.fia.clim, fia.plss.clim,
                   land.clim)
  
  values@data
}

#  This gives us 100 samples from the data for both PLSS and FIA: 
data.tables <-   llply(all.taxa,
                        function(x){
                          lapply(1:5, function(y){
                            do.call(rbind.data.frame, 
                              lapply(1, function(z){
                                aa <- try(get_dens_era(x, y))
                                if(class(aa) == 'try-error'){
                                  aa <- data.frame(data=NA, base=NA, c.ref = NA, cell=NA, pres=NA, climate = NA)
                                }
                                aa}))
                          })
                        }, .progress = 'text')


data.out <- lapply(1:15, function(x)lapply(1:5, function(y){
  aa <- data.tables[[x]][[y]]
  aa$taxon <- all.taxa[x]
  aa$climate <- bio.vals[y]
  aa}))

aa <- lapply(data.out, function(x) do.call(rbind.data.frame, x))
bb <- do.call(rbind.data.frame, aa)

w.means <- ddply(bb, .(base, clim, taxon), summarise, 
      x = weighted.mean(climate, data, na.rm=TRUE))

ggplot(bb, aes(x = clim, fill = class)) + 
  geom_density(alpha = 0.2) + 
  facet_grid(taxon~climate, scale = 'free')

#  Now we have the presence or absence of taxa at PLSS and FIA eras
#  what are we looking for?
#  We want to be able to predict the shape of the PDF for a modern taxon distribution
#  along a climate gradient.  We need to generate the probability density functions:
#  1.  Generate the PDFs for the following:
#    a.  PLSS at PLSS climate (base climate envelope)
#    b.  PLSS at FIA  climate (implies spatial stability) - not used
#    c.  FIA  at PLSS climate (the climate space currently occupied) - not used
#    d.  FIA  at FIA  climate (modern envelope)
#    e.  FIA  at FIA absent.  (places in the FIA where there is no PLSS veg, Land Use)

get_dens <- function(x, class, pres, clim){
  x <- na.omit(x)

  if(class == 'FIA' & !pres)
  
  dens <- density(subset(x, class == class)[pres,clim], na.rm=TRUE,
                  from = range(x[,5:6])[1], to =  range(x[,5:6])[2])
  dens$y <- dens$y / sum(dens$y)
  dens
}

fia    <- get_dens(data.tables[[1]][[5]], class = 'FIA',  pres = pres,  clim = 'fia.clim')
plss   <- get_dens(data.tables[[1]][[5]], 'PLSS', pres,  'plss.clim')
plss.m <- get_dens(data.tables[[1]][[5]], 'PLSS', pres,  'fia.clim')
land   <- get_dens(data.tables[[1]][[5]], 'FIA',  !pres, 'fia.clim')

bin.data <- function(x){
  x <- na.omit(x)
  val.range <- c(seq(floor(range(x[,5:6])[1]), ceiling(range(x[,5:6])[2]), by = 10), 
                 ceiling(range(x[,5:6])[2]) + 10)
  
  plss.int <- hist(x=subset(x, class='PLSS')$plss.clim[x$pres], 
                   breaks = val.range, include.lowest = FALSE)$counts
  cc.int   <- hist(x=subset(x, class='PLSS')$fia.clim[x$pres], 
                   breaks = val.range, include.lowest = FALSE)$counts
  fia.int  <- hist(x=subset(x, class='FIA')$fia.clim[x$pres], 
                   breaks = val.range, include.lowest = FALSE)$counts
  lu.int   <- hist(x=subset(x, class='FIA')$fia.clim[!x$pres], 
                   breaks = val.range, include.lowest = FALSE)$counts
  
  pct <- function(x)x/sum(x)
  
  test.fun <- function(a, b){
    sum((pct(fia.int) - pct((pct(plss.int) + a * pct(fia.int) + b * pct(lu.int))))^2)
  }
  
  aa <- as.matrix(do.call(cbind.data.frame,
                          lapply(seq(-3, 3, by = 0.1), 
                                 function(x)do.call(rbind.data.frame,
                                                    lapply(seq(-3, 3, by = 0.1), test.fun, b = x)))))
  colnames(aa) <- 1:ncol(aa)
  aa[aa>1] <- 1
  image(aa)
}

plot((land$y - fia$y))

#  Hemlock test:
dens.diff <- density()


gen.comp <- function(x, y){
  matrix(sample(nrow(plss.tables[[x]][[y]]), 
                nrow(plss.tables[[x]][[y]]) * 2, 
                replace = TRUE),
         ncol = 2)
}

distances <- ldply(1:length(plss.tables),
      function(x){ 
        ldply(1:5, function(y, x){
          plss <- data.frame(taxon = all.taxa[x],
                             climate = y,
                             era     = 'plss',
                             distance = apply(gen.comp(x,y), 1, 
                                              function(z){
                                                hellinger(plss.tables[[1]][[1]][z[1],], 
                                                          plss.tables[[1]][[1]][z[2],])
                                              }))
          fia  <- data.frame(taxon = all.taxa[x],
                             climate = y,
                             era     = 'fia',
                             distance = apply(gen.comp(x, y), 1, 
                                              function(z){
                                                hellinger(fia.tables[[x]][[y]][z[1],], 
                                                          fia.tables[[x]][[y]][z[2],])
                                              }))
          
          inter <- data.frame(taxon = all.taxa[x],
                              climate = y,
                              era     = 'inter',
                              distance = apply(gen.comp(x, y), 1, 
                                               function(z){
                                                 hellinger(plss.tables[[x]][[y]][z[1],], 
                                                           fia.tables[[x]][[y]][z[2],])
                                               }))
          rbind.data.frame(plss, fia, inter)
      }, x = x)})

distances <- subset(distances, 
                    subset=era %in% c('plss', 'inter') & 
                      !taxon %in% c('Other hardwood', 'Atlantic White Cedar',
                                    'Buckeye', 'Ironwood', 'Walnut'))

summary(lm(distance ~ taxon*era*climate, data = distances))

ggplot(distances, aes(factor(c('Pwet', 'Pdry', 'Tdiff', 'Tmax', 'Tmin')[climate]), 
                      distance, fill = era)) + 
  geom_boxplot() + facet_wrap(~taxon) +
  xlab('Climate Variable') + ylab('Hellinger Distance')

optima.plss <- lapply(plss.tables, function(x){lapply(1:5, function(i){
  xaxis <- seq(ranges[i,1], ranges[i,2], length.out = 100)
  weights <- apply(x[[i]], 1, function(y) sum(xaxis*y)/sum(y))
  weights
})})

optima.fia <- lapply(fia.tables, function(x){lapply(1:5, function(i){
  xaxis <- seq(ranges[i,1], ranges[i,2], length.out = 100)
  weights <- apply(x[[i]], 1, function(y) sum(xaxis*y)/sum(y))
  weights
})})

optim.set <- data.frame(taxon = rep(rep(as.character(unique(distances$taxon)), each = 500),2),
           climate = rep(rep(c('Pwet', 'Pdry', 'Tdiff', 'Tmax', 'Tmin'), each = 100), 30),
           era = rep(c('PLSS', 'FIA'), each = 7500),
           optima = c(unlist(optima.plss), unlist(optima.fia)),
           stringsAsFactors=FALSE)

optim.set$taxon[optim.set$taxon %in% 'Cedar/juniper'] <- 'Cedar'
optim.set$taxon[optim.set$taxon %in% "Poplar/tulip poplar"] <- 'Poplar'

optim.set$era <- factor(optim.set$era, levels = c('PLSS', 'FIA'))
optim.set$taxon <- factor(optim.set$taxon, 
                          levels = c('Tamarack', 'Spruce', 'Pine', 'Fir', 'Cedar', 'Hemlock',
                                     'Poplar', 'Birch', 'Ash', 'Basswood', 'Beech', 'Elm', 
                                     'Hickory', 'Maple', 'Oak'))

#  Big Silly Boxplot
ggplot(optim.set, aes(taxon, optima, fill=era)) + geom_boxplot() + 
  facet_wrap(~climate, ncol = 1, scale='free_y') +
  xlab('Taxon') + ylab('Climatic Variable') +
  scale_fill_discrete(breaks=c("PLSS","FIA")) +
  theme(axis.title.x = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.title.y = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        axis.text.y = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.text = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.title = element_text(family = 'serif', 
                                    face = 'bold', size = 18, color = 'black'),
        strip.text =  element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'))
        

taxon <- unique(optim.set$taxon)

#  Summarize Differences:
optim.diff <- data.frame(diff = unlist(optima.fia) - unlist(optima.plss),
                         taxon = rep(as.character(unique(distances$taxon)), each = 500),
                         climate = rep(rep(c('Pwet', 'Pdry', 'Tdiff', 'Tmax', 'Tmin'), each = 100), 15))

optim.diff <- dcast(optim.diff, taxon ~ climate, value.var='diff', fun.aggregate=mean)

rownames(optim.diff) <- optim.diff[,1]; optim.diff <- optim.diff[,-1]

rownames(optim.diff)[rownames(optim.diff) %in% 'Cedar/juniper'] <- 'Cedar'
rownames(optim.diff)[rownames(optim.diff) %in% "Poplar/tulip poplar"] <- 'Poplar'


library(MASS)
library(vegan)
MDS.plot <- metaMDS(apply(optim.diff, 2, scale)+3)

mds.gg <- data.frame(x = MDS.plot$points[,1],
                     y = MDS.plot$points[,2],
                     label = rownames(optim.diff))
clim.gg <- data.frame(x = MDS.plot$species[,1],
                     y = MDS.plot$species[,2],
                     label = colnames(optim.diff))

ggplot(mds.gg, aes(x = x, y = y, label = label)) + 
  geom_vline(xintercept = 0, alpha = 0.3) + 
  geom_hline(yintercept = 0, alpha = 0.3) + 
  geom_text(family = 'serif') +
  geom_text(data = clim.gg, aes(x = x, y = y, label = label), color = 'red', 
            family = 'serif', face = 'bold') +
  xlab('MDS Axis 1') + ylab('MDS Axis 2') +
  theme_bw() +
  theme(axis.title.x = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.title.y = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        axis.text.y = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'))

##  Plot out the curves:
optima.curves.plss <- do.call(rbind.data.frame, lapply(1:length(taxon), function(x){
  do.call(rbind.data.frame,lapply(1:5, function(y){
    output <- data.frame(taxon = taxon[x],
                         kernel = colMeans(plss.tables[[x]][[y]]),
                         climate = c('Pwet', 'Pdry', 'Tdiff', 'Tmax', 'Tmin')[y],
                         values = seq(ranges[y,1], ranges[y,2], length.out = 100),
                         stringsAsFactors = FALSE)
    output
    }))
  }))
optima.curves.fia <- do.call(rbind.data.frame, lapply(1:length(taxon), function(x){
  do.call(rbind.data.frame,lapply(1:5, function(y){
    output <- data.frame(taxon = taxon[x],
                         kernel = colMeans(fia.tables[[x]][[y]]),
                         climate = c('Pwet', 'Pdry', 'Tdiff', 'Tmax', 'Tmin')[y],
                         values = seq(ranges[y,1], ranges[y,2], length.out = 100),
                         stringsAsFactors = FALSE)
    output
  }))
}))

optima.curves <- data.frame(rbind(optima.curves.plss, optima.curves.fia),
                            era = factor(rep(c('PLSS', 'FIA'), each = nrow(optima.curves.plss))))

ggplot(optima.curves[optima.curves$climate == 'Tmax' & optima.curves$taxon %in% c('Oak', 'Basswood', 'Beech'),], 
       aes(x = values/10, y = kernel, color = era)) +
  geom_line(size = 2) + facet_wrap(~taxon, scales = 'free_x') +
  theme_bw() +
  xlab('Mean Maximum Temperature') +
  ylab('Kernel Density') +
  theme(axis.title.x = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.title.y = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        axis.text.y = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.text = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.title = element_text(family = 'serif', 
                                    face = 'bold', size = 18, color = 'black'),
        strip.text =  element_text(family = 'serif', 
                                   face = 'bold', size = 24, color = 'black'))


ggplot(optima.curves[optima.curves$climate == 'Pwet' & optima.curves$taxon %in% c('Oak', 'Basswood', 'Beech'),], 
       aes(x = values/10, y = kernel, color = era)) +
  geom_line(size = 2) + facet_wrap(~taxon, scales = 'free_x') +
  theme_bw() +
  xlab('Precipitation (mm) of Dry Quarter') +
  ylab('Kernel Density') +
  scale_x_log10() +
  theme(axis.title.x = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.title.y = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        axis.text.y = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.text = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.title = element_text(family = 'serif', 
                                    face = 'bold', size = 18, color = 'black'),
        strip.text =  element_text(family = 'serif', 
                                   face = 'bold', size = 24, color = 'black'))


ggplot(optima.curves[optima.curves$taxon %in% c('Hemlock'),], 
       aes(x = values, y = kernel, color = era)) +
  geom_line(size = 2) + facet_wrap(~climate, scales = 'free_x') +
  theme_bw() +
  xlab('Climate Variables') +
  ylab('Kernel Density') +
  theme(axis.title.x = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.title.y = element_text(family = 'serif', 
                                    face = 'bold.italic', size = 24),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        axis.text.y = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.text = element_text(family = 'serif', 
                                   face = 'bold', size = 18, color = 'black'),
        legend.title = element_text(family = 'serif', 
                                    face = 'bold', size = 18, color = 'black'),
        strip.text =  element_text(family = 'serif', 
                                   face = 'bold', size = 24, color = 'black'))

