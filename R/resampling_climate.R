all.taxa <- names(western$var)

all.taxa <- all.taxa[!all.taxa%in% c('Other hardwood', 'Atlantic White Cedar',
  'Buckeye', 'Ironwood', 'Walnut', 'Black gum/sweet gum')]

#For each taxon sample:
#  1.  The FIA data for one of the now.rast years for each of the climate variables.
#  2.  One of the then.rast years & one of the western layers
#  3.  Get vectors of values

hellinger <- function(kernel1, kernel2){
  1/sqrt(2) * sqrt(sum(sqrt(kernel1) - sqrt(kernel2))^2)
}

ranges <- data.frame(min = c(100, 0, 220, 155, -383),
                     max = c(879, 403, 558, 405, 27))

get_dens_plss <- function(taxon, clim.var){
  layer_ncdf <- sample(500, 1)
  layer_clim <- sample(nlayers(then.rast[[1]]), 1)
  
  taxon.vals <- ncvar_get(western, taxon, c(1,1,1), c(-1, -1, -1))
  
  values <- data.frame(x = western.grid$x,
                       y = western.grid$y, 
                       taxon = as.numeric(taxon.vals[,,layer_ncdf]))
  
  coordinates(values) <- ~ x + y
  proj4string(values) <- CRS('+init=epsg:3175')
  
  good.vals <- extract(presence, values)
  
  values <- values[good.vals == 1 & !is.na(good.vals), ]
  
  clims <- do.call(cbind.data.frame, lapply(then.rast, function(x)extract(x[[layer_clim]], values)))
  colnames(clims) <- c('bio16', 'bio17', 'bio4', 'bio5', 'bio6')
  density(clims[,clim.var], weights = values$taxon/sum(values$taxon),
          from = ranges[clim.var,1], to = ranges[clim.var,2], n = 100)$y
}

get_dens_fia <- function(taxon, clim.var){

  layer_clim <- sample(nlayers(now.rast[[1]]), 1)
  
  taxon.vals <- agg.dens[,taxon]
  
  values <- data.frame(x = xyFromCell(base.rast, agg.dens$cell)[,'x'],
                       y = xyFromCell(base.rast, agg.dens$cell)[,'y'], 
                       taxon = agg.dens[,taxon])
  
  values <- na.omit(values)
  
  coordinates(values) <- ~ x + y
  proj4string(values) <- CRS('+init=epsg:3175')
  
  good.vals <- extract(presence, values)
  
  values <- values[good.vals == 1 & !is.na(good.vals), ]
  
  clims <- do.call(cbind.data.frame, lapply(now.rast, function(x)extract(x[[layer_clim]], values)))
  colnames(clims) <- c('Pwet', 'Pdry', 'Tdiff', 'Tmax', 'Tmin')
  density(clims[,clim.var], weights = values$taxon/sum(values$taxon),
          from = ranges[clim.var,1], to = ranges[clim.var,2], n = 100)$y
}

plss.tables <-   llply(all.taxa,
                        function(x){
                          lapply(1:5, function(y){
                            do.call(rbind.data.frame, 
                              lapply(1:100, function(z)get_dens_plss(x, y)))
                          })
                        }, .progress = 'text')

fia.tables <-   llply(all.taxa,
                      function(x){
                        lapply(1:5, function(y){
                          do.call(rbind.data.frame, 
                                  lapply(1:100, function(z)get_dens_fia(x, y)))
                        })
                      }, .progress = 'text')

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
                                    face = 'bold', size = 18, colo r = 'black'),
        strip.text =  element_text(family = 'serif', 
                                   face = 'bold', size = 18, colo r = 'black'))
        

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
