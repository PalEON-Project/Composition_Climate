#  Loading in the Natural Earth datasets:

library(rgeos)

#  Load & clip the base raster.  This is broader than the clipping extent for the 
#  domain, and it helps speed up the use of the natural earth layers later in the
#  analysis.
#  Natural earth data is stored in the Maps folder, and comes from: www.naturalearthdata.com
#  Here we are using the 10m data.

natural_earth <- function(){
  nat.earth <- crop(stack('data/input/NaturalEarth/BaseRaster/NE1_HR_LC_SR_W_DR.tif'), y=extent(c(-98-1, -83+1, 42-1, 50+1)))
  
  xylimits <-  c(-200000, 1150000, 500000, 1700000)
  
  if(model.proj == '+init=epsg:4326'){
    ext <- c(-98, -83, 42, 50)
  }
  
  if(model.proj == '+init=epsg:3175'){
    ext <- c(-100000, 1050000, 600000, 1600000)
  }
  
  quick.subset <- function(x, y, ext){
    # Subsets a shapefile (x = file, y = layer) to a given extent.
    shape <- readOGR(x, y)
    
    shape@data$id <- rownames(shape@data)
    
    x.f = fortify(shape, region="id")
    x.join = plyr::join(x.f, shape@data, by="id")
    
    aa <- x.join[,1:2]
    coordinates(aa) <- ~ long + lat
    proj4string(aa) <- CRS(proj4string(nat.earth))
    aa <- coordinates(spTransform(aa, CRS(model.proj)))
    
    x.join[,1:2] <- aa
    
    x.subset <- subset(x.join, x.join$long > xylimits[1] & x.join$long < xylimits[2] &
                             x.join$lat > xylimits[3] & x.join$lat < xylimits[4])
    
    x.subset
  }
  
  lakes.subset <- quick.subset('data/input/NaturalEarth/Lakes/ne_10m_lakes.shp',
                               'ne_10m_lakes', ext)
  river.subset <- quick.subset('data/input/NaturalEarth/Rivers/ne_10m_rivers_lake_centerlines.shp',
                               'ne_10m_rivers_lake_centerlines', ext)
  coast.subset <- quick.subset('data/input/NaturalEarth/Coasts/ne_10m_coastline.shp',
                               'ne_10m_coastline', ext)
  
  umw.domain <- spTransform(usa, CRS(model.proj))
  can.domain <- spTransform(canada, CRS(model.proj))
  
  nat.crop <- crop(projectRaster(nat.earth, crs = model.proj), y=extent(ext))
  
  rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                           getValues(nat.crop/255))
  
  rast.table$rgb <- with(rast.table, rgb(NE1_HR_LC_SR_W_DR.1,
                                         NE1_HR_LC_SR_W_DR.2,
                                         NE1_HR_LC_SR_W_DR.3,
                                        1))
  
  list(rast_table = rast.table, 
       base_map = ggplot() + #data = rast.table, aes(x = x, y = y)) +
    #scale_alpha_discrete(range=c(1,0)) +
    geom_path(data=river.subset, aes(x = long, y = lat, group = group), color = 'blue', alpha = 0.1) +
    geom_polygon(data=lakes.subset, aes(x = long, y = lat, group = group), fill = '#ADD8E6') +
    geom_path(data = umw.domain, aes(x = long, y = lat, group = group)) +
    geom_path(data = can.domain, aes(x = long, y = lat, group = group)) +
    coord_cartesian(xlim = ext[1:2], ylim=ext[3:4]) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    xlab('') + ylab(''))
}