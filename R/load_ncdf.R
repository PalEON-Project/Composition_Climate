original.points <- read.csv('data//westerncompv0.3.csv')
presence <- data.frame(x = original.points$x,
                       y = original.points$y,
                       point = rowSums(!is.na(original.points[,3:ncol(original.points)]))>0)

coordinates(presence) <- ~ x + y
proj4string(presence) <- CRS('+init=epsg:3175')
presence <- rasterize(presence, base.rast, field = 'point')

western <- nc_open('data/input//PLScomposition_western_0.2-release.nc', write = FALSE)
#eastern <- nc_open('data/input//PLScomposition_eastern_0.2-release.nc', write = FALSE)

taxa <- unique(c(names(western$var)))#, names(eastern$var)))

#ll.ncdf <- list(long = c(ncvar_get(western, 'x'), ncvar_get(eastern, 'x')),
#                lat  = c(ncvar_get(western, 'y'), ncvar_get(eastern, 'y')))

western.grid <- expand.grid(x = ncvar_get(western, 'x'), y = rev(ncvar_get(western, 'y')))
#eastern.grid <- expand.grid(x = ncvar_get(eastern, 'x'), y = rev(ncvar_get(eastern, 'y')))

plss.density <- read.csv('data/input/plss_density_alb_v0.9-4.csv', row.names = 1)

