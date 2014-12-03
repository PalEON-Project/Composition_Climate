#  Load & compile the FIA data:

fia.trans <- read.csv('data/FIA_conversion.csv', stringsAsFactors = FALSE)
fia.trans[fia.trans == ''] <- NA

#  Load FIA composition data:
get.fia.stack <- function(x, type){
  stack.files <- list.files(x, full = TRUE)
  stack.fia <- stack(stack.files)
  
  fia.vals <- getValues(stack.fia)
  name.end <- regexpr(paste('_', type, sep=''), colnames(fia.vals), fixed=TRUE)
  colnames(fia.vals) <- gsub("_$","", substr(colnames(fia.vals), 1, name.end-1))
  colnames(fia.vals) <- gsub('_', '.', colnames(fia.vals))
  
  fia.vals <- as.data.frame(fia.vals)
  
  fia.vals$x <- xyFromCell(stack.fia, 1:ncell(stack.fia))[,1]
  fia.vals$y <- xyFromCell(stack.fia, 1:ncell(stack.fia))[,2]
  
  fia.vals
}

fia.dens <- get.fia.stack('data/albers/density/', type = 'density')

re.agg <- function(x){
  # This takes the fia data (by species) and aggregates it to PalEON taxa.
  orig <- 1:nrow(x)
  
  x.test <- x[,!colnames(x) %in% c('x', 'y')]
  
  x.test$orig <- orig
  
  x.melt <- melt(x.test, na.rm=TRUE, value.name= 'value', id = c('orig'))
  
  x.melt$equiv <- fia.trans$PalEON[match(gsub('.', ' ', x.melt$variable, fixed = TRUE), 
                                         gsub('.', ' ', fia.trans$scientific, fixed = TRUE))]
  
  x.melt <- x.melt[x.melt$value > 0,]
  
  x.cast <- dcast(x.melt, orig ~ equiv, 
                  sum, drop = FALSE, 
                  value.var = 'value')
  
  x.cast$x <- x$x[x.cast$orig]
  x.cast$y <- x$y[x.cast$orig]
  
  x.cast$cell <- extract(setValues(base.rast, 1:ncell(base.rast)), x.cast[,c('x', 'y')])
  
  x.cast <- x.cast[!is.na(x.cast$cell),!is.na(colnames(x.cast))]
  x.melt_2 <- melt(x.cast[ ,!colnames(x.cast)%in% c('orig', 'x', 'y')], 
                   na.rm = TRUE, 
                   value.name= 'value', 
                   id = c('cell'))
  
  x.melt_2 <- x.melt_2[x.melt_2$value > 0,]
  
  x.cast <-  dcast(x.melt_2, 
                   cell ~ variable, 
                   mean, 
                   drop = FALSE, 
                   value.var = 'value')
  
  if(any(! unique(fia.trans$PalEON) %in% colnames(x.cast))){
    missing <- unique(fia.trans$PalEON)[!unique(fia.trans$PalEON) %in% colnames(x.cast) &
                                          !is.na(unique(fia.trans$PalEON))]
    new.cols <- length(missing)
    append.it <- matrix(ncol = new.cols, nrow = nrow(x.cast))
    colnames(append.it) <- missing
    x.cast <- cbind(x.cast, append.it)
  }
  
  x.cast[,c('cell', unique(fia.trans$PalEON)[!is.na(unique(fia.trans$PalEON))])]
}

agg.dens <- re.agg(fia.dens)
