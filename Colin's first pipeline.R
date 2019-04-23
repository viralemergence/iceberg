
library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

currents <- list.files('D:/ICEBERG/RawENMs/PPM/BinaryLU',
                       pattern='.tif', full.names = TRUE)

shortnames <- list.files('D:/ICEBERG/RawENMs/PPM/BinaryLU',
                       pattern='.tif', full.names = FALSE)

  
VeloxList <- lapply(currents, velox)
  
RasterLista <- lapply(1:length(VeloxList), function(a){
    VeloxList[[a]]$as.RasterLayer(band = 1) 
})
  

blank <- matrix(0,360*2,720*2)
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")
    
RasterListb <- lapply(1:length(RasterLista), function(a){
      print(a)
      raster::resample(RasterLista[[a]], blank, method = 'ngb')
})


RasterListb <- lapply(1:length(RasterLista), function(a){
  print(a)
  raster::resample(RasterLista[[a]], blank, method = 'ngb')
})

    
RasterBrick <- raster::brick(RasterListb)
names(RasterBrick) <- shortnames %>% str_remove(".tif$")
s <- sum(RasterBrick,na.rm = FALSE)

library(rasterVis); library(maps)
colors <- colorRampPalette(brewer.pal(10,"YlGn"))
levelplot(s, col.regions=colors)
map('world',add=TRUE)


########## 

library(ncdf4); setwd("D:/ICEBERG")
ncin <- nc_open("states.nc")
get <- function(varname) {
  soilm_array <- ncvar_get(ncin,varname,start=c(1,1,1166))
  slice <- raster(t(soilm_array))
  raster::extent(slice) <- c(-180,180,-90,90)
  return(slice)
}
lutypes <- names(ncin$var)[1]
landuse2017 <- stack(lapply(lutypes,get))