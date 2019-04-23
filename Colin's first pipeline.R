
library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

currents <- list.files('D:/ICEBERG/RawENMs/PPM/BinaryLU',
                       pattern='.tif', full.names = TRUE)
  
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
    
    RasterBrick <- raster::brick(RasterListb)
    names(RasterBrick) <- Files %>% str_remove(".tif$")
    