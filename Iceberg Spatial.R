
# Parsing Cory's maps ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

IcebergAdjList <- list()

PredReps <- c("Currents", paste0("Futures", 1:2))

blank <- matrix(0,360*2,720*2) # proper resolution
# blank <- matrix(0, 360, 720) # quick and dirty resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

CORES = 1

for(x in 1:length(PredReps)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/",PredReps[x]))
  
  Velox = T
  
  if(Velox){
    
    VeloxList <- lapply(Files, function(a){
      if(which(Files==a) %% 500==0) print(a)
      print(a)
      r1 <- velox(paste(paste0("Iceberg Input Files/",PredReps[x]), a, sep = '/'))
    })
    
    RasterLista <- lapply(1:length(VeloxList), function(a){
      
      if(a %% 500==0) print(Files[a])
      
      VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
      
    })
    
  } else{
    RasterLista <- lapply(Files, function(a){
      
      if(a %% 500==0) print(a)
      
      raster(paste(paste0("Iceberg Input Files/",PredReps[x]), a, sep = '/'))
      
    })
    
  }
  
  if(x==1){
    
    RasterListb <- lapply(1:length(RasterLista), function(a){
      
      if(a %% 500==0) print(a)
      
      testraster <- RasterLista[[a]]
      testraster <- raster::resample(testraster, blank, method = 'ngb')
      
      return(testraster)
      
    })
    
  } else {
    
    RasterListb <- RasterLista
    
  }
  
  names(RasterListb) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  Project = T
  
  if(Project){
  
  RasterListc <- mclapply(RasterListb[1:round(length(RasterListb)/2)], function(a){ 
    crs(a) <- "+proj=longlat +datum=WGS84"
    projectRaster(a, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                  method = "bilinear")
  }, mc.cores = CORES)
  
  RasterListc[(length(RasterListc)+1):length(RasterListb)] <- mclapply(RasterListb[(length(RasterListc)+1):length(RasterListb)], function(a){ 
    crs(a) <- "+proj=longlat +datum=WGS84"
    projectRaster(a, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                  method = "bilinear")
  }, mc.cores = CORES)
  
  #RasterListc <- mclapply(RasterListb, function(a){
  #  crs(a) <- "+proj=longlat +datum=WGS84"
  #  projectRaster(a, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
  #                method = "bilinear")
  #}, mc.cores = 1)
  
  } else {
    RasterListc <- RasterListb
  }
  names(RasterListc) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  IcebergAdj <- PairsWisely(RasterListc)
  
  saveRDS(IcebergAdj, file = paste0(paste0("Iceberg Output Files/", PredReps[x]),"/IcebergRangeAdj.rds"))
  
  rownames(IcebergAdj) <- colnames(IcebergAdj) <- rownames(IcebergAdj) %>% str_replace('[.]',"_")
  
  IcebergAdjList[[PredReps[x]]] <- IcebergAdj
  
}

remove(list(RasterLista, RasterListb, RasterListc))

save(IcebergAdjList, file = "Iceberg Output Files/IcebergAdjList.Rdata")

RasterBrick <- raster::brick(RasterListc)

BrickAdj <- PairsWisely(RasterBrick)
save(BrickAdj, file = "TryingFutures1BrickAdj.Rdata")
ListAdj <- PairsWisely(RasterListc)
save(ListAdj, file = "TryingFutures1ListAdj.Rdata")

