
# Parsing Cory's maps ####

# Rscript "Iceberg Code/Writing Projected Iceberg Rasters.R"

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

DropSpecies <- c("Lobodon carcinophaga","Leptonychotes weddellii","Enhydra lutris","Histriophoca fasciata") %>% 
  str_replace(" ","_")

IcebergAdjList <- list()

PredReps <- c("Currents", paste0("Futures", 1:4))

blank <- matrix(0,360*2,720*2) # proper resolution
# blank <- matrix(0, 360, 720) # quick and dirty resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

CORES = 1


for(x in 1:length(PredReps)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/",PredReps[x])) %>% setdiff(DropSpecies)
  
  Velox = T
  
  if(Velox){
    
    VeloxList <- lapply(Files, function(a){
      if(which(Files==a) %% 500==0) print(a)
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
  
  Project = F
  
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
    
  } else {
    RasterListc <- RasterListb
  }
  
  names(RasterListc) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  RasterListc <- RasterListc[setdiff(names(RasterListc), DropSpecies)]
  
  RasterListc %>% names %>% lapply(function(a){
    print(paste0(PredReps[x], a))
    writeRaster(RasterListc[[a]], file = paste0("Iceberg Output Files/MercatorProjUnc/",PredReps[x],"/",a,".tif"))
  })
  
}

stop()

FUCK

for(x in 1:length(PredReps)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/",PredReps[x])) %>% setdiff(DropSpecies)
  
  Velox = T
  
  if(Velox){
    
    VeloxList <- lapply(Files, function(a){
      if(which(Files==a) %% 500==0) print(a)
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
      testraster <- raster::resample(testraster, blank, method = 'ngb')*AreaRaster
      
      return(testraster)
      
    })
    
  } else {
    
    RasterListb <- lapply(RasterLista, function(a) a*AreaRaster)
    
  }
  
  names(RasterListb) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  Project = F
  
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
    
  } else {
    RasterListc <- RasterListb
  }
  
  names(RasterListc) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  RasterListc <- RasterListc[setdiff(names(RasterListc), DropSpecies)]
  
  RasterListc %>% names %>% lapply(function(a){
    print(paste0(PredReps[x], a))
    writeRaster(RasterListc[[a]], file = paste0("Iceberg Output Files/MercatorProj/",PredReps[x],"/",a,".tif"))
  })
  
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
    
  } else {
    RasterListc <- RasterListb
  }
  names(RasterListc) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  RasterListc <- RasterListc[setdiff(names(RasterListc), DropSpecies)]
  
  RasterListc %>% names %>% lapply(function(a){
    print(paste0(PredReps[x], a))
    writeRaster(RasterListc[[a]], file = paste0("Iceberg Output Files/MollweideProj/",PredReps[x],"/",a,".tif"))
  })
  
}