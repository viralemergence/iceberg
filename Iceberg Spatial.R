
# Making range adjacency matrices ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

IcebergAdjList <- list()

PredReps <- c("Currents", paste0("Futures", 1:4))

x = 1

for(x in 1:length(PredReps)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/",PredReps[x]))
  
  VeloxList <- lapply(Files, function(a){
    print(which(Files==a))
    print(a)
    r1 <- velox(paste(paste0("Iceberg Input Files/",PredReps[x]), a, sep = '/'))
  })
  
  RasterLista <- lapply(1:length(VeloxList), function(a){
    
    print(Files[a])
    
    VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
    
  })
  
  Method = "resample"
  
  if(Method == "resample"){
    
    # Using resample (quicker but maybe doesn't actually work??) ####
    
    ### fix blank
    
    blank <- matrix(0,360*2,720*2)
    blank <- raster(blank)
    extent(blank) <- c(-180,180,-90,90)
    projection(blank) <- CRS("+proj=longlat +datum=WGS84")
    
    RasterListb <- lapply(1:length(RasterLista), function(a){
      
      print(a)
      
      testraster <- RasterLista[[a]]
      testraster <- raster::resample(testraster, blank, method = 'ngb')
      
    })
    
    RasterBrick <- raster::brick(RasterListb)
    names(RasterBrick) <- Files %>% str_remove(".tif$")
    
    # Colin says don't do this: crs(RasterBrick) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    RasterBrick <- projectRaster(RasterBrick, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                                 method = "bilinear")
    
    IcebergAdj <- PairsWisely(RasterBrick)
    
  } else {
    
    # Using rastertopolygons and SF package (which definitely works but takes ages) ####
    
    RasterListb <- lapply(1:length(RasterLista), function(a){
      
      print(Files[a])
      
      RasterLista[[a]] %>% rasterToPolygons(dissolve = T)
      
    })
    
    RasterListc <- st_as_sf(bind(RasterListb))
    
    RasterListc <- st_transform(RasterListc, 
                                "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 
    
    UpperFiles <- Files
    substr(UpperFiles, 1,1) <- substr(UpperFiles, 1,1) %>% toupper()
    
    RasterListc$Binomial <- UpperFiles
    
    IcebergRaster <- raster(RasterListc, res = 25000)
    
    RasterBrick <- fasterize(RasterListc, IcebergRaster, by = "Binomial")
    
    # The next stage ####
    IcebergAdj <- PairsWisely(RasterBrick)
    
  }
  
  save(IcebergAdj, file = paste0(paste0("Iceberg Output Files/",PredReps[x]),"IcebergRangeAdj.Rdata"))
  save(RasterBrick, file = paste0(paste0("Iceberg Output Files/",PredReps[x]),"IcebergRasterBrick.Rdata"))
  
  IcebergAdjList[[x]] <- IcebergAdj
  
}

