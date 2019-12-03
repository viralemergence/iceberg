
# Graskby ####

# Rscript "Iceberg Graskby.R"

library(sf); library(tidyverse); library(raster)

# Rscript "Final Iceberg Code/1_Iceberg Spatial.R" ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

CORES = 65

t1 <- Sys.time()

PredReps <- c("Currents", paste0("Futures", 1:4))

RCPs <- c("2.6","4.5","6.0","8.5")

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

Species <- 
  paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$")) %>% 
  sort

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Currents") %>% list.files(full.names = T)) ->
  CurrentFiles

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$")) ->
  names(CurrentFiles)

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Futures") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Futures") %>% list.files(full.names = T)) ->
  FutureFiles

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Futures") %>% list.files() %>% str_remove(".rds$")) ->
  names(FutureFiles)

CurrentCDFList <- FutureCDFList <- list()

Species <- SpeciesList %>% unlist %>% sort


Mammal_Shapes <- st_read("~/ShapeFiles")

Mammal_Shapes$Binomial = str_replace(Mammal_Shapes$binomial, " ", "_")
Mammal_Shapes <- Mammal_Shapes[order(Mammal_Shapes$Binomial),]

buffer2 <- function(r, dist) {
  
  if("raster"%in%class(r)|"RasterLayer"%in%class(r)){
    projR <- projectRaster(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    projRb <- raster::buffer(projR, dist)
    projRb <- projectRaster(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
    projRb[!is.na(projRb)] <- 1
  }
  
  if("sf"%in%class(r)){
    projR <- st_transform(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    projRb <- st_buffer(projR, dist)
    projRb <- st_transform(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
  }
  
  return(projRb)
  
}

GraskBy <- function(ShapeFile, Column, Distance = 1000){
  
  Species <- ShapeFile[[Column]] %>% unique %>% sort
  
  BufferedList <- mclapply(Species, function(a){
    
    print(a)
    
    ShapeFile[ShapeFile[[Column]] == a,] %>% buffer2(., Distance)
    
  })
  
  names(BufferedList) <- Species
  
  return(BufferedList)
  
}

IUCNBuffers <- GraskBy(Mammal_Shapes, "Binomial", Distance = 1e5)

save(IUCNBuffers, file = "Iceberg Input Files/IUCNBuffers.Rdata")

# IUCNBuffers <- GraskBy(Mammal_Shapes, "Binomial", Distance = 10^6)
