

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
  paste0("Iceberg Input Files/","MaxEnt","/GretCDF_NoIUCN/Currents") %>% list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF_NoIUCN/Currents") %>% list.files() %>% str_remove(".rds$")) %>% 
  sort

paste0("Iceberg Input Files/","MaxEnt","/GretCDF_NoIUCN/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF_NoIUCN/Currents") %>% list.files(full.names = T)) ->
  CurrentFiles

paste0("Iceberg Input Files/","MaxEnt","/GretCDF_NoIUCN/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF_NoIUCN/Currents") %>% list.files() %>% str_remove(".rds$")) ->
  names(CurrentFiles)

paste0("Iceberg Input Files/","MaxEnt","/GretCDF_NoIUCN/Futures") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF_NoIUCN/Futures") %>% list.files(full.names = T)) ->
  FutureFiles

paste0("Iceberg Input Files/","MaxEnt","/GretCDF_NoIUCN/Futures") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF_NoIUCN/Futures") %>% list.files() %>% str_remove(".rds$")) ->
  names(FutureFiles)

CurrentCDFList <- FutureCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

PipelineReps <- LETTERS[1:4]

Species <- names(MammalStackFull) %>% intersect(Species)

SubSums <- rep(0, length(values(blank)))

for(i in 1:length(Species)){
  
  print(i)
  
  r1 <- values(MammalStackFull[[Species[i]]])
  r1[is.na(r1)] <- 0
  SubSums <- SubSums + r1
  
}

CurrentsGridDF <- readRDS("~/Albersnet/Iceberg Output Files/CurrentsGridDF.rds")

CurrentsGridDF$IUCN <- (SubSums/AreaValues)[-Sea]

CurrentsGridDF %>% gather("Key", "Value", Climate, ClimateLandUse, IUCN) %>% ggplot(aes(X,Y,fill = Value)) + geom_tile() + facet_wrap(~Key)

CurrentsGridDF %>% gather("Key", "Value", Climate, ClimateLandUse, IUCN) %>% ggplot(aes(X,Y,fill = Value)) + geom_tile() + facet_grid(Key~.) + coord_fixed()


