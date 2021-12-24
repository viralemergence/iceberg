
# Checking the effect of IUCN Clipping

# Iceberg pre-GAM spatial processing ####

# Rscript "01_Iceberg ENM Currents.R"

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)

t1 <- Sys.time()

#print("Dropping Species!")

#source("Iceberg Code/Iceberg Greg ENM Code/00_Iceberg Species Dropping.R")

print("Doing Currents!")

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% list.files(full.names = T)) ->
  FullFiles

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>% str_split("__") %>% map_chr(2) %>%
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% 
           list.files() %>% str_remove(".rds$") %>% str_split("__") %>% map_chr(2)) ->
  names(FullFiles)

"Iceberg Input Files/GretCDF/Currents" %>% list.files %>% str_remove_all(".rds$") ->

Species

Files <- FullFiles[Species]

PredReps <- c("Currents", paste0("Futures", 1:4))

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

# Land Use Data
iucndat <- read.csv('Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Input Files/IUCN_LUH_conversion_table.csv')

iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){
  
  Habitats %>% filter(name == a) %>% pull(LUH)
  
}) -> HabitatList

names(HabitatList) <- Species

landuse2017 <- brick('Iceberg Input Files/landuse2017.grd')

# Continents ####
print("Continents!")

ContinentRaster <- raster("Iceberg Input Files/continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- 
  lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# IUCN ranges for continent clipping ####
print("IUCN!")

load("~/LargeFiles/MammalStackFullMercator.Rdata")
load("Iceberg Input Files/IUCNBuffers.Rdata")

IUCNSp <- names(MammalStackFull) %>% intersect(Species)
MammalStackFull <- MammalStackFull[IUCNSp]

# Dispersals ####

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal_Corrected.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

ToBuffer <- intersect(Species, Dispersals$Scientific_name)

# Adding in exception for bats ####

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) ->
  BatSpecies

ToBuffer <- setdiff(ToBuffer, BatSpecies)

BatSpecies <- intersect(BatSpecies, Species)

# 01_Processing Currents ####

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

Processed <- paste0("Iceberg Input Files/GretCDF/","IUCNCheck") %>% 
  list.files %>% str_remove(".rds$")

ToProcess <- setdiff(Species, Processed)

ToProcess#  <- Species

if(length(ToProcess)>0){
  
  print("The rasters!")
  
  lapply(1:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    SubFiles <- Files[[Sp]]# %>% list.files(full.names = T)
    
    if(length(SubFiles)>0){
      
      RasterLista <- raster(SubFiles[1])
      
      # 02_Resampling rasters ####
      
      RasterLista <- raster::resample(RasterLista, blank, method = 'ngb')
      
      GretCDF <- readRDS(file = paste0("Iceberg Input Files/GretCDF/Currents/", Sp, ".rds")) %>%
        as.matrix %>% as.data.frame
      
      GretCDF[,"PreClipClim"] <- values(RasterLista)[-Sea]
      
      GretCDF %>% 
        saveRDS(file = paste0("Iceberg Input Files/GretCDF/IUCNCheck/", Sp, ".rds"))
      
    }
    
  }) #, mc.preschedule = F, mc.cores = 1)
  
}
