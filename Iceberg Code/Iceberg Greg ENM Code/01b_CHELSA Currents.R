
# CHELSA Currents

# Rscript "~/Albersnet/Iceberg Code/Iceberg Greg ENM Code/01b_CHELSA Currents.R"

# Rscript "Iceberg Code/Iceberg Greg ENM Code/01b_CHELSA Currents.R"

# setwd(here::here())

CORES <- 45

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)
library(fs)

t1 <- Sys.time()

#print("Dropping Species!")

setwd(paste0("~/Albersnet/Iceberg Files/", "CHELSA"))

# source("~/Albersnet/Iceberg Code/Iceberg Greg ENM Code/00_Iceberg Species Dropping.R")

print("Doing Currents!")

# paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
#   list.files(full.names = T) %>% 
#   append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% list.files(full.names = T)) ->
#   FullFiles

# "Rasters" %>% 
#   list.files(full.names = T, pattern = "PPM") %>%
#   map(~list.files(.x)) %>% map(length)

# paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
#   list.files() %>% str_remove(".rds$") %>% #str_split("__") %>% map_chr(2) %>%
#   append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% 
#            list.files() %>% str_remove(".rds$")) -> # %>% str_split("__") %>% map_chr(2)) ->
#   names(FullFiles)

# PredReps <- c("Currents", paste0("Futures", 1:4))

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

CopyFolders2 <- "FinalRasters" %>% list.files()

Rep <- "presen"

# for(Rep in CopyFolders2[7:length(CopyFolders2)]){

# if(0){

# for(Rep in "presen"){

# print(Rep)

FullFiles <- "FinalRasters/" %>% paste0(Rep) %>% list.files(full.names = T)

Species <- FullFiles %>% str_split("/") %>% map_chr(last) %>% str_remove(".tif$")

names(FullFiles) <- Species

Files <- FullFiles[Species]

#     mclapply(1:length(Species), function(i){
#       
#       print(Species[i])
#       
#       RasterLista <- raster(Files[[Species[i]]])
#       
#       RasterLista <- raster::resample(RasterLista, blank, method = 'ngb')
#       
#       writeRaster(RasterLista, 
#                   filename = FullFiles[[Species[i]]],
#                   overwrite = T)
#       
#     }, mc.cores = CORES)
#     
#   }
#   
#   stop()
#   
# }

# Grid areas
AreaRaster <- raster("LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# Land Use Data
iucndat <- read.csv('IucnHabitatData.csv')
convcodes <- read.csv('IUCN_LUH_conversion_table.csv')

iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){

  Habitats %>% filter(name == a) %>% pull(LUH)

}) -> HabitatList

names(HabitatList) <- Species

# landuse2017 <- brick('landuse2017.grd')

LandUse1995 <- "LandUses/histlulc1995.grd" %>% brick

# Continents ####
print("Continents!")

ContinentRaster <- raster("continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- 
  lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# IUCN ranges for continent clipping ####
print("IUCN!")

load("~/LargeFiles/MammalStackFullMercator.Rdata")
load("IUCNBuffers.Rdata")

IUCNSp <- names(MammalStackFull) %>% intersect(Species)
MammalStackFull <- MammalStackFull[IUCNSp]

# Dispersals ####

Dispersals <- read.csv("Data for dispersal_Corrected.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

ToBuffer <- intersect(Species, Dispersals$Scientific_name)

# Adding in exception for bats ####

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder%in%c("Cetacea", "Sirenia")|
                    hFamily%in%c("Phocidae", "Odobenidae", "Otariidae")) %>% pull(Sp) ->
  MarineSp

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) ->
  BatSpecies

ToBuffer <- setdiff(ToBuffer, BatSpecies)

BatSpecies <- intersect(BatSpecies, Species)

# 01_Processing Currents ####

Rep <- "presen"

FullFiles <- "FinalRasters/" %>% paste0(Rep) %>% list.files(full.names = T)

Species <- FullFiles %>% str_split("/") %>% map_chr(last) %>% str_remove(".tif$")

names(FullFiles) <- Species

Files <- FullFiles[Species]

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

dir_create("Iceberg Input Files/GretCDF/Currents")

Processed <- paste0("Iceberg Input Files/GretCDF/","Currents") %>% 
  list.files %>% str_remove(".rds$")

ToProcess <- setdiff(Species, Processed)

# ToProcess <- Species

if(length(ToProcess)>0){
  
  print("The rasters!")
  
  mclapply(1:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    SubFiles <- Files[[Sp]]# %>% list.files(full.names = T) %>% sort
    
    # if(!Sp %in% RangeBagSp){
    #   
    #   SubFiles <- SubFiles[str_detect(SubFiles, "TP05")]
    #   
    # }
    
    # SubFiles %>% print
    
    # if(length(SubFiles)>0){
    
    RasterLista <- 
      raster(SubFiles)
    
    # 02_Resampling rasters ####
    
    # RasterLista <- raster::resample(RasterLista, blank, method = 'ngb')
    
    GretCDF <- data.frame(
      
      X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
      Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol),
      
      Climate = as.numeric(!is.na(values(RasterLista)))
      
    )# %>% mutate(PreClipClim = Climate)
    
    # IUCN Buffer clipping ####
    
    if(Sp%in%IUCNSp){
      
      if(length(nrow(IUCNBuffers[[Sp]])>0)){
        
        if(nrow(IUCNBuffers[[Sp]])>0){
          
          sf1 <- st_cast(IUCNBuffers[[Sp]], "MULTIPOLYGON")
          r1 <- fasterize::fasterize(sf1, blank)
          
          IUCNValues <- values(r1)
          
          if(sum(IUCNValues, na.rm = T)>0){
            
            GretCDF$IUCN <- 0
            GretCDF[which(!is.na(IUCNValues)),"IUCN"] <- 1
            
          }else{
            
            GretCDF$IUCN <- 1
            
          }
          
        } else{
          
          GretCDF$IUCN <- 1
          
        }
        
      } else{
        
        GretCDF$IUCN <- 1
        
      }
      
    } else{
      
      GretCDF$IUCN <- 1
      
    }
    
    GretCDF[GretCDF$IUCN == 0, c("Climate")] <- 0
    
    # Continent clipping ####
    
    if(Sp%in%IUCNSp){
      
      r1 <- MammalStackFull[[Sp]]
      SpWhich <- which(!is.na(values(r1)))
      ContinentsInhabited <- unique(values(ContinentRaster)[SpWhich])
      
      if(Sp == "Canis_lupus"){
        
        ContinentsInhabited <- 
          ContinentsInhabited %>% setdiff(1)
        
      }
      
      if(length(ContinentsInhabited)>0){
        
        GretCDF$Continent <- 0
        
        GretCDF[unlist(ContinentWhich[ContinentsInhabited]),"Continent"] <- 1
        
      }else{
        
        GretCDF$Continent <- 1
        
      }
      
    } else{
      
      GretCDF$Continent <- 1
      
    }
    
    GretCDF[GretCDF$Continent == 0, c("Climate")] <- 0
    
    # Land Use filters #####
    
    SpHabitat <- HabitatList[[Sp]] %>%
      as.character %>%
      str_split("[.]") %>%
      unlist %>% unique %>%
      na.omit
    
    if(length(SpHabitat)>0){

      LandUse1995[[SpHabitat]] %>% getValues ->

        ValueDF

      if(length(SpHabitat)==1){

        Habitable <- as.numeric(ValueDF==1)
        Habitable[is.na(Habitable)] <- 0

      }else{

        Habitable <- as.numeric((ValueDF %>% rowSums(na.rm = T))>0)

      }

    }else{

      Habitable <- rep(1, nrow(GretCDF))

    }
    
    GretCDF$LandUse <- Habitable
    
    # GretCDF$LandUse <- 1
    
    GretCDF$ClimateLandUse <- as.numeric(rowSums(GretCDF[,c("Climate","LandUse")])>1)
    
    # Dispersals ####
    
    if(Sp%in%ToBuffer){
      
      Dist <- Dispersals %>% filter(Scientific_name == Sp) %>% pull(disp50)*1000/50
      
      Dist20 <- Dist*20
      Dist50 <- Dist*50
      Dist80 <- Dist*80
      
      for(j in c("Climate", "ClimateLandUse")){
        
        if(sum(GretCDF[,j]) > 0){
          
          if(j == "ClimateLandUse"&all(GretCDF$Climate==GretCDF$ClimateLandUse)){
            
            GretCDF[,"Buffer20ClimateLandUse"] <- GretCDF[,"Buffer20Climate"]
            
            GretCDF[,"Buffer50ClimateLandUse"] <- GretCDF[,"Buffer50Climate"]
            
            GretCDF[,"Buffer80ClimateLandUse"] <- GretCDF[,"Buffer80Climate"]
            
          }else{
            
            r1 <- blank
            
            values(r1) <- c(1)[ifelse(GretCDF[,j] == 1, 1, NA)]
            r2 <- buffer2(r1, Dist20)
            r3 <- resample(r2, blank, method = "ngb")
            
            GretCDF[,paste0("Buffer20",j)] <- values(r3)
            GretCDF[,paste0("Buffer20",j)][is.na(GretCDF[,paste0("Buffer20",j)])] <- 0
            
            r2 <- buffer2(r1, Dist50)
            r3 <- resample(r2, blank, method = "ngb")
            
            GretCDF[,paste0("Buffer50",j)] <- values(r3)
            GretCDF[,paste0("Buffer50",j)][is.na(GretCDF[,paste0("Buffer50",j)])] <- 0
            
            r2 <- buffer2(r1, Dist80)
            r3 <- resample(r2, blank, method = "ngb")
            
            GretCDF[,paste0("Buffer80",j)] <- values(r3)
            GretCDF[,paste0("Buffer80",j)][is.na(GretCDF[,paste0("Buffer80",j)])] <- 0
            
          }
          
        }else{
          
          GretCDF[,paste0("Buffer20",j)] <- 
            GretCDF[,paste0("Buffer50",j)] <- 
            GretCDF[,paste0("Buffer80",j)] <- 0
          
        }
      }
      
    } else {
      
      GretCDF[,paste0("Buffer", 
                      rep(c("20", "50", "80"), each = 2),
                      c("Climate", "ClimateLandUse"))] <- 1
      
    }
    
    GretCDF %>% slice(-Sea) %>% 
      dplyr::select(-LandUse) %>% 
      as.matrix %>% as("dgCMatrix") %>% 
      saveRDS(file = paste0("Iceberg Input Files/GretCDF/Currents/", Sp, ".rds"))
    
    # }
    
  }, mc.preschedule = F, mc.cores = CORES)
  
  # stop()
  
  t2 <- Sys.time()
  
  print(t2 - t1)
  
}else print("None to process!")

setwd(here::here())

"Iceberg Files/CHELSA/FinalRasters" %>% dir_info(recurse = T) %>% 
  arrange(desc(modification_time)) %>% 
  pull(path) %>% length

