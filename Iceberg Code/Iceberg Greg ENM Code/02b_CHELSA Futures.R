
# Iceberg pre-GAM spatial processing ####

# Rscript "Iceberg Code/Iceberg Greg ENM Code/02b_CHELSA Futures.R"

t1 <- Sys.time()

CORES <- 50

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr)
library(fs)

setwd(paste0("~/Albersnet/Iceberg Files/", "CHELSA"))

paste0("~/Albersnet/Iceberg Files/", 
       "CHELSA/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

# PredReps <- c("Currents", paste0("Futures", 1:4))

PredReps <- 
  paste0("~/Albersnet/Iceberg Files/", 
         "CHELSA", "/FinalRasters") %>% 
  list.files

PredReps <- PredReps[!str_detect(PredReps, "Currents")] %>% setdiff("pr")

# Blanks
blank <- matrix(0,360*2, 720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# 01_Processing Currents ####

CurrentFiles <- "Iceberg Input Files/GretCDF/Currents" %>% list.files(full.names = T)

names(CurrentFiles) <- 
  Species <-
  "Iceberg Input Files/GretCDF/Currents" %>% list.files %>% 
  str_remove(".rds$")

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

# Continents ####
print("Continents!")

ContinentRaster <- raster("continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# # Land Use Data
iucndat <- read.csv('IucnHabitatData.csv')
convcodes <- read.csv('IUCN_LUH_conversion_table.csv')

# Land
iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){
  
  Habitats %>% filter(name == a) %>% pull(LUH)
  
}) -> HabitatList

names(HabitatList) <- Species

LandUseList <- 
  "LandUses" %>% 
  dir_ls(regex = "ssp..........grd") %>% 
  map(brick)

names(LandUseList) <- 
  "LandUses" %>% 
  list.files(pattern = "ssp..........grd") %>% 
  str_remove(".grd$") %>% 
  str_remove("lulc")

GCMs <- 
  "FinalRasters" %>% 
  list.files %>% substr(1, 2) %>% 
  unique %>% 
  sort

FocalGCM <- GCMs[1]

for(FocalGCM in GCMs){
  
  print(FocalGCM)
  
  PredReps2 <- PredReps[str_detect(PredReps, FocalGCM)]
  
  dir_create(paste0("Iceberg Input Files/GretCDF/", FocalGCM))
  
  Processed <- paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>% 
    list.files %>% str_remove(".rds$")
  
  # (ToProcess <- setdiff(Species, Processed)) %>% length
  
  ToProcess <- Species
  
  Files <- paste0("Iceberg Input Files/GretCDF/", "Currents", "/", ToProcess, ".tif")
  names(Files) <- ToProcess
  
  i = 1
  
  print(paste0("To process:", length(ToProcess)))
  
  mclapply(i:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    # 02_Resampling rasters ####
    
    RasterLista <- lapply(c("presen", PredReps2), function(a){
      
      # print(a)
      
      SubFiles <- paste0("~/Albersnet/Iceberg Files/", 
                         "CHELSA/FinalRasters/", a, "/", Sp, ".tif")
      
      # if(file.exists(SubFiles)) 
      
      raster(SubFiles)
      
    })
    
    names(RasterLista) <- c("Current", PredReps2)
    
    GretCDF <- data.frame(
      
      X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
      Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol)
      
    ) 
    
    GretCDF[,paste0("Climate.", c("Current", PredReps2))] <- 
      RasterLista %>% 
      lapply(function(a) as.numeric(!is.na(values(a)))) %>% 
      bind_cols
    
    # Land Use filters #####
    
    SpHabitat <- HabitatList[[Sp]] %>%
      as.character %>%
      str_split("[.]") %>%
      unlist %>% unique %>%
      na.omit
    
    names(LandUseList) <- PredReps2
    
    lapply(PredReps2, function(a){
      
      FocalYear <- a %>% substr(5, 6)
      
      if(length(SpHabitat)>0){
        
        LandUseList[[a]][[SpHabitat]] %>% getValues ->
          
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
      
      Habitable %>% return
      
    }) %>% bind_cols ->
      GretCDF[,paste0("LandUse.", PredReps2)]
    
    # GretCDF[,paste0("LandUse.", PredReps)] <- 1
    
    # Importing currents grid ####
    
    CurrentsGretCDF <- 
      readRDS(paste0("~/Albersnet/Iceberg Files/", 
                     "CHELSA/Iceberg Input Files/GretCDF/", 
                     "Currents", "/",
                     Sp, ".rds"))
    
    GretCDF <- GretCDF %>% slice(-Sea)
    
    for(Years in c(20, 50, 80)){
      
      GretCDF[,c("Continent", paste0("Buffer", Years, c("Climate","ClimateLandUse")))] <- 
        
        CurrentsGretCDF[,c("Continent", paste0("Buffer", Years, c("Climate","ClimateLandUse")))]
      
    }
    
    GretCDF[GretCDF$Continent == 0, paste0("Climate.", PredReps)] <- 0
    
    GretCDF[,paste0("ClimateLandUse.", PredReps2)] <-

      lapply(PredReps2, function(a){

        as.numeric(rowSums(GretCDF[,paste0(c("Climate.", "LandUse."),a)])>1)

      }) %>% bind_cols
    
    PredReps2 %>% lapply(function(a){
      
      FocalYear <- a %>% substr(5, 6) %>% as.numeric %>% add(9)
      
      List1 <- lapply(paste0("Buffer",
                             c(paste0(FocalYear, "Climate"),
                               paste0(FocalYear, "ClimateLandUse"))), 
                      function(b){

                                 as.numeric(rowSums(GretCDF[, c(b,
                                                                paste0(substr(b, 9, # THIS NUMBER NEEDS TO CHANGE
                                                                              nchar(b)),
                                                                       ".",
                                                                       a))])>1)

                               })

      names(List1) <- paste0(c("BufferClimate", "BufferClimateLandUse"),
                             ".",
                             a)
      
      # as.numeric((GretCDF[,c(paste0("Climate.", a),
      #                        paste0("Buffer", FocalYear, "Climate"))] %>% 
      #               
      #               rowSums)>1) -> List1
      
      return(List1)
      
    }) %>% bind_cols() -> FillDF
    
    # names(FillDF) <- paste0("BufferClimate.", PredReps2)
    
    GretCDF %>% bind_cols(FillDF) ->
      GretCDF
    
    GretCDF %>% 
      dplyr::select(setdiff(colnames(GretCDF), colnames(CurrentsGretCDF))) %>%
      as.matrix %>% as("dgCMatrix") %>% 
      saveRDS(file = paste0("Iceberg Input Files/GretCDF/", FocalGCM, "/", Sp, ".rds"))
    
  }, mc.preschedule = F, mc.cores = CORES)
  
  t2 = Sys.time()
  
  t2 - t1
  
}

setwd(here::here())

# source("~/Albersnet/Iceberg Code/Iceberg Greg ENM Code/02b_CHELSA Futures.R")