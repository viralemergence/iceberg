# Rscript "Iceberg Code/X_Historical Comparison.R"

# 1_Processing rasters ####

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)

t1 <- Sys.time()

CORES <- 20

#print("Dropping Species!")

#source("Iceberg Code/Iceberg Greg ENM Code/00_Iceberg Species Dropping.R")

print("Doing Currents!")

paste0("Iceberg Files/","Historical/HistoricalRasters") %>% 
  list.files(full.names = T) -> 
  FullFiles

paste0("Iceberg Files/","Historical/HistoricalRasters") %>% 
  list.files() %>% str_remove(".tif$") %>% 
  str_split("__") %>% map_chr(2) %>% 
  str_remove_all("_mean_noBias_1e.06_0_1e.06_0_all_all_none_all_maxnet_none_equalWeights$") ->
  names(FullFiles)

Species <- 
  names(FullFiles) %>% 
  intersect("Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents" %>% 
              list.files %>% str_remove(".rds$"))

Files <- FullFiles[Species]

PredReps <- c("Currents", paste0("Futures", 1:4))

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Files/Climate1/Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("Iceberg Files/Climate1/Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# Land Use Data
iucndat <- read.csv('Iceberg Files/Climate1/Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Files/Climate1/Iceberg Input Files/IUCN_LUH_conversion_table.csv')

iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){
  
  Habitats %>% filter(name == a) %>% pull(LUH)
  
}) -> HabitatList

names(HabitatList) <- Species

# landuse2017 <- brick('Iceberg Files/Climate1/Iceberg Input Files/landuse2017.grd')

# landuse2015 <- brick('Iceberg Files/Historical/landuse2015era.tif')

# landuse1991 <- brick('Iceberg Files/Historical/landuse1991era.tif')

landuse1991 <- brick("Iceberg Files/Historical/landuse1991era.grd")

# Continents ####

print("Continents!")

ContinentRaster <- raster("Iceberg Files/Climate1/Iceberg Input Files/continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- 
  lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)==a))

names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# IUCN ranges for continent clipping ####
print("IUCN!")

load("~/LargeFiles/MammalStackFullMercator.Rdata")
load("Iceberg Files/Climate1/Iceberg Input Files/IUCNBuffers.Rdata")

IUCNSp <- names(MammalStackFull) %>% intersect(Species)
MammalStackFull <- MammalStackFull[IUCNSp]

# Dispersals ####

Dispersals <- read.csv("Iceberg Files/Climate1/Iceberg Input Files/Data for dispersal_Corrected.csv", header = T)

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

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

Processed <- paste0("Iceberg Files/Historical/GretCDF") %>% 
  list.files %>% str_remove(".rds$")

# ToProcess <- setdiff(Species, Processed)

ToProcess <- Species

if(length(ToProcess)>0){
  
  print("The rasters!")
  
  mclapply(1:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    SubFiles <- Files[[Sp]]# %>% list.files(full.names = T)
    
    if(length(SubFiles)>0){
      
      RasterLista <- raster(SubFiles[1])
      
      # 02_Resampling rasters ####
      
      RasterLista <- raster::resample(RasterLista, blank, method = 'ngb')
      
      GretCDF <- data.frame(
        X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
        Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol),
        
        Climate = as.numeric(!is.na(values(RasterLista)))
        
      ) %>% mutate(PreClipClim = Climate)
      # 
      # # IUCN Buffer clipping ####
      # 
      # if(Sp%in%IUCNSp){
      #   
      #   if(nrow(IUCNBuffers[[Sp]])>0){
      #     
      #     sf1 <- st_cast(IUCNBuffers[[Sp]], "MULTIPOLYGON")
      #     r1 <- fasterize::fasterize(sf1, blank)
      #     
      #     IUCNValues <- values(r1)
      #     
      #     if(length(IUCNValues)>0){
      #       
      #       GretCDF$IUCN <- 0
      #       GretCDF[which(!is.na(IUCNValues)),"IUCN"] <- 1
      #       
      #     }else{
      #       
      #       GretCDF$IUCN <- 1
      #       
      #     }
      #     
      #   } else{
      #     
      #     GretCDF$IUCN <- 1
      #     
      #   }
      #   
      #   
      # } else{
      #   
      #   GretCDF$IUCN <- 1
      #   
      # }
      # 
      # GretCDF[GretCDF$IUCN == 0, c("Climate")] <- 0
      # 
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
        
        landuse1991[[SpHabitat]] %>% getValues ->
          
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
      
      GretCDF$ClimateLandUse <- as.numeric(rowSums(GretCDF[,c("Climate","LandUse")])>1)
      
      # # Dispersals ####
      
      if(Sp%in%ToBuffer){
        
<<<<<<< HEAD
        Dist <- Dispersals %>% 
          filter(Scientific_name == Sp) %>% 
          pull(disp50)*1000/2 # DIVIDE BY TWO BECAUSE HALF THE TIME
=======
        Dist <- Dispersals %>% filter(Scientific_name == Sp) %>% pull(disp50)*1000
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814
        
        j = "Climate"
        
        for(j in c("Climate","ClimateLandUse")){
          
          if(j == "ClimateLandUse"&all(GretCDF$Climate==GretCDF$ClimateLandUse)){
            
            GretCDF[,"BufferClimateLandUse"] <- GretCDF[,"BufferClimate"]
            
          }else{
            
            r1 <- blank
            
            values(r1) <- c(1)[ifelse(GretCDF[,j]==1,1,NA)]
            r2 <- buffer2(r1, Dist)
            r3 <- resample(r2, blank, method = "ngb")
            
            GretCDF[,paste0("Buffer",j)] <- values(r3)
            GretCDF[,paste0("Buffer",j)][is.na(GretCDF[,paste0("Buffer",j)])] <- 0
            
          }
          
        }
        
      } else {
        
        GretCDF[,paste0("Buffer",c("Climate", "ClimateLandUse"))] <- 1
        
      }
      
      GretCDF %>% slice(-Sea) %>% 
        as.matrix %>% as("dgCMatrix") %>% 
        saveRDS(file = paste0("Iceberg Files/Historical/GretCDF/", Sp, ".rds"))
      
    }
    
  }, mc.preschedule = F, mc.cores = 10)
  
  # stop()
  
  t2 <- Sys.time()
  
  print(t2 - t1)
  
}else print("None to process!")

# 2_Making Pairwise ####

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

CurrentCDFList <- list()

if(0){
  
<<<<<<< HEAD
  Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))
  
  CurrentFiles <- CurrentFiles[Species]
  
  CurrentCDFList <- 
=======
  1:length(Species) %>% 
  
  mclapply(function(a){
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814
    
    1:length(Species) %>% 
    
    mclapply(function(a){
      
      Sp = Species[a]
      
      print(CurrentFiles[a])
      
      readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
        as.data.frame() %>% dplyr::select(Climate, ClimateLandUse)
      
    }, mc.preschedule = F, mc.cores = CORES)
  
  object.size(CurrentCDFList)/(10^9)
  
  names(CurrentCDFList) <- Species
  
  # CurrentCDFList %>% map("Climate") %>% 
  #   map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  #   ValueDF
  # 
  # print("Calculating overlap!")
  # 
  # RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  # 
  # saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", "HistoricalRangeAdj.rds"))
  
  CurrentCDFList %>% map("ClimateLandUse") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  NewAdj <- RangeAdj
  
  InsertSpecies <- setdiff(colnames(ValueDF), rownames(NewAdj))
  
  if(length(InsertSpecies)>0){
    
<<<<<<< HEAD
    NewAdj <- NewAdj %>% data.frame()
    NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
    NewAdj <- NewAdj %>% as.matrix
    
    RangeAdj <- NewAdj[colnames(ValueDF), colnames(ValueDF)]
    
  }
  
  saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", "HistoricalLandUseRangeAdj.rds"))
  
}

# stop()

# Doing the currents ####

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck")  %>% 
  list.files(full.names = T) ->
  CurrentFiles2

Species ->
  names(CurrentFiles2)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents")  %>% 
  list.files(full.names = T) ->
  CurrentFiles3

Species ->
  names(CurrentFiles3)

CurrentCDFList <- list()

Species <- Species %>% sort %>% 
  intersect(names(CurrentFiles))  %>% 
  intersect(names(CurrentFiles2)) %>% 
  intersect(names(CurrentFiles3))
=======
    readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% dplyr::select(Climate, ClimateLandUse)
    
  }, mc.preschedule = F, mc.cores = CORES)

object.size(CurrentCDFList)/(10^9)

names(CurrentCDFList) <- Species

# CurrentCDFList %>% map("Climate") %>% 
#   map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
#   ValueDF
# 
# print("Calculating overlap!")
# 
# RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
# 
# saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", "HistoricalRangeAdj.rds"))

CurrentCDFList %>% map("ClimateLandUse") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", "HistoricalLandUseRangeAdj.rds"))
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814

stop()

# Importing the non-IUCN Clip ones ####

if(0){
  
<<<<<<< HEAD
  mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles2[a])
    
    DF <- readRDS(CurrentFiles2[[Sp]]) %>% 
      as.matrix %>% 
      as.data.frame() %>% 
      dplyr::select(PreClipClim) %>% 
      mutate(BufferPreClipClimLandUse = PreClipClim)
    
    # Clip by buffer? ####
    
    # if(HistoricalBufferClip){
    
    DF2 <- readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(BufferClimateLandUse)
    
    DF$BufferPreClipClimLandUse[!DF2==1] <- 0 # This one is non-IUCN-clip, dispersal-clip
    
    DF$IUCNClippedClimateLandUse <- 
      DF$DispersalIUCNClippedClimateLandUse <- 
      readRDS(CurrentFiles3[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(ClimateLandUse)
    
    DF$DispersalIUCNClippedClimateLandUse[!DF2==1] <- 0 # This one is IUCN-clip, dispersal-clip
    
    # }
    
    return(DF)
    
  }, mc.preschedule = F, mc.cores = CORES)

object.size(CurrentCDFList)/(10^9)

names(CurrentCDFList) <- Species

CurrentCDFList %>% map("DispersalIUCNClippedClimateLandUse") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% 
  bind_cols() %>% 
  as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                "DispersalIUCNClippedCurrentsRangeAdjLandUse.rds"))

CurrentCDFList %>% map("BufferPreClipClimLandUse") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                "DispersalPreClipRangeAdjLandUse.rds"))

# Importing the non-IUCN Clip ones ####

if(0){
  
  paste0("~/Albersnet/Iceberg Files/", 
         "Historical/GretCDF") %>% 
    list.files() %>% 
    str_remove(".rds$") %>% sort ->
    Species
  
  paste0("~/Albersnet/Iceberg Files/", 
         "Historical/GretCDF") %>% 
    list.files(full.names = T) ->
    CurrentFiles
  
  Species ->
    names(CurrentFiles)
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/IUCNCheck") %>% 
    list.files() %>% 
    str_remove(".rds$") %>% sort ->
    Species
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/IUCNCheck")  %>% 
    list.files(full.names = T) ->
    CurrentFiles2
  
  Species ->
    names(CurrentFiles2)
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/Currents") %>% 
    list.files() %>% 
    str_remove(".rds$") %>% sort ->
    Species
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/Currents")  %>% 
    list.files(full.names = T) ->
    CurrentFiles3
  
  Species ->
    names(CurrentFiles3)
  
  CurrentCDFList <- list()
  
  Species <- Species %>% sort %>% 
    intersect(names(CurrentFiles))  %>% 
    intersect(names(CurrentFiles2)) %>% 
    intersect(names(CurrentFiles3))
  
  CurrentFiles2 <- CurrentFiles2[Species]
  CurrentFiles3 <- CurrentFiles3[Species]
  
  HistoricalBufferClip <- T
  
  CurrentCDFList <- 
    
=======
  paste0("~/Albersnet/Iceberg Files/", 
         "Historical/GretCDF") %>% 
    list.files() %>% 
    str_remove(".rds$") %>% sort ->
    Species
  
  paste0("~/Albersnet/Iceberg Files/", 
         "Historical/GretCDF") %>% 
    list.files(full.names = T) ->
    CurrentFiles
  
  Species ->
    names(CurrentFiles)
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/IUCNCheck") %>% 
    list.files() %>% 
    str_remove(".rds$") %>% sort ->
    Species
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/IUCNCheck")  %>% 
    list.files(full.names = T) ->
    CurrentFiles2
  
  Species ->
    names(CurrentFiles2)
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/Currents") %>% 
    list.files() %>% 
    str_remove(".rds$") %>% sort ->
    Species
  
  paste0("~/Albersnet/Iceberg Files/",
         "Climate1/Iceberg Input Files/",
         "GretCDF/Currents")  %>% 
    list.files(full.names = T) ->
    CurrentFiles3
  
  Species ->
    names(CurrentFiles3)
  
  CurrentCDFList <- list()
  
  Species <- Species %>% sort %>% 
    intersect(names(CurrentFiles))  %>% 
    intersect(names(CurrentFiles2)) %>% 
    intersect(names(CurrentFiles3))
  
  CurrentFiles2 <- CurrentFiles2[Species]
  CurrentFiles3 <- CurrentFiles3[Species]
  
  HistoricalBufferClip <- T
  
  CurrentCDFList <- 
    
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814
    mclapply(1:length(Species), function(a){
      
      Sp = Species[a]
      
      print(CurrentFiles2[a])
      
      DF <- readRDS(CurrentFiles2[[Sp]]) %>% as.matrix %>% 
        as.data.frame() %>% dplyr::select(PreClipClim) %>% 
        mutate(BufferPreClipClim = PreClipClim)
      
      # Clip by buffer? ####
      
      # if(HistoricalBufferClip){
      
      DF2 <- readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
        as.data.frame() %>% pull(BufferClimate)
      
      DF$BufferPreClipClim[!DF2==1] <- 0
      
      DF$IUCNClippedClimate <- 
        DF$DispersalIUCNClippedClimate <- 
        readRDS(CurrentFiles3[[Sp]]) %>% as.matrix %>% 
        as.data.frame() %>% pull(Climate)
      
      DF$DispersalIUCNClippedClimate[!DF2==1] <- 0
      
      
      # }
      
      return(DF)
      
    }, mc.preschedule = F, mc.cores = CORES)
  
  object.size(CurrentCDFList)/(10^9)
  
  names(CurrentCDFList) <- Species
  
  CurrentCDFList %>% map("IUCNClippedClimate") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                  "IUCNClippedCurrentsRangeAdj.rds"))
  
  CurrentCDFList %>% map("BufferPreClipClim") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                  "PreClipRangeAdjDispersalClipped.rds"))
  
  CurrentCDFList %>% map("PreClipClim") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                  "PreClipRangeAdj.rds"))
  
  CurrentCDFList %>% map("DispersalIUCNClippedClimate") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                  "DispersalIUCNClippedClimateRangeAdj.rds"))
  
  "Iceberg Files/Historical" %>% list.files(pattern = "RangeAdj", full.names = T) %>% 
    map(readRDS) -> RangeAdjList
  
  names(RangeAdjList) <- 
    "Iceberg Files/Historical" %>% list.files(pattern = "RangeAdj")
  
  RangeAdjList %>% map(c(unlist, Prev))
  RangeAdjList %>% map(c(unlist, mean))
  
}

# 3_Making data frame ####
# Final Iceberg Code/Rscript "2_Iceberg Data Import.R" ####

Unlist1 <- function(x) unlist(x, recursive = F)

# Viral data import ###

library(igraph); library(magrittr); library(dplyr); library(ggplot2); require(RCurl); library(readr);
library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(cowplot)

"Iceberg Files/Historical" %>% list.files(pattern = "RangeAdj", full.names = T) %>% 
  map(readRDS) -> RangeAdjList

names(RangeAdjList) <- 
  "Iceberg Files/Historical" %>% list.files(pattern = "RangeAdj")

AssocsBase <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% data.frame()
HostTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/hosts.csv") %>% data.frame()
VirusTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/viruses.csv") %>% data.frame()

names(AssocsBase)[1:2] <- c("Virus", "Host")
AssocsBase <- mutate(AssocsBase, Virus = as.factor(Virus), Host = as.factor(Host))

AssocsBase2 <- AssocsBase
AssocsBase2 <- droplevels(AssocsBase[!AssocsBase$Host == "Homo_sapiens"&
                                       !AssocsBase$Virus == "Rabies_virus",])

# Making bipartite projections ####

AssocsTraits <- AssocsBase2[,1:2]

m <- table(AssocsTraits)
M <- as.matrix(m)

bipgraph <- graph.incidence(M, weighted = NULL)

Hostgraph <- bipartite.projection(bipgraph)$proj2

HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(AssocsBase2$Host)
Remove <- which(rowSums(HostAdj)==diag(HostAdj))
HostAdj <- HostAdj[-Remove,-Remove]

# Deriving metrics from the networks ####

Hosts <- data.frame(Sp = names(V(Hostgraph)),
                    Degree = degree(Hostgraph),
                    Eigenvector = eigen_centrality(Hostgraph)$vector,
                    Kcore = coreness(Hostgraph),
                    Between = betweenness(Hostgraph))

Hosts <- merge(Hosts, HostTraits, by.x = "Sp", by.y = "hHostNameFinal", all.x = T)
Hosts <- Hosts %>% dplyr::rename(hDom = hWildDomFAO)

Domestics <- Hosts[Hosts$hDom == "domestic", "Sp"]
Wildlife <- Hosts[Hosts$hDom == "wild", "Sp"]

AssocsTraits <- merge(AssocsTraits, HostTraits, by.x = "Host", by.y = "hHostNameFinal", all.x = T)

AssocsTraits$Domestic <- ifelse(AssocsTraits$Host%in%Domestics,1,0)
AssocsTraits$Wildlife <- ifelse(AssocsTraits$Host%in%Wildlife,1,0)

ZoonoticViruses <- AssocsBase %>% filter(Host == "Homo_sapiens") %>% dplyr::select(Virus)

Hosts <- Hosts %>% 
  mutate(
    Domestic = ifelse(Sp %in% Domestics, 1, 0),
    Wildlife = ifelse(Sp %in% Wildlife, 1, 0),
    hZoonosisCount = c(table(AssocsTraits[AssocsTraits$Virus%in%ZoonoticViruses$Virus,"Host"])),
    Records = c(table(AssocsTraits$Host))
  )

#devtools::install_github("gfalbery/ggregplot")
library(ggregplot); library(ggplot2); library(RColorBrewer)

ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

theme_set(theme_cowplot())

# Adding in all mammal supertree ####

if(!file.exists("Iceberg Files/Climate1/Iceberg Input Files/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
  write.csv(FullSTMatrix, file = "Iceberg Files/Climate1/Iceberg Input Files/FullSTMatrix.csv", row.names = F)
  
} else{FullSTMatrix <- as.matrix(read.csv("Iceberg Files/Climate1/Iceberg Input Files/FullSTMatrix.csv", header = T)) }

# Making Viral Associations and Polygons ####

VirusAssocs <- apply(M, 1, function(a) names(a[a>0]))

# Creating final dataset

# Replacing absent names in the full ST matrix ####

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

FinalHostNames <- 
  
  RangeAdjList[c("HistoricalLandUseRangeAdj.rds", 
                 "DispersalPreClipRangeAdjLandUse.rds")] %>% 
  # list %>%  
  map(rownames) %>% 
  # %>% map(rownames) %>% 
  append(list(rownames(HostAdj),
              colnames(FullSTMatrix))) %>% 
  reduce(intersect)

FHN <- FinalHostNames %>% setdiff(NonEutherianSp); length(FHN)

AllMammals <- 
  RangeAdjList %>% map(rownames) %>% 
  append(list(colnames(FullSTMatrix))) %>% 
  reduce(intersect)

AllMammals %<>% sort

AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

NameReplace <- c(
  "Micaelamys_namaquensis",
  "Akodon_paranaensis",
  "Bos_frontalis",
  "Bos_grunniens",
  "Bubalus_arnee", # Absent
  "Capra_hircus",
  "Hexaprotodon_liberiensis",
  "Equus_burchellii",
  "Oryzomys_alfaroi" ,
  "Oryzomys_laticeps",
  "Oryzomys_megacephalus",
  "Callithrix_argentata",
  "Miniopterus_schreibersii",
  "Myotis_ricketti",
  "Oryzomys_albigularis",
  "Ovis_aries",
  "Piliocolobus_badius",
  "Piliocolobus_rufomitratus" ,
  "Lycalopex_gymnocercus" ,
  "Rhinolophus_hildebrandtii",
  "Oryzomys_angouya",
  "Mops_condylurus",
  "Chaerephon_plicatus",
  "Chaerephon_pumilus",
  "Taurotragus_oryx")

names(NameReplace) <- AbsentHosts

rownames(FullSTMatrix) <- colnames(FullSTMatrix) <- sapply(colnames(FullSTMatrix), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

tFullSTMatrix <- 1 - 
  (FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp] - 
     min(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp]))/
  max(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp])

tSTMatrix <- tFullSTMatrix

# Going Ahead ####

rownames(Hosts) = Hosts$Sp

FinalHostNames <- FHN

FinalHostNames %>% setdiff(NonEutherianSp)

FHN <- FinalHostNames; length(FHN)

# UpperHosts <- # Removing diagonals, as they're uninformative
#   which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostMatrixdf <- data.frame(#Virus = c(HostAdj[FHN, FHN]),
  Historical = c(RangeAdjList$HistoricalLandUseRangeAdj.rds[FHN, FHN]),
<<<<<<< HEAD
  Space = c(RangeAdjList$DispersalPreClipRangeAdjLandUse.rds[FHN, FHN]),
  # Space = c(RangeAdjList$PreClipRangeAdjDispersalClipped.rds[FHN, FHN]),
=======
  Space = c(RangeAdjList$PreClipRangeAdjDispersalClipped.rds[FHN, FHN]),
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814
  Phylo = c(tFullSTMatrix[FHN, FHN]),
  Virus = c(HostAdj[FHN, FHN]),
  Sp = as.character(rep(FHN, each = length(FHN))),
  Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites"
                   #"LongMean", "LatMean")
)

HostMatrixdf[,HostMatrixVar] <- Hosts[HostMatrixdf$Sp, HostMatrixVar]
HostMatrixdf[,paste0(HostMatrixVar,".Sp2")] <- Hosts[HostMatrixdf$Sp2, HostMatrixVar]
HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] <- HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] %>% mutate(hAllZACites = 1167, hDiseaseZACites = 115)

HostMatrixdf <- HostMatrixdf %>% mutate(
  
  hOrder = Hosts[HostMatrixdf$Sp,"hOrder"],
  hFamily = Hosts[HostMatrixdf$Sp,"hFamily"],
  hDom = Hosts[HostMatrixdf$Sp,"hDom"]
  
)

HostMatrixdf$Space0 <- ifelse(HostMatrixdf$Space == 0, "No Overlap", "Overlap")
HostMatrixdf$Cites <- log(HostMatrixdf$hAllZACites + 1)
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp2 + 1)
HostMatrixdf$MinCites <- apply(HostMatrixdf[,c("hAllZACites", "hAllZACites.Sp2")],1, function(a) min(a, na.rm = T))

HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$MinDCites <- apply(HostMatrixdf[,c("hDiseaseZACites", "hDiseaseZACites.Sp2")],1, function(a) min(a, na.rm = T))
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp2 + 1)

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

Remove1 <- FinalHostMatrix %>% group_by(Sp) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% dplyr::select(Sp)
Remove2 <- FinalHostMatrix %>% group_by(Sp2) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% dplyr::select(Sp2)

Remove3 <- which(table(c((FinalHostMatrix %>% filter(Phylo < 0.25) %>% dplyr::select(Sp, Sp2))$Sp %>% as.character(),
                         (FinalHostMatrix %>% filter(Phylo < 0.25) %>% dplyr::select(Sp, Sp2))$Sp2 %>% as.character()))>20) %>% 
  names

RemoveSp <- intersect(Remove1$Sp, Remove2$Sp2)

FinalHostMatrix <- FinalHostMatrix %>% filter(!Sp%in%RemoveSp&!Sp2%in%RemoveSp)

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

FinalHostMatrix <- FinalHostMatrix %>% arrange(Sp, Sp2)

# Let's try this a second time ####

FHN <- levels(FinalHostMatrix$Sp)

HostMatrixdf <- data.frame(#Virus = c(HostAdj[FHN, FHN]),
  Historical = c(RangeAdjList$HistoricalLandUseRangeAdj.rds[FHN, FHN]),
<<<<<<< HEAD
  Space = c(RangeAdjList$DispersalPreClipRangeAdjLandUse[FHN, FHN]),
  # Space = c(RangeAdjList$PreClipRangeAdjDispersalClipped.rds[FHN, FHN]),
=======
  Space = c(RangeAdjList$PreClipRangeAdjDispersalClipped.rds[FHN, FHN]),
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814
  Phylo = c(tFullSTMatrix[FHN, FHN]),
  Virus = c(HostAdj[FHN, FHN]),
  Sp = as.character(rep(FHN, each = length(FHN))),
  Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites"
                   #"LongMean", "LatMean")
)

HostMatrixdf[,HostMatrixVar] <- Hosts[HostMatrixdf$Sp, HostMatrixVar]
HostMatrixdf[,paste0(HostMatrixVar,".Sp2")] <- Hosts[HostMatrixdf$Sp2, HostMatrixVar]
HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] <- HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] %>% mutate(hAllZACites = 1167, hDiseaseZACites = 115)

HostMatrixdf <- HostMatrixdf %>% mutate(
  
  hOrder = Hosts[HostMatrixdf$Sp,"hOrder"],
  hFamily = Hosts[HostMatrixdf$Sp,"hFamily"],
  hDom = Hosts[HostMatrixdf$Sp,"hDom"]
  
)

HostMatrixdf$Space0 <- ifelse(HostMatrixdf$Space == 0, "No Overlap", "Overlap")
HostMatrixdf$Cites <- log(HostMatrixdf$hAllZACites + 1)
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp2 + 1)
HostMatrixdf$MinCites <- apply(HostMatrixdf[,c("hAllZACites", "hAllZACites.Sp2")],1, function(a) min(a, na.rm = T))

HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$MinDCites <- apply(HostMatrixdf[,c("hDiseaseZACites", "hDiseaseZACites.Sp2")],1, function(a) min(a, na.rm = T))
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp2 + 1)

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

Remove1 <- FinalHostMatrix %>% group_by(Sp) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% dplyr::select(Sp)
Remove2 <- FinalHostMatrix %>% group_by(Sp2) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% dplyr::select(Sp2)

Remove3 <- which(table(c((FinalHostMatrix %>% filter(Phylo < 0.25) %>% dplyr::select(Sp, Sp2))$Sp %>% as.character(),
                         (FinalHostMatrix %>% filter(Phylo < 0.25) %>% dplyr::select(Sp, Sp2))$Sp2 %>% as.character()))>20) %>% 
  names

RemoveSp <- intersect(Remove1$Sp, Remove2$Sp2)

FinalHostMatrix <- FinalHostMatrix %>% filter(!Sp%in%RemoveSp&!Sp2%in%RemoveSp)

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

FinalHostMatrix <- FinalHostMatrix %>% arrange(Sp, Sp2)

AllMammals <- 
  RangeAdjList %>% map(rownames) %>% 
  append(list(colnames(FullSTMatrix))) %>% 
  reduce(intersect) %>% 
  sort

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Historical = c(RangeAdjList$HistoricalLandUseRangeAdj.rds[AllMammals, AllMammals]),
<<<<<<< HEAD
  Space = c(RangeAdjList$DispersalPreClipRangeAdjLandUse.rds[AllMammals, AllMammals]),
  # Space = c(RangeAdjList$PreClipRangeAdjDispersalClipped.rds[AllMammals, AllMammals]),
=======
  Space = c(RangeAdjList$PreClipRangeAdjDispersalClipped.rds[AllMammals, AllMammals]),
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
) %>% 
  mutate(Gz = as.numeric(Space==0)) %>% droplevels

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals, AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

N = nrow(AllMammaldf); N

# Running the model ####

# Running Frequentist GAMS

# Rscript "3_Iceberg GAMs.R"

library(mgcv); library(tidyverse); library(ggregplot); library(MASS); library(cowplot); library(colorspace)
library(ggregplot); library(parallel); library(igraph); 
library(Matrix); library(ROCR)

# source("Iceberg Greg GAMM Code/2_Iceberg Data Import.R")

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

Covar <- c("s(Phylo, by = ordered(Gz))",
           "t2(Phylo, Space, by = ordered(!Gz))",
           "MinCites",
           "Domestic",
           "Spp")

BAMList <- DataList <- PPList <- list()

r = 1

Pipeline <- "A"

print(Pipeline)

DataList[[Pipeline]] <- FinalHostMatrix[!NARows(FinalHostMatrix, "VirusBinary"),] %>% droplevels

# DataList[[Pipeline]]$Space <- DataList[[Pipeline]][,paste0("Space",Pipeline)]
DataList[[Pipeline]]$Space <- DataList[[Pipeline]]$Historical

DataList[[Pipeline]] %<>% mutate(Gz = as.numeric(Space == 0))

DataList[[Pipeline]]$Sp <- factor(DataList[[Pipeline]]$Sp, levels = sort(union(DataList[[Pipeline]]$Sp,DataList[[Pipeline]]$Sp2)))
DataList[[Pipeline]]$Sp2 <- factor(DataList[[Pipeline]]$Sp2, levels = sort(union(DataList[[Pipeline]]$Sp,DataList[[Pipeline]]$Sp2)))

DataList[[Pipeline]] <- DataList[[Pipeline]] %>% slice(order(Sp, Sp2))

MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Pipeline]]) %>% as.matrix
MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Pipeline]]) %>% as.matrix

SppMatrix = MZ1 + MZ2

DataList[[Pipeline]]$Spp <- SppMatrix
DataList[[Pipeline]]$Cites <- rowSums(log(DataList[[Pipeline]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
DataList[[Pipeline]]$MinCites <- apply(log(DataList[[Pipeline]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
DataList[[Pipeline]]$Domestic <- ifelse(rowSums(cbind(2- DataList[[Pipeline]]$hDom %>% as.factor %>% as.numeric,
                                                      2- DataList[[Pipeline]]$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)

PPList[[Pipeline]] <- list(Spp = list(rank = nlevels(DataList[[Pipeline]]$Sp), 
                                      diag(nlevels(DataList[[Pipeline]]$Sp))))

Formula = as.formula(paste0("VirusBinary", 
                            " ~ ",
                            paste(Covar, collapse = " + ")
))

BAMList[[Pipeline]] <- bam(Formula,
                           data = DataList[[Pipeline]], 
                           family = binomial(),
                           paraPen = PPList[[Pipeline]], select = T
)

save(DataList, PPList, BAMList, file = paste0("Iceberg Files/Historical/","BAMList.Rdata"))

FitList <- PostList <- DrawList <- list()

Model <- BAMList[[Pipeline]]

print(Pipeline)

# Model Checking ####

SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
SpCoef <- Model$coef[SpCoefNames]

# Effects ####

SpaceRange <- seq(from = 0,
                  to = 1,
                  length = 101) %>% 
  c(mean(DataList[[Pipeline]]$Space))

PhyloRange <- seq(from = 0,
                  to = 1,
                  length = 101)  %>% 
  c(mean(DataList[[Pipeline]]$Phylo))

FitList[[Pipeline]] <- expand.grid(Space = SpaceRange,
                                   Phylo = PhyloRange,
                                   MinCites = mean(DataList[[Pipeline]]$MinCites),
                                   Domestic = 0
) %>%
  mutate(SpaceQuantile = ifelse(Space == last(unique(Space)), "1.5%",
                                ifelse(Space == 0, "0%",
                                       ifelse(Space == 0.25, "25%",
                                              ifelse(Space == 0.5, "50%", NA)))),
         
         PhyloQuantile = ifelse(Phylo == last(unique(Phylo)), "0.1",
                                ifelse(Phylo == 0, "0",
                                       ifelse(Phylo == 0.25, "0.25",
                                              ifelse(Phylo == 0.5, "0.5", NA)))),
         Gz = as.numeric(Space==0))

FitList[[Pipeline]]$Spp <- matrix(0 , nrow = nrow(FitList[[Pipeline]]), ncol = length(SpCoef))

FitPredictions  <- predict.gam(Model, 
                               newdata = FitList[[Pipeline]], 
                               se.fit = T)

FitList[[Pipeline]][,c("Fit","Lower", "Upper")] <- logistic(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))

save(FitList, file = paste0("Iceberg Files/Historical/","FitList.Rdata"))

# Model Output Figure ####

pdf("Iceberg Files/Historical/GAMOutput.pdf", width = 9, height = 8)

plot_grid(FitList[[Pipeline]] %>% 
            filter(!is.na(SpaceQuantile)) %>%
            ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + theme_cowplot() + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
            geom_line(aes(group = as.factor(Space))) +
            labs(y = "Viral sharing probability", x = "Phylogenetic similarity", 
                 colour = "Geographic overlap", fill = "Geographic overlap") +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
            scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
            theme(legend.position = c(0.1, 0.8), 
                  legend.title = element_text(size = 10),
                  legend.background = element_rect(colour = "white")) +
            geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
          
          FitList[[Pipeline]] %>% 
            filter(!is.na(PhyloQuantile)) %>%
            ggplot(aes(Space, Fit, colour = PhyloQuantile)) +  theme_cowplot() + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
            geom_line(aes(group = as.factor(Phylo))) +
            labs(y = "Viral sharing probability", x = "Geographic overlap", 
                 colour = "Phylogenetic similarity", fill = "Phylogenetic similarity") +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
            scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
            theme(legend.position = c(0.1, 0.8), 
                  legend.title = element_text(size = 10),
                  legend.background = element_rect(colour = "white")) +
            geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01),
          
          FitList[[Pipeline]] %>% 
            filter(!Phylo == last(unique(Phylo)),
                   !Space == last(unique(Space))) %>%
            ggplot(aes(Space, Phylo)) +  theme_cowplot() + 
            geom_tile(aes(fill = Fit)) + 
            labs(x = "Geographic overlap", 
                 y = "Phylogenetic similarity",
                 fill = "Viral sharing\nprobability") +
            #ggtitle("Tensor Field") +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            theme(legend.position = "bottom",
                  legend.title = element_text(size = 10)) +
            geom_contour(aes(z = Fit), colour = "white", alpha = 0.8) + 
            metR::geom_text_contour(aes(z = Fit), colour = "white", size = 2.5, hjust = 0.5, vjust = 1.1, check_overlap = T) +
            scale_fill_continuous_sequential(palette = "ag_GrnYl",
                                             limits = c(0,1),
                                             breaks = c(0,0.5,1)),
          
          DataList[[Pipeline]] %>%
            ggplot(aes(Space, Phylo)) +  theme_cowplot() + 
            labs(x = "Geographic overlap", 
                 y = "Phylogenetic similarity") +
            #ggtitle("Data Distribution") +
            scale_fill_continuous_sequential(palette = "Heat 2", breaks = c(0:2*5)) +
            lims(x = c(0,1), y = c(0,1)) +
            coord_fixed() +
            theme(legend.position = "bottom") +
            geom_hex(aes(fill = stat(log(count)))),
          
          nrow = 2, 
          rel_heights = c(1,1.23), 
          labels = "AUTO") 

dev.off()

# 4_Iceberg Prediction ####

# CORES <- 60

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger)

# source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")

# PipelineReps <- LETTERS[1:4]

# load(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/",
#             "BAMList.Rdata"))

# IcebergAdjList <- readRDS("Iceberg Output Files/IcebergAdjList.rds")

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

tFullSTMatrix <- 1 - (FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp] - 
                        min(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp]))/
  max(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp])

# Making the prediction data frame ####

AllMammals <- reduce(map(RangeAdjList, rownames), 
                     intersect) %>% 
  intersect(rownames(FullSTMatrix)) %>%
  setdiff(NonEutherianSp) %>% 
  sort()

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals, each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals, length(AllMammals))),
  Phylo = c(tFullSTMatrix[AllMammals, AllMammals])
) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp")) %>% 
  droplevels()

# SpaceVars <- paste0(paste("Space", PredReps, sep = "."),
#                     rep(PipelineReps, each = length(PredReps)))
# 
# SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), 
#                       rep(PipelineReps, each = length(PredReps)))
# 
# names(SpaceVars) <- names(SharingVars) <- paste0(PredReps, rep(PipelineReps, each = length(PredReps)))

SpaceVars <- c("Historical", "Current")
SharingVars <- c("Sharing.Historical", "Sharing.Current")

AllMammalMatrix[,SpaceVars] <-
  RangeAdjList[c("HistoricalLandUseRangeAdj.rds", "DispersalPreClipRangeAdjLandUse.rds")] %>% 
  
  lapply(function(a){
    
    c(as.matrix(a[AllMammals,AllMammals]))
    
  }) %>% bind_cols()

UpperMammals <- which(upper.tri(tFullSTMatrix[AllMammals, AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

N = nrow(AllMammaldf)

Model <- BAMList[[1]]

SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
SpCoef <- Model$coef[SpCoefNames]

FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")

AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)

AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), 
                               log(FinalHostMatrix$MinCites.Sp2+1)), na.rm = T)

AllMammaldf$Domestic <- 0

Pipeline = "A"

print("Prediction Effects!")

Random = F

PredReps <- c("Historical", "Currents")

# for(Pipeline in PipelineReps){

print(Pipeline)

Model <- BAMList[[1]]

SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
SpCoef <- Model$coef[SpCoefNames]

x <- 1

for(x in 1:2){
  
  print(PredReps[x])
  
  AllMammaldf$Space = AllMammaldf[,SpaceVars[x]]
  AllMammaldf$Gz = as.numeric(AllMammaldf$Space==0)
  
  AllPredictions1b <- predict.bam(Model, 
                                  newdata = AllMammaldf, # %>% select(-Spp),
                                  type = "terms",
                                  exclude = "Spp")
  
  AllIntercept <- attr(AllPredictions1b, "constant")
  
  AllPredictions <- AllPredictions1b %>% as.data.frame
  
  AllPredictions[,"Intercept"] <- AllIntercept
  
  # if(Random){
  #   
  #   AllPredList <- mclapply(1:100, function(a){
  #     
  #     AllPredictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
  #       sample(SpCoef, N, replace = T)
  #     
  #     AllPredSums <- logistic(rowSums(AllPredictions))
  #     
  #     return(AllPredSums)
  #     
  #   }, mc.cores = CORES)
  #   
  #   PredSums <- AllPredList %>% bind_cols %>% rowSums
  #   
  #   AllMammaldf[,paste0(SharingVars[paste0(PredReps[x], Pipeline)])] <- PredSums/length(AllPredList)
  #   
  # }else{      
  #   
  
  PredSums <- logistic(rowSums(AllPredictions))
  
  AllMammaldf[,SharingVars[x]] <- PredSums
  
  #   
}

AllMammaldf <- AllMammaldf %>% dplyr::select(-Spp)

AllMammaldf[, paste0("Delta", "Space")] <-
  
  AllMammaldf$Current - AllMammaldf$Historical

AllMammaldf[, paste0("Delta", "Sharing")] <-
  
  AllMammaldf$Sharing.Current - AllMammaldf$Sharing.Historical

# Making new encounters ####

NewEncountersList <- 
  
  AllMammaldf[AllMammaldf[,("Historical")]==0&
                AllMammaldf[,("Current")]>0,]
<<<<<<< HEAD


# Making old encounters ####

OldEncountersList <- 
  AllMammaldf[AllMammaldf[,("Historical")]>0&
                AllMammaldf[,("Current")]==0,]

saveRDS(AllMammaldf, file = "Iceberg Files/Historical/AllMammaldf.rds")
saveRDS(NewEncountersList, file = paste0("Iceberg Files/Historical/NewEncounters.rds"))
saveRDS(OldEncountersList, file = paste0("Iceberg Files/Historical/OldEncounters.rds"))


# 5_Making Maps ####

# NewEncountersList %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique ->
# InvolvedSpecies

# Rscript "Final Iceberg Code/5_Iceberg Mapping.R"

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools)

AreaRaster <- raster("Iceberg Files/Climate1/Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# IcebergAdjList <- readRDS(paste0("Iceberg Files/Historical/", "IcebergAdjList.rds"))

NewEncountersList <- readRDS(paste0("Iceberg Files/Historical/", "NewEncounters.rds"))

AllMammaldf <- readRDS(paste0("Iceberg Files/Historical/", "AllMammaldf.rds"))

# SpaceVars <- paste0(paste("Space", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))
# SharingVars <- paste0(paste("Sharing", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

# names(SpaceVars) <- names(SharingVars) <- 
#  paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

# Iceberg General Maps ####

GridList <- NEGridList <- list()

CORES = 10

t1 <- Sys.time()

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

NewEncountersList %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique ->
  Species

paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(FutureFiles)

# Adding bat exception ####

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) %>%
  intersect(Species) ->
  BatSpecies

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

# PipelineReps <- LETTERS[1:4]

# Importing files ####

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck")  %>% 
  list.files(full.names = T) ->
  CurrentFiles2

Species ->
  names(CurrentFiles2)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents")  %>% 
  list.files(full.names = T) ->
  CurrentFiles3

Species ->
  names(CurrentFiles3)

CurrentCDFList <- list()

Species <- Species %>% sort %>% 
  intersect(names(CurrentFiles))  %>% 
  intersect(names(CurrentFiles2)) %>% 
  intersect(names(CurrentFiles3))

CurrentFiles2 <- CurrentFiles2[Species]
CurrentFiles3 <- CurrentFiles3[Species]

Species <- NewEncountersList %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique %>% 
  intersect(Species) %>% sort

HistoricalBufferClip <- T

# Insert CDFList

CurrentCDFList <- 
  
  mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles2[a])
    
    DF <- readRDS(CurrentFiles2[[Sp]]) %>% 
      as.matrix %>% 
      as.data.frame() %>% 
      dplyr::select(PreClipClim) %>% 
      mutate(BufferPreClipClimLandUse = PreClipClim)
    
    # Clip by buffer? ####
    
    # if(HistoricalBufferClip){
    
    DF2 <- readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(BufferClimateLandUse)
    
    DF$BufferPreClipClimLandUse[!DF2==1] <- 0 # This one is non-IUCN-clip, dispersal-clip
    
    DF$IUCNClippedClimateLandUse <- 
      DF$DispersalIUCNClippedClimateLandUse <- 
      readRDS(CurrentFiles3[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(ClimateLandUse)
    
    DF$DispersalIUCNClippedClimateLandUse[!DF2==1] <- 0 # This one is IUCN-clip, dispersal-clip
    
    # }
    
    return(DF)
    
  }, mc.preschedule = F, mc.cores = CORES)

object.size(CurrentCDFList)/(10^9)
=======

>>>>>>> cfb2e9804a2e253608955586fde1c95715142814

names(CurrentCDFList) <- Species

<<<<<<< HEAD
# New Encounters ####

NewIntersectsManual <- 
  # CurrentCDFList[[1]] %>% dplyr::select(X, Y)
  readRDS(CurrentFiles[[1]]) %>% as.matrix %>% 
  as.data.frame() %>% dplyr::select(X, Y)

# if(file.exists("Iceberg Files/Historical/NewIntersects.rds")){
#   
#   NewIntersectsManual <- readRDS("Iceberg Output Files/NewIntersects.rds")
#   
#   PipelineReps <- PipelineReps %>% setdiff(c("A", "B"))
#   
# }

NewEncounters <- NewEncountersList#[[Pipeline]][[PredReps[x]]]

CurrentCDFList %>%
  map("BufferPreClipClimLandUse") %>% 
  bind_cols %>% as.data.frame() ->
  ValueDF

ValueDF %<>% 
  mutate_all(~ifelse(is.na(.x), 0, .x))

names(ValueDF) <- Species

OverlapSums <- OverlapSharingSums <- DeltaOverlapSharing <- rep(0, nrow(ValueDF))
NoBatOverlapSums <- NoBatOverlapSharingSums <- NoBatDeltaOverlapSharing <- rep(0, nrow(ValueDF))

y = 1

for(y in y:nrow(NewEncounters)){
  
  if(y %% 10000==0) print(y)
  
  Sp1 <- NewEncounters[y,"Sp"]
  Sp2 <- NewEncounters[y,"Sp2"]
  
  SubSums <- as.numeric(rowSums(ValueDF[,c(Sp1,Sp2)])>1)
  
  OverlapSums <- OverlapSums + SubSums
  # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
  
  OverlapSharingSums <- OverlapSharingSums +
    
    SubSums*NewEncounters[y, paste0("Sharing")]
  
  DeltaOverlapSharing <- DeltaOverlapSharing +
    
    SubSums*NewEncounters[y, paste0("DeltaSharing")]
  
  if(!Sp1%in%BatSpecies&!Sp2%in%BatSpecies){
    
    print("Bat!")
    
    NoBatOverlapSums <- NoBatOverlapSums + SubSums
    # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
    
    NoBatOverlapSharingSums <- NoBatOverlapSharingSums +
      
      SubSums*NewEncounters[y, paste0("Sharing.",PredReps[x],Pipeline)]
    
    NoBatDeltaOverlapSharing <- NoBatDeltaOverlapSharing +
      
      SubSums*NewEncounters[y, paste0("DeltaSharing.",PredReps[x],Pipeline)]
    
    
  }
  
}

NewIntersectsManual[,paste0("Overlap")] <-
  OverlapSums

NewIntersectsManual[,paste0("OverlapSharing")] <-
  OverlapSharingSums

NewIntersectsManual[,paste0("DeltaOverlapSharing")] <-
  DeltaOverlapSharing

NewIntersectsManual[,paste0("NoBatOverlap.",PredReps[x],Pipeline)] <-
  NoBatOverlapSums

NewIntersectsManual[,paste0("NoBatOverlapSharing.",PredReps[x],Pipeline)] <-
  NoBatOverlapSharingSums

NewIntersectsManual[,paste0("NoBatDeltaOverlapSharing.",PredReps[x],Pipeline)] <-
  NoBatDeltaOverlapSharing

saveRDS(NewIntersectsManual, file = "Iceberg Files/Historical/NewIntersects.rds")

library(cowplot); library(colorspace); library(patchwork)
library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot);
library(RColorBrewer); library(fs); library(colorspace); library(glue); library(mgcv)
library(patchwork); library(ggregplot); library(ggplot2); library(RColorBrewer)

theme_set(theme_cowplot())

(NewIntersectsManual %>% 
    ggplot(aes(X, Y, fill = Overlap)) +
    geom_tile() + coord_sf() + 
    scale_fill_continuous_sequential(AlberPalettes[[1]])) /
  
  (NewIntersectsManual %>% 
     ggplot(aes(X, Y, fill = DeltaOverlapSharing)) +
     geom_tile() + coord_sf() + 
     scale_fill_continuous_sequential(AlberPalettes[[1]]))

UniversalBlank <- raster("Iceberg Files/Climate1/Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

DarkSpectral <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))

GregRast <- function(DF, x, Projection = NULL) {
  
  r <- UniversalBlank
  
  values(r)[Sea] <- NA
  values(r)[-Sea] <- DF[,x]
  
  r
  
}

conts <- raster('Iceberg Files/Climate1/Iceberg Input Files/continents-madagascar.tif')
africa <- conts; africa[!(africa == 1)] <- NA; africa <- africa-1

africashp <- rasterToPolygons(africa, dissolve = TRUE)
contsshp <- rasterToPolygons(conts, dissolve = TRUE)

pdf("Iceberg Files/Historical/Historical.pdf", width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(0, 0.5, 0, 5))

NewIntersectsManual %>% 
  GregRast("Overlap") %>% 
  plot(axes = FALSE, box = FALSE,
       legend.shrink = 0.8, 
       legend.width = 2, 
       axis.args = list(cex.axis = 1.2),
       col = DarkSpectral,
       legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

NewIntersectsManual %>% 
  GregRast("DeltaOverlapSharing") %>% 
  plot(axes = FALSE, box = FALSE, legend.shrink = 0.8, 
       #legend.shrink = 0.8, 
       legend.width = 2, 
       axis.args = list(cex.axis = 1.2),
       col = DarkSpectral,
       legend.args=list(text='Viral sharing events', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

dev.off()


# Stats Colin wants ####

# i'm thinking probably just number first encounters +/- sd

NewEncounters %>% #filter(Sp %in% BatSpecies|Sp2 %in% BatSpecies) %>% 
  nrow()

# and maybe # viral sharings?

NewEncounters$DeltaSharing %>% sum

NewEncounters$DeltaSpace %>% sum

# maybe # first encounter bat-anything

NewEncounters %>% filter(Sp %in% BatSpecies|Sp2 %in% BatSpecies) %>% 
  nrow()

# and maybe mean % range change


=======
OldEncountersList <- 
  AllMammaldf[AllMammaldf[,("Historical")]>0&
                AllMammaldf[,("Current")]==0,]

saveRDS(AllMammaldf, file = "Iceberg Files/Historical/AllMammaldf.rds")
saveRDS(NewEncountersList, file = paste0("Iceberg Files/Historical/NewEncounters.rds"))
saveRDS(OldEncountersList, file = paste0("Iceberg Files/Historical/OldEncounters.rds"))


# 5_Making Maps ####

# NewEncountersList %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique ->
# InvolvedSpecies

# Rscript "Final Iceberg Code/5_Iceberg Mapping.R"

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools)

AreaRaster <- raster("Iceberg Files/Climate1/Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# IcebergAdjList <- readRDS(paste0("Iceberg Files/Historical/", "IcebergAdjList.rds"))

NewEncountersList <- readRDS(paste0("Iceberg Files/Historical/", "NewEncounters.rds"))

AllMammaldf <- readRDS(paste0("Iceberg Files/Historical/", "AllMammaldf.rds"))

# SpaceVars <- paste0(paste("Space", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))
# SharingVars <- paste0(paste("Sharing", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

# names(SpaceVars) <- names(SharingVars) <- 
#  paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

# Iceberg General Maps ####

GridList <- NEGridList <- list()

CORES = 10

t1 <- Sys.time()

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

NewEncountersList %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique ->
  Species

paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(FutureFiles)

# Adding bat exception ####

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) %>%
  intersect(Species) ->
  BatSpecies

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

# PipelineReps <- LETTERS[1:4]

# Importing files ####

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck")  %>% 
  list.files(full.names = T) ->
  CurrentFiles2

Species ->
  names(CurrentFiles2)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents")  %>% 
  list.files(full.names = T) ->
  CurrentFiles3

Species ->
  names(CurrentFiles3)

CurrentCDFList <- list()

Species <- Species %>% sort %>% 
  intersect(names(CurrentFiles))  %>% 
  intersect(names(CurrentFiles2)) %>% 
  intersect(names(CurrentFiles3))

CurrentFiles2 <- CurrentFiles2[Species]
CurrentFiles3 <- CurrentFiles3[Species]

Species <- NewEncountersList %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique %>% 
  intersect(Species) %>% sort

HistoricalBufferClip <- T

CurrentCDFList <- 
  
  mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles2[a])
    
    DF <- readRDS(CurrentFiles2[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% dplyr::select(PreClipClim) %>% 
      mutate(BufferPreClipClim = PreClipClim)
    
    # Clip by buffer? ####
    
    # if(HistoricalBufferClip){
    
    DF2 <- readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(BufferClimate)
    
    DF$BufferPreClipClim[!DF2==1] <- 0
    
    # DF$IUCNClippedClimate <- 
    #   DF$DispersalIUCNClippedClimate <- 
    #   readRDS(CurrentFiles3[[Sp]]) %>% as.matrix %>% 
    #   as.data.frame() %>% pull(Climate)
    # 
    # DF$DispersalIUCNClippedClimate[!DF2==1] <- 0
    
    
    # }
    
    DF %>% dplyr::select(BufferPreClipClim) %>% return
    
  }, mc.preschedule = F, mc.cores = CORES)

names(CurrentCDFList) <- Species

# NewEncounters ####

NewIntersectsManual <- #CurrentCDFList[[1]] %>% dplyr::select(X, Y)
  readRDS(CurrentFiles[[1]]) %>% as.matrix %>% 
  as.data.frame() %>% dplyr::select(X, Y)

if(file.exists("Iceberg Files/Historical/NewIntersects.rds")){
  
  NewIntersectsManual <- readRDS("Iceberg Output Files/NewIntersects.rds")
  
  PipelineReps <- PipelineReps %>% setdiff(c("A", "B"))
  
}

NewEncounters <- NewEncountersList#[[Pipeline]][[PredReps[x]]]

CurrentCDFList %>%
  # map("BufferPreClipClim") %>% 
  bind_cols %>% as.data.frame() ->
  ValueDF

ValueDF %<>% 
  mutate_all(~ifelse(is.na(.x), 0, .x))
>>>>>>> cfb2e9804a2e253608955586fde1c95715142814

names(ValueDF) <- Species

OverlapSums <- OverlapSharingSums <- DeltaOverlapSharing <- rep(0, nrow(ValueDF))
NoBatOverlapSums <- NoBatOverlapSharingSums <- NoBatDeltaOverlapSharing <- rep(0, nrow(ValueDF))

y = 1

for(y in y:nrow(NewEncounters)){
  
  if(y %% 10000==0) print(y)
  
  Sp1 <- NewEncounters[y,"Sp"]
  Sp2 <- NewEncounters[y,"Sp2"]
  
  SubSums <- as.numeric(rowSums(ValueDF[,c(Sp1,Sp2)])>1)
  
  OverlapSums <- OverlapSums + SubSums
  # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
  
  OverlapSharingSums <- OverlapSharingSums +
    
    SubSums*NewEncounters[y, paste0("Sharing")]
  
  DeltaOverlapSharing <- DeltaOverlapSharing +
    
    SubSums*NewEncounters[y, paste0("DeltaSharing")]
  
  # if(!Sp1%in%BatSpecies&!Sp2%in%BatSpecies){
  #   
  #   NoBatOverlapSums <- NoBatOverlapSums + SubSums
  #   # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
  #   
  #   NoBatOverlapSharingSums <- NoBatOverlapSharingSums +
  #     
  #     SubSums*NewEncounters[y, paste0("Sharing.",PredReps[x],Pipeline)]
  #   
  #   NoBatDeltaOverlapSharing <- NoBatDeltaOverlapSharing +
  #     
  #     SubSums*NewEncounters[y, paste0("DeltaSharing.",PredReps[x],Pipeline)]
  #   
  #   
  # }
  
}

NewIntersectsManual[,paste0("Overlap")] <-
  OverlapSums

NewIntersectsManual[,paste0("OverlapSharing")] <-
  OverlapSharingSums

NewIntersectsManual[,paste0("DeltaOverlapSharing")] <-
  DeltaOverlapSharing

# NewIntersectsManual[,paste0("NoBatOverlap.",PredReps[x],Pipeline)] <-
#   NoBatOverlapSums
# 
# NewIntersectsManual[,paste0("NoBatOverlapSharing.",PredReps[x],Pipeline)] <-
#   NoBatOverlapSharingSums
# 
# NewIntersectsManual[,paste0("NoBatDeltaOverlapSharing.",PredReps[x],Pipeline)] <-
#   NoBatDeltaOverlapSharing

saveRDS(NewIntersectsManual, file = "Iceberg Files/Historical/NewIntersects.rds")

library(cowplot); library(colorspace); library(patchwork)
library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot);
library(RColorBrewer); library(fs); library(colorspace); library(glue); library(mgcv)
library(patchwork); library(ggregplot); library(ggplot2); library(RColorBrewer)

theme_set(theme_cowplot())

(NewIntersectsManual %>% 
    ggplot(aes(X, Y, fill = Overlap)) +
    geom_tile() + coord_sf() + 
    scale_fill_continuous_sequential(AlberPalettes[[1]])) /
  
  (NewIntersectsManual %>% 
     ggplot(aes(X, Y, fill = DeltaOverlapSharing)) +
     geom_tile() + coord_sf() + 
     scale_fill_continuous_sequential(AlberPalettes[[1]]))

UniversalBlank <- raster("Iceberg Files/Climate1/Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

DarkSpectral <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))

GregRast <- function(DF, x, Projection = NULL) {
  
  r <- UniversalBlank
  
  values(r)[Sea] <- NA
  values(r)[-Sea] <- DF[,x]
  
  r
  
}

conts <- raster('Iceberg Files/Climate1/Iceberg Input Files/continents-madagascar.tif')
africa <- conts; africa[!(africa == 1)] <- NA; africa <- africa-1

africashp <- rasterToPolygons(africa, dissolve = TRUE)
contsshp <- rasterToPolygons(conts, dissolve = TRUE)

pdf("Iceberg Files/Historical/Historical.pdf", width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(0, 0.5, 0, 5))

NewIntersectsManual %>% 
  GregRast("Overlap") %>% 
  plot(axes = FALSE, box = FALSE,
       legend.shrink = 0.8, 
       legend.width = 2, 
       axis.args = list(cex.axis = 1.2),
       col = DarkSpectral,
       legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

NewIntersectsManual %>% 
  GregRast("DeltaOverlapSharing") %>% 
  plot(axes = FALSE, box = FALSE, legend.shrink = 0.8, 
       #legend.shrink = 0.8, 
       legend.width = 2, 
       axis.args = list(cex.axis = 1.2),
       col = DarkSpectral,
       legend.args=list(text='Viral sharing events', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

dev.off()
