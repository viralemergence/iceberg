# Rscript "Iceberg Code/Iceberg Greg GAMM Code/1b_CHELSA Spatial.R" ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel); library(fs)

CORES = 50

setwd("~/Albersnet/Iceberg Files/CHELSA")

dir_create("Output Files")
dir_create("Output Files/Areas")
dir_create("Output Files/Ranges")

GCMs <- 
  paste0("~/Albersnet/Iceberg Files/", 
         "CHELSA", "/FinalRasters") %>% 
  list.files() %>% 
  substr(1, 2) %>% unique %>% sort %>% setdiff("pr")

PredReps <- 
  paste0("~/Albersnet/Iceberg Files/", 
         "CHELSA", "/FinalRasters") %>% 
  list.files() %>% 
  setdiff("presen")

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180, 180, -90, 90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# PanTheria 

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder%in%c("Cetacea", "Sirenia")|
                    hFamily%in%c("Phocidae", "Odobenidae", "Otariidae")) %>% pull(Sp) ->
  MarineSp

# Grid areas
AreaRaster <- raster("LandArea.asc")
AreaValues <- raster::values(AreaRaster)

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

Species <-
  "FinalRasters" %>% dir_ls() %>% 
  map(function(a){
    
    a %>% 
      dir_ls() %>% 
      str_split("/") %>% 
      map_chr(last) %>% 
      str_remove(".tif$")
    
  }) %>% reduce(intersect) %>% 
  setdiff(MarineSp)

CurrentFiles <- CurrentFiles[Species]

CurrentCDFList <- FutureCDFList <- list()

IcebergAdjList <- list()

# Currents range adj ####

if(0){
  
  CurrentCDFList <- mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles[a])
    
    readRDS(CurrentFiles[[Sp]]) %>% 
      as.matrix %>% 
      as.data.frame() %>% 
      dplyr::select(Climate, ClimateLandUse)
    
  }, mc.preschedule = F, mc.cores = CORES)
  
  object.size(CurrentCDFList)/(10^9)
  
  names(CurrentCDFList) <- Species
  
  CurrentCDFList %>% map("Climate") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% 
    bind_cols() %>% as.data.frame() ->
    ValueDF
  
  ValueDF %>% 
    colSums(na.rm = T) %>% 
    as.data.frame %>% 
    rownames_to_column %>% 
    rename(Area = 2, Sp = 1) %>% 
    saveRDS("Output Files/Areas/CurrentClimateAreas.rds")
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  dir_create("Output Files")
  
  saveRDS(RangeAdj, file = paste0("Output Files/Ranges/", "CurrentsRangeAdj",".rds"))
  
  IcebergAdjList$Currents <- RangeAdj
  
  # Land Use ####
  
  CurrentCDFList %>% map("ClimateLandUse") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% 
    bind_cols() %>% as.data.frame() ->
    ValueDF
  
  ValueDF %>% 
    colSums(na.rm = T) %>% 
    as.data.frame %>% 
    rownames_to_column %>% 
    rename(Area = 2, Sp = 1) %>% 
    saveRDS("Output Files/Areas/CurrentClimateLandUseAreas.rds")
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  dir_create("Output Files")
  
  saveRDS(RangeAdj, file = paste0("Output Files/Ranges/", "CurrentsLandUseRangeAdj",".rds"))
  
  IcebergAdjList$CurrentsLandUse <- RangeAdj
  
}

# Futures ####

print("Doing the futures!")

# rm(CurrentCDFList)

# Looping ####

FocalGCM <- GCMs[2]

for(FocalGCM in GCMs[-2]){
  
  print(FocalGCM)
  
  paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>% 
    list.files(full.names = T) ->
    FutureFiles
  
  paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>% 
    list.files() %>% 
    str_remove(".rds$") ->
    names(FutureFiles)
  
  # FutureCDFList <- list()
  # 
  # FutureCDFList[1:length(Species)] <- list(list())
  
  if(1){
    
    FutureCDFList <- 
      
      mclapply(1:length(Species), function(a){
        
        # FutureCDFList <- list()
        
        # for(a in 1:length(Species)){
        
        Sp = Species[a]
        
        print(FutureFiles[a])
        
        # FutureCDFList[[a]] <-
        readRDS(FutureFiles[[Sp]]) %>% 
          as.matrix %>% 
          as.data.frame() %>% 
          dplyr::select(contains("BufferClimate"))
        
      }, mc.preschedule = F, mc.cores = CORES)
    
    names(FutureCDFList) <- Species
    
    # print(paste0(length(FutureCDFList[map_lgl(FutureCDFList, ~class(.x) == "try-error")]), " futurecdf try errors!"))
    # FutureCDFList <- FutureCDFList[map_lgl(FutureCDFList, ~class(.x) != "try-error")]
    # Species <- names(FutureCDFList)
    
    # object.size(FutureCDFList)/(10^9)
    
    PredReps2 <- PredReps[str_detect(PredReps, paste0("^", FocalGCM))]
    
    PredReps2 %>% lapply(function(a){
      
      FutureCDFList %>% map(paste0("BufferClimateLandUse.", a)) %>% 
        map(function(b) b*AreaValues[-Sea]) %>% bind_cols() %>% as.data.frame() ->
        ValueDF
      
      ValueDF %>% 
        colSums(na.rm = T) %>% 
        as.data.frame %>% 
        rownames_to_column %>% 
        rename(Area = 2, Sp = 1) %>% 
        saveRDS(paste0("Output Files/Areas/CLUDAreas.", a, ".rds"))
      
      RangeAdj <- PairsWisely(ValueDF, Area = T)
      
      saveRDS(RangeAdj, file = paste0("Output Files/Ranges/", "CLUDRangeAdj.", a,".rds"))
      
      return(RangeAdj)
      
    }) -> IcebergAdjList2
    
    # IcebergAdjList %<>% append(IcebergAdjList2)
    
    PredReps2 %>% lapply(function(a){
      
      FutureCDFList %>% map(paste0("BufferClimate.", a)) %>% 
        map(function(b) b*AreaValues[-Sea]) %>% bind_cols() %>% as.data.frame() ->
        ValueDF
      
      ValueDF %>% 
        colSums(na.rm = T) %>% 
        as.data.frame %>% 
        rownames_to_column %>% 
        rename(Area = 2, Sp = 1) %>% 
        saveRDS(paste0("Output Files/Areas/CDAreas.", a, ".rds"))
      
      RangeAdj <- PairsWisely(ValueDF, Area = T)
      
      saveRDS(RangeAdj, file = paste0("Output Files/Ranges/", "CDRangeAdj.", a,".rds"))
      
      return(RangeAdj)
      
    }) -> IcebergAdjList2
    
    rm(FutureCDFList)
    
    # IcebergAdjList %<>% append(IcebergAdjList2)
    
  }
  
  FutureCDFList <- 
    
    mclapply(1:length(Species), function(a){
      
      # FutureCDFList <- list()
      
      # for(a in 1:length(Species)){
      
      Sp = Species[a]
      
      print(FutureFiles[a])
      
      # FutureCDFList[[a]] <-
      readRDS(FutureFiles[[Sp]]) %>% 
        as.matrix %>% 
        as.data.frame() %>% 
        dplyr::select(matches("^Climate"))
      
    }, mc.preschedule = F, mc.cores = CORES)
  
  names(FutureCDFList) <- Species
  
  PredReps2 <- PredReps[str_detect(PredReps, paste0("^", FocalGCM))]
  
  PredReps2 %>% lapply(function(a){
    
    FutureCDFList %>% map(paste0("ClimateLandUse.", a)) %>% 
      map(function(b) b*AreaValues[-Sea]) %>% bind_cols() %>% as.data.frame() ->
      ValueDF
    
    ValueDF %>% 
      colSums(na.rm = T) %>% 
      as.data.frame %>% 
      rownames_to_column %>% 
      rename(Area = 2, Sp = 1) %>% 
      saveRDS(paste0("Output Files/Areas/CLUAreas.", a, ".rds"))
    
    RangeAdj <- PairsWisely(ValueDF, Area = T)
    
    saveRDS(RangeAdj, file = paste0("Output Files/Ranges/", "CLURangeAdj.", a,".rds"))
    
    return(RangeAdj)
    
  }) -> IcebergAdjList2
  
  # IcebergAdjList %<>% append(IcebergAdjList2)
  
  PredReps2 %>% lapply(function(a){
    
    FutureCDFList %>% map(paste0("Climate.", a)) %>% 
      map(function(b) b*AreaValues[-Sea]) %>% bind_cols() %>% as.data.frame() ->
      ValueDF
    
    ValueDF %>% 
      colSums(na.rm = T) %>% 
      as.data.frame %>% 
      rownames_to_column %>% 
      rename(Area = 2, Sp = 1) %>% 
      saveRDS(paste0("Output Files/Areas/CAreas.", a, ".rds"))
    
    RangeAdj <- PairsWisely(ValueDF, Area = T)
    
    saveRDS(RangeAdj, file = paste0("Output Files/Ranges/", "CRangeAdj.", a,".rds"))
    
    return(RangeAdj)
    
  }) -> IcebergAdjList2
  
  # IcebergAdjList %<>% append(IcebergAdjList2)
  # 
  # names(IcebergAdjList) <- c("CurrentsC", "CurrentsCLU", 
  #                            paste0(PredReps2, "CLUD"),
  #                            paste0(PredReps2, "CD"),
  #                            paste0(PredReps2, "CLU"),
  #                            paste0(PredReps2, "C"))
  # 
  # CurrentSpecies <- rownames(IcebergAdjList[[1]])
  # 
  # for(x in 2:length(IcebergAdjList)){
  #   
  #   NewAdj <- IcebergAdjList[[x]]
  #   InsertSpecies <- setdiff(CurrentSpecies, rownames(NewAdj))
  #   
  #   if(length(InsertSpecies)>0){
  #     
  #     NewAdj <- NewAdj %>% data.frame()
  #     NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
  #     NewAdj <- NewAdj %>% as.matrix
  #     
  #     IcebergAdjList[[x]] <- NewAdj[CurrentSpecies, CurrentSpecies]
  #   }
  # }
  # 
  # saveRDS(IcebergAdjList, file = paste0("Output Files/", FocalGCM, "IcebergAdjList.rds"))
  # 
  # IcebergAdjList <- IcebergAdjList[1:2]
  
}

setwd(here::here())

IcebergAdjList %>% map(nrow)

# PairsDF <- 
#   c("Currents", PredReps) %>% 
#   
#   map(function(a){
#     
#     DF <- IcebergAdjList[[a]] %>% reshape2::melt()
#     
#     names(DF)[3] <- a
#     
#     DF %>% return
#     
#   }) %>% 
#   reduce(~full_join(.x, .y)) %>% 
#   rename(Sp1 = Var1, Sp2 = Var2)
# 
# NonEncounters <- PairsDF %>% filter(Currents == 0)
# 
# EncounterList <- 
#   c(26, 70, 85) %>% map(function(a){
#     
#     New11 <- NonEncounters[NonEncounters[,paste0("gf", a, 11)] > 0, c("Sp1", "Sp2")] %>% 
#       dplyr::select(Sp1, Sp2)
#     
#     New41 <- NonEncounters[NonEncounters[,paste0("gf", a, 41)] > 0, c("Sp1", "Sp2")] %>% 
#       dplyr::select(Sp1, Sp2)
#     
#     # New41 %<>% full_join(New11, by = c("Sp1", "Sp2"))
#     
#     New71 <- NonEncounters[NonEncounters[,paste0("gf", a, 71)] > 0, c("Sp1", "Sp2")] %>% 
#       dplyr::select(Sp1, Sp2)
#     
#     # New71 %<>% full_join(New41, by = c("Sp1", "Sp2"))
#     
#     list(New11 = New11,
#          New41 = New41,
#          New71 = New71)
#     
#   })
# 
# EncounterList %>% saveRDS("Output Files/EncounterList.rds")
