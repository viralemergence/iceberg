# Rscript "Final Iceberg Code/1_Iceberg Spatial.R" ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel); library(fs)

CORES = 20

setwd(paste0("~/Albersnet/Iceberg Files/", "CHELSA"))

PredReps <- 
  paste0("~/Albersnet/Iceberg Files/", 
         "CHELSA", "/FinalRasters") %>% list.files(pattern = "^gf")

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180, 180, -90, 90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

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

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% 
  str_remove(".rds$") ->
  names(FutureFiles)

CurrentCDFList <- FutureCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) %>% 
  intersect(names(FutureFiles))

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

IcebergAdjList <- list()

if(file.exists(paste0("Iceberg Output Files/", "CurrentsRangeAdj",".rds"))){
  
  RangeAdj <- readRDS(paste0("Iceberg Output Files/", "CurrentsRangeAdj",".rds"))
  
  IcebergAdjList$Currents <- RangeAdj
  
  # if(file.exists(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "A",".rds"))&
  #    file.exists(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "B",".rds"))){
  #   
  #   print("Loading RangeAdj Files! Phew")
  #   
  #   IcebergAdjList$A <- IcebergAdjList$C <- IcebergAdjList$B <- IcebergAdjList$D <- list()
  #   
  #   IcebergAdjList$A$Currents <- 
  #     
  #     IcebergAdjList$C$Currents <- 
  #     
  #     readRDS(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "A",".rds"))
  #   
  #   IcebergAdjList$B$Currents <- 
  #     
  #     IcebergAdjList$D$Currents <- 
  #     
  #     readRDS(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "B",".rds"))
  #   
  
} else {
  
  # Species <- Species[1:100]
  
  CurrentCDFList <- mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles[a])
    
    readRDS(CurrentFiles[[Sp]]) %>% 
      as.matrix %>% 
      as.data.frame() %>% 
      dplyr::select(Climate)
    
  }, mc.preschedule = F, mc.cores = CORES)
  
  object.size(CurrentCDFList)/(10^9)
  
  names(CurrentCDFList) <- Species
  
  # for(Pipeline in LETTERS[1:4]){
  #   
  #   print(Pipeline)
  #   
  #   IcebergAdjList[[Pipeline]] <- list()
  #   
  #   if(Pipeline == "A"){
  #     
  #     CurrentVar <- "ClimateLandUse"
  #     FuturesVar <- "BufferClimateLandUse"
  #     
  #   }
  #   
  #   if(Pipeline == "B"){
  #     
  #     CurrentVar <- "Climate"
  #     FuturesVar <- "BufferClimate"
  #     
  #   }
  #   
  #   if(Pipeline == "C"){
  #     
  #     CurrentVar <- "ClimateLandUse"
  #     FuturesVar <- "ClimateLandUse"
  #     
  #   }
  #   
  #   if(Pipeline == "D"){
  #     
  #     CurrentVar <- "Climate"
  #     FuturesVar <- "Climate"
  #     
  #   }
  #   
  #   IcebergAdjList[[Pipeline]] <- list()
  #   
  #   if(Pipeline%in%LETTERS[c(1,2)]){
  #     
  #     print("Getting values!")
  
  CurrentCDFList %>% map("Climate") %>% 
    map(function(a) a*(AreaValues[-Sea])) %>% 
    bind_cols() %>% as.data.frame() ->
    ValueDF
  
  print("Calculating overlap!")
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  dir_create("Iceberg Output Files")
  
  saveRDS(RangeAdj, file = paste0("Iceberg Output Files/", "CurrentsRangeAdj",".rds"))
  
  IcebergAdjList$Currents <- RangeAdj
  
}

# }else{
#   
#   IcebergAdjList$Currents <- RangeAdj
#   
#   #       IcebergAdjList[[Pipeline]]$Currents <- IcebergAdjList[[which(LETTERS == Pipeline)-2]]$Currents
#   #       
#   #     }
#   #   }
#   
# }

# Futures ####

print("Doing the futures!")

rm(CurrentCDFList)

# FutureCDFList <- list()
# 
# FutureCDFList[1:length(Species)] <- list(list())

FutureCDFList <- mclapply(1:length(Species), function(a){
  
  # FutureCDFList <- list()
  
  # for(a in 1:length(Species)){
  
  Sp = Species[a]
  
  print(FutureFiles[a])
  
  FutureCDFList[[a]] <-
    readRDS(FutureFiles[[Sp]]) %>% 
    as.matrix %>% 
    as.data.frame() %>% 
    dplyr::select(contains("BufferClimate"))
  
}, mc.preschedule = F, mc.cores = CORES)

names(FutureCDFList) <- Species

object.size(FutureCDFList)/(10^9)

# for(Pipeline in LETTERS[1:4]){
#   
#   FuturesVar <- ifelse(Pipeline == "A", "BufferClimateLandUse", 
#                        ifelse(Pipeline == "B", "BufferClimate", 
#                               ifelse(Pipeline == "C","ClimateLandUse", 
#                                      "Climate")))

PredReps %>% lapply(function(a){
  
  FutureCDFList %>% map(paste0("BufferClimate.", a)) %>% 
    map(function(b) b*AreaValues[-Sea]) %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  RangeAdj <- PairsWisely(ValueDF, Area = T)
  
  saveRDS(RangeAdj, file = paste0("Iceberg Output Files/", "RangeAdj.", a,".rds"))
  
  return(RangeAdj)
  
}) -> IcebergAdjList2

IcebergAdjList %<>% append(IcebergAdjList2)

names(IcebergAdjList) <- c("Currents", PredReps)

CurrentSpecies <- rownames(IcebergAdjList[[1]])

for(x in 2:length(IcebergAdjList)){
  
  NewAdj <- IcebergAdjList[[x]]
  InsertSpecies <- setdiff(CurrentSpecies, rownames(NewAdj))
  
  if(length(InsertSpecies)>0){
    
    NewAdj <- NewAdj %>% data.frame()
    NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
    NewAdj <- NewAdj %>% as.matrix
    
    IcebergAdjList[[x]] <- NewAdj[CurrentSpecies, CurrentSpecies]
  }
}

saveRDS(IcebergAdjList, file = paste0("Iceberg Output Files/","IcebergAdjList.rds"))

IcebergAdjList %>% map(nrow)

PairsDF <- 
  c("Currents", PredReps) %>% 
  
  map(function(a){
    
    DF <- IcebergAdjList[[a]] %>% reshape2::melt()
    
    names(DF)[3] <- a
    
    DF %>% return
    
  }) %>% 
  reduce(~full_join(.x, .y)) %>% 
  rename(Sp1 = Var1, Sp2 = Var2)

NonEncounters <- PairsDF %>% filter(Currents == 0)

EncounterList <- 
  c(26, 70, 85) %>% map(function(a){
    
    New11 <- NonEncounters[NonEncounters[,paste0("gf", a, 11)] > 0, c("Sp1", "Sp2")] %>% 
      dplyr::select(Sp1, Sp2)
    
    New41 <- NonEncounters[NonEncounters[,paste0("gf", a, 41)] > 0, c("Sp1", "Sp2")] %>% 
      dplyr::select(Sp1, Sp2)
    
    # New41 %<>% full_join(New11, by = c("Sp1", "Sp2"))
    
    New71 <- NonEncounters[NonEncounters[,paste0("gf", a, 71)] > 0, c("Sp1", "Sp2")] %>% 
      dplyr::select(Sp1, Sp2)
    
    # New71 %<>% full_join(New41, by = c("Sp1", "Sp2"))
    
    list(New11 = New11,
         New41 = New41,
         New71 = New71)
    
  })

EncounterList %>% saveRDS("Iceberg Output Files/EncounterList.rds")