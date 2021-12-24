
# 5_Iceberg Making Map Data ####

# Rscript "Final Iceberg Code/5_Iceberg Mapping.R"

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools)

PredReps <- c("Currents", paste0("Futures", 1:4))

if(CoryClimateReps[CR] == "gf"){
  
  PredReps <- c("Currents", paste0("Futures", 1:4))[c(1, 2, 4)]
  
}

PipelineReps <- LETTERS[1:4]

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

IcebergAdjList <- readRDS(paste0("Iceberg Output Files/", "IcebergAdjList.rds"))

NewEncountersList <- readRDS(paste0("Iceberg Output Files/", "NewEncounters.rds"))

AllMammaldf <- readRDS(paste0("Iceberg Output Files/", "AllMammaldf.rds"))

SpaceVars <- paste0(paste("Space", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

SharingVars <- paste0(paste("Sharing", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- 
  
  paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

# Iceberg General Maps ####

GridList <- NEGridList <- list()

CORES = 60

t1 <- Sys.time()

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

Species <- 
  paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>%
  sort

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

FutureCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))

# Adding bat exception ####

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) %>%
  intersect(Species) ->
  BatSpecies

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

PipelineReps <- LETTERS[1:4]

# Currents pipeline ####

print("Currents!")

if(file.exists("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/CurrentsGridDF.rds")){
  
  GridDF <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/CurrentsGridDF.rds")
  
} else{
  
  Sharing <- AllMammaldf %>% dplyr::select(Sp, Sp2, starts_with("Sharing.Currents"))
  
  # Importing GretCDFs ####
  
  FocalSp = Species[1]
  
  GridDF <- readRDS(CurrentFiles[[FocalSp]]) %>% 
    as.matrix %>% 
    as.data.frame() %>% 
    dplyr::select(X, Y, Climate, ClimateLandUse)
  
  SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
  
  GridDF[,paste0("Sharing.Currents",LETTERS[1:2])] <- cbind(
    sum(SubSharing$Sharing.CurrentsA)*GridDF[,"ClimateLandUse"],
    sum(SubSharing$Sharing.CurrentsB)*GridDF[,"Climate"]
  )
  
  for(a in 2:length(Species)){
    
    FocalSp = Species[a]
    
    if(a %% 1000 == 0)  print(CurrentFiles[a])
    
    SubGretCDF <- readRDS(CurrentFiles[[FocalSp]]) %>% as.matrix %>% as.data.frame()# %>% dplyr::select(Climate, ClimateLandUse)
    
    SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
    
    SubGretCDF[,paste0("Sharing.Currents",LETTERS[1:2])] <- cbind(
      sum(SubSharing$Sharing.CurrentsA)*SubGretCDF[,"ClimateLandUse"],
      sum(SubSharing$Sharing.CurrentsB)*SubGretCDF[,"Climate"]
    )
    
    GridDF[,c("Climate", "ClimateLandUse",paste0("Sharing.Currents",LETTERS[1:2]))] <- 
      GridDF[,c("Climate", "ClimateLandUse",paste0("Sharing.Currents",LETTERS[1:2]))] + 
      SubGretCDF[,c("Climate", "ClimateLandUse",paste0("Sharing.Currents",LETTERS[1:2]))]
    
  }
  
  saveRDS(GridDF, file = "~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/CurrentsGridDF.rds")
  
}

# Futures pipeline ####

print("Futures!")

FutureCDFList <- list()

Sharing <- AllMammaldf %>% 
  dplyr::select(Sp, Sp2, 
                starts_with("Sharing.Futures"))

# Importing GretCDFs ####

# Overall ####

FocalSp = Species[1]

GridDF <- readRDS(FutureFiles[[FocalSp]]) %>% 
  as.matrix %>% 
  as.data.frame()

FutureCDFList[[FocalSp]] <- GridDF %>% dplyr::select(X, Y, contains("Futures"), -starts_with("LandUse"))

SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)

for(x in 2:length(PredReps)){
  
  GridDF[,paste0("Sharing.",PredReps[x],PipelineReps)] <- cbind(
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"A")]])*GridDF[,paste0("BufferClimateLandUse.", PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"B")]])*GridDF[,paste0("BufferClimate.", PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"C")]])*GridDF[,paste0("ClimateLandUse.", PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"D")]])*GridDF[,paste0("Climate.", PredReps[x])]
    
  )
}

for(a in 2:length(Species)){
  
  FocalSp = Species[a]
  
  if(a %% 1000 == 0)  print(FutureFiles[a])
  
  SubGretCDF <- readRDS(FutureFiles[[FocalSp]]) %>% 
    as.matrix %>% 
    as.data.frame()
  
  FutureCDFList[[FocalSp]] <- 
    SubGretCDF %>% 
    dplyr::select(X, Y, contains("Futures"), -starts_with("LandUse"))
  
  SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
  
  for(x in 2:length(PredReps)){
    
    SubGretCDF[,paste0("Sharing.",PredReps[x],PipelineReps)] <- cbind(
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"A")]])*SubGretCDF[,paste0("BufferClimateLandUse.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"B")]])*SubGretCDF[,paste0("BufferClimate.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"C")]])*SubGretCDF[,paste0("ClimateLandUse.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"D")]])*SubGretCDF[,paste0("Climate.",PredReps[x])]
    )
    
  }
  
  GridDF[,2:ncol(GridDF)] <- GridDF[,2:ncol(GridDF)] + SubGretCDF[,2:ncol(SubGretCDF)]
  
}

saveRDS(GridDF, file = "Iceberg Output Files/FuturesGridDF.rds")

# No bats ####

print("No bats!")

NBSpecies <- setdiff(Species, BatSpecies)

FocalSp = NBSpecies[1]

GridDF <- FutureCDFList[[FocalSp]] 

SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)

for(x in 2:length(PredReps)){
  
  GridDF[,paste0("Sharing.",PredReps[x],PipelineReps)] <- cbind(
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"A")]])*GridDF[,paste0("BufferClimateLandUse.",PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"B")]])*GridDF[,paste0("BufferClimate.",PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"C")]])*GridDF[,paste0("ClimateLandUse.",PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"D")]])*GridDF[,paste0("Climate.",PredReps[x])]
    
  )
}

for(a in 2:length(NBSpecies)){
  
  FocalSp = NBSpecies[a]
  
  if(a %% 1000 == 0)  print(FutureFiles[a])
  
  SubGretCDF <- FutureCDFList[[FocalSp]]
  
  SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
  
  for(x in 2:length(PredReps)){
    
    SubGretCDF[,paste0("Sharing.",PredReps[x],PipelineReps)] <- cbind(
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"A")]])*SubGretCDF[,paste0("BufferClimateLandUse.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"B")]])*SubGretCDF[,paste0("BufferClimate.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"C")]])*SubGretCDF[,paste0("ClimateLandUse.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"D")]])*SubGretCDF[,paste0("Climate.",PredReps[x])]
    )
    
  }
  
  GridDF[,2:ncol(GridDF)] <- GridDF[,2:ncol(GridDF)] + SubGretCDF[,2:ncol(SubGretCDF)]
  
}

saveRDS(GridDF, file = "Iceberg Output Files/NoBatFuturesGridDF.rds")

# NewEncounters ####

print("New Encounters!")

FutureCDFList <- list()

for(i in Species){
  
  print(i)
  
  FutureCDFList[[i]] <- readRDS(FutureFiles[[i]]) %>% 
    as.matrix %>% 
    as.data.frame()
  
}

FuturesBase <- c(A = "BufferClimateLandUse",
                 B = "BufferClimate",
                 C = "ClimateLandUse",
                 D = "Climate")

# NewIntersectsManual <- list()

NewIntersectsManual <- FutureCDFList[[1]] %>% dplyr::select(X, Y)

if(file.exists("Iceberg Output Files/NewIntersects.rds")){
  
  NewIntersectsManual <- readRDS("Iceberg Output Files/NewIntersects.rds")
  
  PipelineReps <- PipelineReps %>% setdiff(c("A", "B"))

}

for(Pipeline in PipelineReps){
  
  print(Pipeline)
  
  for(x in 2:length(PredReps)){
    
    print(PredReps[x])
    
    NewEncounters <- NewEncountersList[[Pipeline]][[PredReps[x]]]
    
    FutureCDFList %>%
      map(paste0(FuturesBase[Pipeline], ".", PredReps[x])) %>% 
      bind_cols %>% as.data.frame() ->
      ValueDF
    
    OverlapSums <- OverlapSharingSums <- DeltaOverlapSharing <- rep(0, nrow(ValueDF))
    NoBatOverlapSums <- NoBatOverlapSharingSums <- NoBatDeltaOverlapSharing <- rep(0, nrow(ValueDF))
    
    for(y in 1:nrow(NewEncounters)){
      
      if(y %% 50000==0) print(y)
      
      Sp1 <- NewEncounters[y,"Sp"]
      Sp2 <- NewEncounters[y,"Sp2"]
      
      SubSums <- as.numeric(rowSums(ValueDF[,c(Sp1,Sp2)])>1)
      
      OverlapSums <- OverlapSums + SubSums
      # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
      
      OverlapSharingSums <- OverlapSharingSums +
        
        SubSums*NewEncounters[y, paste0("Sharing.",PredReps[x],Pipeline)]
      
      DeltaOverlapSharing <- DeltaOverlapSharing +
        
        SubSums*NewEncounters[y, paste0("DeltaSharing.",PredReps[x],Pipeline)]
      
      if(!Sp1%in%BatSpecies&!Sp2%in%BatSpecies){
        
        NoBatOverlapSums <- NoBatOverlapSums + SubSums
        # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
        
        NoBatOverlapSharingSums <- NoBatOverlapSharingSums +
          
          SubSums*NewEncounters[y, paste0("Sharing.",PredReps[x],Pipeline)]
        
        NoBatDeltaOverlapSharing <- NoBatDeltaOverlapSharing +
          
          SubSums*NewEncounters[y, paste0("DeltaSharing.",PredReps[x],Pipeline)]
        
        
      }
      
    }
    
    NewIntersectsManual[,paste0("Overlap.",PredReps[x],Pipeline)] <-
      OverlapSums
    
    NewIntersectsManual[,paste0("OverlapSharing.",PredReps[x],Pipeline)] <-
      OverlapSharingSums
    
    NewIntersectsManual[,paste0("DeltaOverlapSharing.",PredReps[x],Pipeline)] <-
      DeltaOverlapSharing
    
    NewIntersectsManual[,paste0("NoBatOverlap.",PredReps[x],Pipeline)] <-
      NoBatOverlapSums
    
    NewIntersectsManual[,paste0("NoBatOverlapSharing.",PredReps[x],Pipeline)] <-
      NoBatOverlapSharingSums
    
    NewIntersectsManual[,paste0("NoBatDeltaOverlapSharing.",PredReps[x],Pipeline)] <-
      NoBatDeltaOverlapSharing
    
    saveRDS(NewIntersectsManual, file = "Iceberg Output Files/NewIntersects.rds")
    
  }
}
