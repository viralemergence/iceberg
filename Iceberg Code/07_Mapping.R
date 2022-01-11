
# 5_Iceberg Making Map Data ####

# Rscript "Iceberg Code/07_Mapping.R"

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools); library(fs)
library(magrittr)

setwd(here::here())
setwd("Iceberg Files/CHELSA")

PredReps <- 
  "Output Files/Predictions" %>% list.files %>% 
  # str_split("[.]") %>% map_chr(2)
  str_remove_all("Predictions|.rds$")

AreaRaster <- raster("LandArea.asc")
AreaValues <- raster::values(AreaRaster)

if(0){
  
  print("Importing adj list!")
  
  IcebergAdjList <- #readRDS("Iceberg Output Files/IcebergAdjList.rds")
    list(readRDS("Output Files/CurrentsLandUseRangeAdj.rds"),
         readRDS("Output Files/CurrentsRangeAdj.rds")) %>% 
    append(
      
      "Output Files/Ranges" %>% dir_ls() %>% map(readRDS)
      
    )
  
  names(IcebergAdjList) <- 
    c("CLUCurrents", "CCurrents") %>% 
    c("Output Files/Ranges" %>% list.files %>% 
        str_remove(".rds$") %>% 
        str_remove("RangeAdj"))
  
  CurrentSpecies <- Species <- rownames(IcebergAdjList$CCurrents)
  
  for(x in 1:length(IcebergAdjList)){
    
    NewAdj <- IcebergAdjList[[x]]
    
    InsertSpecies <- setdiff(CurrentSpecies, rownames(NewAdj))
    
    if(length(InsertSpecies)>0){
      
      NewAdj <- NewAdj %>% data.frame()
      NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
      NewAdj <- NewAdj %>% as.matrix
      
      IcebergAdjList[[x]] <- NewAdj[CurrentSpecies, CurrentSpecies]
    }
  }
  
  NewEncounterList <- 
    IcebergAdjList[3:length(IcebergAdjList)] %>% 
    names() %>% 
    map(function(a){
      
      if(substr(a, 2, 3) == "LU"){
        
        Witch <- (IcebergAdjList$CLUCurrents == 0 & IcebergAdjList[[a]] > 0) %>% 
          which() %>% 
          intersect(which(lower.tri(IcebergAdjList[[a]])))
        
      }else{
        
        Witch <- (IcebergAdjList$CCurrents == 0 & IcebergAdjList[[a]] > 0) %>% 
          which() %>% 
          intersect(which(lower.tri(IcebergAdjList[[a]])))
        
      }
      
      IcebergAdjList[[a]] %>% 
        reshape2::melt() %>% 
        slice(Witch) %>% 
        return
      
    })
  
  names(NewEncounterList) <-
    IcebergAdjList[3:length(IcebergAdjList)] %>% 
    names()
  
  saveRDS(NewEncounterList, file = "Output Files/NewEncounterList.rds")
  
  rm(IcebergAdjList)
  
}

NewEncounterList <- readRDS("Output Files/NewEncounterList.rds")

Species <- CurrentSpecies <- readRDS("Output Files/CurrentsRangeAdj.rds") %>% 
  rownames

# Iceberg General Maps ####

GridList <- NEGridList <- list()

CORES = 35

t1 <- Sys.time()

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

paste0("~/Albersnet/Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("~/Albersnet/Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

# paste0("Iceberg Input Files/GretCDF/Futures") %>% 
#   list.files(full.names = T) ->
#   FutureFiles
# 
# paste0("Iceberg Input Files/GretCDF/Futures") %>% 
#   list.files() %>% str_remove(".rds$") ->
#   names(FutureFiles)

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
# FutureFiles <- FutureFiles[Species]

# Currents pipeline ####

print("Currents!")

if(file.exists("~/Albersnet/Iceberg Files/CHELSA/Output Files/CurrentsGridDF.rds")){
  
  GridDF <- readRDS("~/Albersnet/Iceberg Files/CHELSA/Output Files/CurrentsGridDF.rds")
  
} else{
  
  SharingPairs <- #AllMammaldf %>% dplyr::select(Sp, Sp2, starts_with("Sharing.Currents"))
    readRDS("Output Files/SpeciesPairs.rds")
  
  Sharing <- #AllMammaldf %>% dplyr::select(Sp, Sp2, starts_with("Sharing.Currents"))
    readRDS("Output Files/CurrentSharing.rds")
  
  SharingValues <- SharingPairs %>% bind_cols(Sharing)
  
  # Importing GretCDFs ####
  
  FocalSp = Species[1]
  
  GridDF <- readRDS(CurrentFiles[[FocalSp]]) %>% 
    as.matrix %>% 
    as.data.frame() %>% 
    dplyr::select(X, Y, Climate, ClimateLandUse)
  
  SubSharing <- SharingValues %>% filter(Sp %in% FocalSp | 
                                           Sp2 %in% FocalSp)
  
  GridDF[,c("Sharing.CLUCurrents",  "Sharing.CCurrents")] <- cbind(
    
    sum(SubSharing$Sharing.CLUCurrents)*GridDF[,"ClimateLandUse"],
    
    sum(SubSharing$Sharing.CCurrents)*GridDF[,"Climate"]
    
  )
  
  for(a in 2:length(Species)){
    
    FocalSp = Species[a]
    
    if(a %% 100 == 0){
      
      print(CurrentFiles[a])
      
      (GridDF %>% ggplot(aes(X, Y, fill = Sharing.CLUCurrents)) + geom_tile() + coord_sf()) %>% 
        plot
      
    }
    
    SubGretCDF <- readRDS(CurrentFiles[[FocalSp]]) %>% 
      as.matrix %>% as.data.frame()# %>% dplyr::select(Climate, ClimateLandUse)
    
    SubSharing <- SharingValues %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
    
    SubGretCDF[,c("Sharing.CLUCurrents",  "Sharing.CCurrents")] <- cbind(
      
      sum(SubSharing$Sharing.CLUCurrents)*SubGretCDF[,"ClimateLandUse"],
      
      sum(SubSharing$Sharing.CCurrents)*SubGretCDF[,"Climate"]
      
    )
    
    GridDF[,c("Sharing.CLUCurrents",  "Sharing.CCurrents")] <- 
      GridDF[,c("Sharing.CLUCurrents",  "Sharing.CCurrents")] + 
      SubGretCDF[,c("Sharing.CLUCurrents",  "Sharing.CCurrents")]
    
  }
  
  saveRDS(GridDF, file = "~/Albersnet/Iceberg Files/CHELSA/Output Files/CurrentsGridDF.rds")
  
}

# stop()

# Futures pipeline ####

print("Futures!")

SharingPairs <- readRDS("Output Files/SpeciesPairs.rds")

Sharing <- "Output Files/Predictions" %>% dir_ls %>% 
  map(~.x %>% readRDS %>% 
        dplyr::select(starts_with("Sharing"))) %>% 
  bind_cols()

Rep <- PredReps[1]

GCMs <- PredReps %>% str_split("[.]") %>% map_chr(2) %>% substr(1, 2) %>% unique %>% sort

FocalGCM <- GCMs[1]

if(0){
  
  for(FocalGCM in GCMs[4]){
    
    print(FocalGCM)
    
    FutureCDFList <- list()
    
    paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>%
      list.files(full.names = T) ->
      FutureFiles
    
    paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>%
      list.files() %>% str_remove(".rds$") ->
      names(FutureFiles)
    
    PredReps2 <- PredReps[str_detect(PredReps, FocalGCM)]# %>% rev
    
    # for(FocalSp in Species){
    #   
    #   print(FocalSp)
    #   
    #   FutureCDFList[[FocalSp]] <-
    #     readRDS(FutureFiles[[FocalSp]]) %>%
    #     as.matrix %>%
    #     as.data.frame() %>%
    #     # dplyr::select(matches("Climate[.]"))# %>%
    #     # dplyr::select(matches("^Climate"))# %>%
    #     # dplyr::select(matches("^Climate[.]")) %>%
    #     # dplyr::select(matches("^ClimateLandUse"))# %>%
    #     dplyr::select(matches("BufferClimate[.]"))# %>%
    #   # dplyr::select(matches("BufferClimateLandUse"))
    #   
    # }
    
    FutureCDFList <-
      Species %>%
      mclapply(function(FocalSp){
        
        readRDS(FutureFiles[[FocalSp]]) %>%
          as.matrix %>%
          as.data.frame() %>%
          # dplyr::select(matches("Climate[.]"))# %>%
          # dplyr::select(matches("^Climate"))# %>%
          dplyr::select(matches("^Climate[.]"))# %>%
        # dplyr::select(matches("^ClimateLandUse"))# %>%
        # dplyr::select(matches("BufferClimate[.]"))# %>%
        # dplyr::select(matches("BufferClimateLandUse"))
        
      }, mc.cores = CORES)
    
    print("mclapply done!")
    
    names(FutureCDFList) <- Species
    
    FutureCDFList %<>% 
      map(function(a){
        
        a %>% 
          rename_all(~str_replace_all(.x, c("BufferClimateLandUse" = "CLUD", 
                                            "BufferClimate" = "CD", 
                                            "ClimateLandUse" = "CLU",
                                            "Climate" = "C")))
        
      })
    
    PredReps2 <- FutureCDFList[[1]] %>% names() %>% intersect(PredReps2)
    
    for(Rep in PredReps2){
      
      print(Rep)
      
      FocalSharing <- 
        Sharing[,paste0("Sharing.", Rep)] %>% 
        bind_cols(SharingPairs, Sharing = .)
      
      Pipeline <- Rep %>% str_split("[.]") %>% map_chr(1)
      
      # Importing GretCDFs ####
      
      if(!file.exists(paste0("Output Files/IntersectDFs/", Rep, ".NewIntersects.rds"))){
        
        print("New Encounters!")
        
        NewEncounters <- NewEncounterList[[Rep]]
        
        # CurrentsGridDF <- 
        #   readRDS(CurrentFiles[[FocalSp]]) %>% 
        #   as.matrix %>% 
        #   as.data.frame() %>% 
        #   dplyr::select(X, Y)
        
        FocalSp = Species[1]
        
        NewIntersectsManual <- 
          readRDS(CurrentFiles[[FocalSp]]) %>% 
          as.matrix %>% 
          as.data.frame() %>% 
          dplyr::select()
        
        FutureCDFList %>%
          map(Rep) %>% 
          bind_cols %>% as.data.frame() ->
          ValueDF
        
        OverlapSums <- OverlapSharingSums <- DeltaOverlapSharing <- rep(0, nrow(ValueDF))
        NoBatOverlapSums <- NoBatOverlapSharingSums <- NoBatDeltaOverlapSharing <- rep(0, nrow(ValueDF))
        
        for(y in 1:nrow(NewEncounters)){
          
          if(y %% 50000==0) print(y)
          
          Sp1 <- NewEncounters[y, "Var1"]
          Sp2 <- NewEncounters[y, "Var2"]
          
          SubSums <- as.numeric(rowSums(ValueDF[,c(Sp1,Sp2)])>1)
          
          OverlapSums <- OverlapSums + SubSums
          # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
          
          OverlapSharingSums <- OverlapSharingSums +
            
            SubSums*NewEncounters[y,  "value"]
          
          # DeltaOverlapSharing <- DeltaOverlapSharing +
          #   
          #   SubSums*NewEncounters[y, paste0("DeltaSharing.",Rep)]
          
          if(!Sp1%in%BatSpecies&!Sp2%in%BatSpecies){
            
            NoBatOverlapSums <- NoBatOverlapSums + SubSums
            # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
            
            NoBatOverlapSharingSums <- NoBatOverlapSharingSums +
              
              SubSums*NewEncounters[y, "value"]
            
            # NoBatDeltaOverlapSharing <- NoBatDeltaOverlapSharing +
            #   
            #   SubSums*NewEncounters[y, paste0("DeltaSharing.",Rep)]
            
            
          }
          
        }
        
        NewIntersectsManual[,paste0("Overlap.", Rep)] <-
          OverlapSums
        
        NewIntersectsManual[,paste0("OverlapSharing.", Rep)] <-
          OverlapSharingSums
        
        # NewIntersectsManual[,paste0("DeltaOverlapSharing.",Rep)] <-
        #   DeltaOverlapSharing
        
        NewIntersectsManual[,paste0("NoBatOverlap.", Rep)] <-
          NoBatOverlapSums
        
        NewIntersectsManual[,paste0("NoBatOverlapSharing.", Rep)] <-
          NoBatOverlapSharingSums
        
        # NewIntersectsManual[,paste0("NoBatDeltaOverlapSharing.",Rep)] <-
        #   NoBatDeltaOverlapSharing
        
        saveRDS(NewIntersectsManual, 
                file = paste0("Output Files/IntersectDFs/", Rep, ".NewIntersects.rds"))
        
      }
    }
    
    rm(FutureCDFList)
    
  }
  
}

for(FocalGCM in rev(GCMs)){
  
  print(FocalGCM)
  
  FutureCDFList <- list()
  
  paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>%
    list.files(full.names = T) ->
    FutureFiles
  
  paste0("Iceberg Input Files/GretCDF/", FocalGCM) %>%
    list.files() %>% str_remove(".rds$") ->
    names(FutureFiles)
  
  PredReps2 <- PredReps[str_detect(PredReps, FocalGCM)] %>% rev
  
  for(FocalSp in Species){
    
    print(FocalSp)
    
    FutureCDFList[[FocalSp]] <- 
      readRDS(FutureFiles[[FocalSp]]) %>% 
      as.matrix %>% 
      as.data.frame() %>% 
      # dplyr::select(matches("^Climate[.]"))# %>%
      # dplyr::select(matches("^ClimateLandUse"))# %>%
      dplyr::select(matches("BufferClimate[.]"))# %>%
    # dplyr::select(matches("BufferClimateLandUse"))
    
  }
  
  FutureCDFList %<>% 
    map(function(a){
      
      a %>% 
        rename_all(~str_replace_all(.x, c("BufferClimateLandUse" = "CLUD", 
                                          "BufferClimate" = "CD", 
                                          "ClimateLandUse" = "CLU",
                                          "Climate" = "C")))
    })
  
  PredReps2 <- FutureCDFList[[1]] %>% names() %>% intersect(PredReps2)
  
  for(Rep in PredReps2){
    
    print(Rep)
    
    FocalSharing <- 
      Sharing[,paste0("Sharing.", Rep)] %>% 
      bind_cols(SharingPairs, Sharing = .)
    
    Pipeline <- Rep %>% str_split("[.]") %>% map_chr(1)
    
    # Overall ####
    
    if(!file.exists(paste0("Output Files/GridDFs/", Rep, ".GridDF.rds"))){
      
      FocalSp = Species[1]
      
      # GridDF <- readRDS(FutureFiles[[FocalSp]]) %>% 
      #   as.matrix %>% 
      #   as.data.frame()
      # 
      # FutureCDFList[[FocalSp]] <- GridDF# %>% dplyr::select(X, Y)
      
      CurrentsGridDF <- 
        readRDS(CurrentFiles[[FocalSp]]) %>% 
        as.matrix %>% 
        as.data.frame() %>% 
        dplyr::select(X, Y)
      
      GridDF <- FutureCDFList[[FocalSp]] %>% 
        dplyr::select(Rep) %>% 
        bind_cols(CurrentsGridDF, .)
      
      SubSharing <- FocalSharing %>% 
        filter(Sp%in%FocalSp|Sp2%in%FocalSp) %>% 
        rename(Value = 3)
      
      GridDF[,paste0("Sharing.", Rep)] <- 
        sum(SubSharing$Value)*GridDF[,Rep]
      
      for(a in 2:length(Species)){
        
        FocalSp = Species[a]
        
        if(a %% 100 == 0){
          
          print(FutureFiles[a])
          
          # (GridDF %>% ggplot(aes(X, Y, fill = Sharing.CLUCurrents)) + geom_tile() + coord_sf()) %>% 
          #   plot
          
        }
        
        SubGretCDF <- FutureCDFList[[FocalSp]]
        
        SubSharing <- FocalSharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp) %>% 
          rename(Value = 3)
        
        SubGretCDF[,paste0("Sharing.", Rep)] <- 
          sum(SubSharing$Value)*SubGretCDF[,Rep]
        
        GridDF[,paste0("Sharing.", Rep)] <-
          GridDF[,paste0("Sharing.", Rep)] +
          SubGretCDF[,paste0("Sharing.", Rep)]
        
        GridDF[,Rep] <-
          GridDF[,Rep] +
          SubGretCDF[,Rep]
        
      }
      
      saveRDS(GridDF, file = paste0("Output Files/GridDFs/", Rep, ".GridDF.rds"))
      
    }
    
    # No bats ####
    
    if(!file.exists(paste0("Output Files/GridDFs/", Rep, ".GridDF_NB.rds"))){
      
      print("No bats!")
      
      NBSpecies <- setdiff(Species, BatSpecies)
      
      FocalSp = NBSpecies[1]
      
      CurrentsGridDF <- 
        readRDS(CurrentFiles[[FocalSp]]) %>% 
        as.matrix %>% 
        as.data.frame() %>% 
        dplyr::select(X, Y)
      
      GridDF <- FutureCDFList[[FocalSp]] %>% 
        dplyr::select(Rep) %>% 
        bind_cols(CurrentsGridDF, .)
      
      SubSharing <- FocalSharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp) %>% 
        rename(Value = 3)
      
      GridDF[,paste0("Sharing.", Rep)] <- 
        sum(SubSharing$Value)*GridDF[,Rep]
      
      for(a in 2:length(NBSpecies)){
        
        FocalSp = Species[a]
        
        if(a %% 100 == 0){
          
          print(FutureFiles[a])
          
          # (GridDF %>% ggplot(aes(X, Y, fill = Sharing.CLUCurrents)) + geom_tile() + coord_sf()) %>% 
          #   plot
          
        }
        
        SubGretCDF <- FutureCDFList[[FocalSp]]
        
        SubSharing <- FocalSharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp) %>% 
          rename(Value = 3)
        
        SubGretCDF[,paste0("Sharing.", Rep)] <- 
          sum(SubSharing$Value)*SubGretCDF[,Rep]
        
        GridDF[,paste0("Sharing.", Rep)] <-
          GridDF[,paste0("Sharing.", Rep)] +
          SubGretCDF[,paste0("Sharing.", Rep)]
        
        GridDF[,Rep] <-
          GridDF[,Rep] +
          SubGretCDF[,Rep]
        
      }
      
      saveRDS(GridDF, file = paste0("Output Files/GridDFs/", Rep, ".GridDF_NB.rds"))
      
    }
  }
  
  rm(FutureCDFList)
  
}

stop()

FileList <- 
  "Iceberg Files/CHELSA/Output Files/IntersectDFs" %>% 
  dir_ls

# c("C", "CD", "CLU", "CLUD") %>% 
#   extract(3:4) %>% 
#   map(function(a){
#     
#     FileList[FileList %>% 
#       str_split("[.]") %>% 
#       map_chr(1) %>% 
#       str_split("/") %>% 
#       map_chr(last) == a] %>% 
#       map(readRDS) %>% 
#       reduce(full_join)
#     
#   })

FullDF <- 
  FileList %>% 
  map(readRDS) %>% 
  bind_cols()

CurrentsGridDF <- 
  readRDS(CurrentFiles[[1]]) %>% 
  as.matrix %>% 
  as.data.frame() %>% 
  dplyr::select(X, Y)



FullDF %>% 
  bind_cols(CurrentsGridDF) %>% 
  gather("Key", "Value", 
         -c(X:Y)) %>% 
  separate("Key", sep = "[.]",
           into = c("Variable",
                    "Pipeline", 
                    "Rep")) %>% 
  mutate(GCM = substr(Rep, 1, 2),
         RCP = substr(Rep, 3, 4),
         Year = substr(Rep, 5, 6)) %>% 
  group_by(Pipeline, 
           RCP, 
           Year) %>% 
  summarise_at("Value", ~mean(.x, na.rm = T))




FullDF %>% 
  dplyr::select()

"Output Files/GridDFs" %>% 
  list.files %>% str_split("[.]") %>% 
  map_chr(1) %>% table

"Output Files/IntersectDFs" %>% 
  list.files %>% str_split("[.]") %>% 
  map_chr(1) %>% 
  table

