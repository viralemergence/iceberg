
# 5_Iceberg Making Map Data ####

library(tidyverse); library(Matrix); library(parallel); library(mgcv); 
library(MCMCglmm); library(SpRanger); library(raster)

# Rscript "Final Iceberg Code/5_Iceberg Mapping.R"

Method = "MaxEnt"

PredReps <- c("Currents", paste0("Futures", 1:4))

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

load(paste0("Iceberg Output Files/",Method,"/IcebergAdjList.Rdata"))
load(paste0("Iceberg Output Files/",Method,"/AllMammaldf.Rdata"))
load(paste0("Iceberg Output Files/",Method,"/NewEncounters.Rdata"))

SpaceVars <- paste("Space",PredReps[1:length(IcebergAdjList)], sep = ".")
SharingVars <- paste("Sharing",PredReps[1:length(IcebergAdjList)], sep = ".")

names(SpaceVars) <- names(SharingVars) <- PredReps

# Iceberg General Maps ####

GridList <- NEGridList <- list()

for(x in 1:length(IcebergAdjList)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/",Method,"/Final/", PredReps[x]))
  
  RasterList <- mclapply(Files, function(a){
    raster(paste(paste0("Iceberg Input Files/",Method,"/Final/", PredReps[x]), a, sep = '/'))
  })
  
  names(RasterList) <- Files %>% str_remove(".tif$")
  
  OverlapSums <- OverlapSharingSums <- rep(0, ncol(RasterList[[1]])*nrow(RasterList[[1]]))
  
  print("Getting richness")
  
  for(i in 1:length(RasterList)){
    
    if(i %% 500 == 0) print(i)
    SubSums <- raster::getValues(RasterList[[i]])
    SubSums[is.na(SubSums)] <- 0
    OverlapSums <- OverlapSums + SubSums
    
    sp = names(RasterList)[i]
    if(sp%in%intersect(AllMammaldf$Sp, AllMammaldf$Sp2)){
      SubSums[SubSums>0] <- (AllMammaldf %>% filter(Sp==sp|Sp2==sp))[,SharingVars[PredReps[x]]] %>% sum
      OverlapSharingSums <- OverlapSharingSums + SubSums
    }
  }
  
  
  GridDF <- data.frame(
    Richness = OverlapSums/AreaValues,
    SharingSum = OverlapSharingSums,
    X = rep(1:ncol(RasterList[[1]]), nrow(RasterList[[1]])),
    Y = rep(nrow(RasterList[[1]]):1, each = ncol(RasterList[[1]]))
  ) %>%
    mutate(SharingMean = SharingSum/Richness) %>%
    mutate(SharingMean = ifelse(is.na(SharingMean), 0, SharingMean))
  
  GridDF$PredRep <- PredReps[x]
  
  GridList[[PredReps[x]]] <- GridDF
  
  saveRDS(GridDF, file = paste0("Iceberg Output Files/", Method, "/", PredReps[x],"/",PredReps[x],"GridDF.rds"))
  
  if(x>1){
    
    print(PredReps[x])
    print("New encounters")
    
    FutureRasters <- RasterList
    
    NewEncounters = NewEncountersList[[x-1]]
    
    NewIntersectsManual <- list()
    
    t1 = Sys.time()
    
    i = 1
    
    Rasterstack = RasterList
    PairList = NewEncounters[,c("Sp","Sp2")]
    Species = "All"
    
    t1 <- Sys.time()
    print("Getting the grid values")
    
    if(class(Rasterstack)=="RasterBrick"){
      
      Valuedf <- data.frame(raster::getValues(Rasterstack)) %>% as.matrix
      
    }else{
      
      Valuedf <- lapply(1:length(Rasterstack), function(a){
        
        getValues(Rasterstack[[a]])
        
      }) %>% bind_cols %>% as.data.frame()
      
    }
    
    colnames(Valuedf) <- names(Rasterstack)
    
    Valuedf[is.na(Valuedf)] <- 0
    Valuedf <- Valuedf %>% as.matrix() %>% as("dgCMatrix")
    
    print(paste0("Data frame size = ", nrow(Valuedf)))
    
    if (Species != "All"){
      Valuedf <- Valuedf[, Species]
    }
    
    if (any(Matrix::colSums(Valuedf) == 0)) {
      print("Removing some species with no ranging data :(")
      Valuedf <- Valuedf[, -which(Matrix::colSums(Valuedf) == 0)]
    }
    
    if(length(intersect(c(PairList[,1], PairList[,2]), colnames(Valuedf)))<ncol(Valuedf)){
      print("Removing some species in the PairList but with no rasters :(")
      
      Species <- intersect(c(PairList[,1], PairList[,2]), colnames(Valuedf))
      
      Valuedf <- Valuedf[,Species]
      PairList <- PairList[PairList[,1]%in%Species&PairList[,2]%in%Species,]
      
    }
    
    print(paste0("Data frame size = ", nrow(Valuedf)))
    
    OverlapSums <- OverlapSharingSums <- rep(0, length(AreaValues))
    
    Zeroes <- list()
    
    print("Getting Richness")
    
    for(i in 1:nrow(PairList)){
      
      if(i %% 1000 == 0) print(i)
      
      Sp1 <- PairList[i,1]
      Sp2 <- PairList[i,2]
      
      SubSums <- Valuedf[,Sp1]*Valuedf[,Sp2]
      
      if(sum(SubSums)==0){
        print(paste(i, "Panic!!! Size = 0"))
        Zeroes[[length(Zeroes)+1]] <- i
      }
      
      OverlapSums <- OverlapSums + SubSums
      
      SubSums[SubSums>0] <- NewEncounters[i,SharingVars[x]]
      
      OverlapSharingSums <- OverlapSharingSums + SubSums
      
    }
    
    print("Getting Sharing")
    
    GridDF <- data.frame(
      OverlapSum = OverlapSums/(AreaValues^2),
      SharingSum = OverlapSharingSums,
      X = rep(1:ncol(RasterList[[1]]), nrow(RasterList[[1]])),
      Y = rep(nrow(RasterList[[1]]):1, each = ncol(RasterList[[1]]))
    ) %>%
      mutate(SharingMean = SharingSum/OverlapSum) %>%
      mutate(SharingMean = ifelse(is.na(SharingMean), 0, SharingMean))
    
    GridDF$PredRep <- PredReps[x]
    
    NEGridList[[PredReps[x]]] <- GridDF
    
    saveRDS(GridDF, file = paste0("Iceberg Output Files/", Method, "/", PredReps[x],"/",PredReps[x],"NewEncountersGridDF.rds"))
    
    remove(NewIntersectsManual)
    
  }
  
}

remove(RasterList)

# Removing the sea!

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

for(i in 1:length(GridList)) GridList[[i]] <- GridList[[i]][-Sea,]

save(GridList, file = "Iceberg Output Files/GridList.Rdata")

for(i in 1:length(NEGridList)) NEGridList[[i]] <- NEGridList[[i]][-Sea,]

save(NEGridList, file = "Iceberg Output Files/NEGridList.Rdata")



for(x in 1:length(PredReps)){
  
  GridList[[PredReps[x]]] <- 
    
    readRDS(paste0("Iceberg Output Files/", Method, "/", PredReps[x],"/",PredReps[x],"GridDF.rds"))
  
  if(x>1){
    
    NEGridList[[PredReps[x]]] <- 
      
      readRDS(paste0("Iceberg Output Files/", Method, "/", PredReps[x],"/",PredReps[x],"NewEncountersGridDF.rds"))
    
  }
}


GridList %>% bind_cols %>% #gather("Key", "Value", (starts_with("Richness"))) %>% na.omit() %>%
  ggplot(aes(X, Y, fill = Richness)) + geom_tile() + coord_fixed() 

GridList %>% bind_cols %>% gather("Key", "Value", (starts_with("Richness"))) %>% na.omit() %>%
  ggplot(aes(X, Y, fill = Value)) + geom_tile() + coord_fixed() +
  facet_wrap(~Key)

NEGridList %>% bind_cols %>% #gather("Key", "Value", (starts_with("OverlapSum"))) %>% # na.omit() %>%
  ggplot(aes(X, Y, fill = OverlapSum)) + geom_tile() + coord_fixed()

NEGridList %>% bind_cols %>% gather("Key", "Value", (starts_with("OverlapSum"))) %>% # na.omit() %>%
  ggplot(aes(X, Y, fill = Value)) + geom_tile() + coord_fixed() +
  facet_wrap(~Key)

