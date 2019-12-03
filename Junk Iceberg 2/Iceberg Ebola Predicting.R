
# Iceberg Ebola Hosts code ####

HP3EbolaHosts <- AssocsBase %>% filter(Virus == "Zaire_ebolavirus")
LauraEbolaHosts <- read.csv("Iceberg Input Files/ZEBOV hosts.csv")

AllSums = PredNetworkList[[1]] 

PredRows <- as.matrix(AllSums)[HP3EbolaHosts$Host,]
PredRowSums <- colSums(PredRows) %>% sort(decreasing = T)

Preddf = data.frame(
  Sp = names(PredRowSums),
  Prob = PredRowSums/length(HP3EbolaHosts$Host)
) %>% mutate(Rank = 1:n(),
             Known = as.numeric(Sp %in%HP3EbolaHosts),
             Predicted = as.numeric(Sp%in%LauraEbolaHosts$bat_species))

Preddf %>% filter(Known==0) %>% SinaGraph("Predicted", "Rank")

EbolaHosts <- union(HP3EbolaHosts$Host, LauraEbolaHosts$bat_species)

EbolaEncounters <- lapply(NewEncountersList, function(a){
  
  a %>% filter(Sp%in%EbolaHosts|Sp2%in%EbolaHosts, !(Sp%in%EbolaHosts&Sp2%in%EbolaHosts))
  
})

names(EbolaEncounters) %>% lapply(function(a){
  
  list(
    Number = length(EbolaEncounters[[a]][,SharingVars[a]]),
    Mean.Sharing = mean(EbolaEncounters[[a]][,SharingVars[a]]),
    Sum.Sharing = sum(EbolaEncounters[[a]][,SharingVars[a]]),
    New.Sharing = sum(EbolaEncounters[[a]][,paste0("Delta",SharingVars[a])])
    
  ) %>% return
})

names(EbolaEncounters) %>% lapply(function(a){
  
  list(
    Number = length(EbolaEncounters[[a]][,paste0(SharingVars[a],2)]),
    Mean.Sharing = mean(EbolaEncounters[[a]][,paste0(SharingVars[a],2)]),
    Sum.Sharing = sum(EbolaEncounters[[a]][,paste0(SharingVars[a],2)]),
    New.Sharing = sum(EbolaEncounters[[a]][,paste0("Delta",paste0(SharingVars[a],2))])
    
  ) %>% return
})

EbolaGridList <- list()

for(x in 1:length(IcebergAdjList)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/Clipped/", PredReps[x]))
  
  Files <- intersect(Files %>% str_replace(" ","_"), paste0(EbolaHosts,".tif"))
  
  RasterListb <- lapply(Files, function(a){
    raster(paste(paste0("Iceberg Input Files/Clipped/",PredReps[x]), a, sep = '/'))
  })
  
  OverlapSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  
  for(i in 1:length(RasterListb)){
    
    if(i %% 500 == 0) print(i)
    SubSums <- raster::getValues(RasterListb[[i]])
    SubSums[is.na(SubSums)] <- 0
    OverlapSums <- OverlapSums + SubSums
    
  }
  
  OverlapSharingSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))
  
  for(i in 1:length(RasterListb)){
    
    sp = names(RasterListb)[i]
    
    if(i %% 500 == 0) print(i)
    if(sp%in%rownames(PredNetworkList[[x]])){
      SubSums <- raster::getValues(RasterListb[[i]])
      SubSums[is.na(SubSums)] <- 0
      SubSums[SubSums>0] <- (AllMammaldf %>% filter(Sp==sp|Sp2==sp))[,paste0(SharingVars[PredReps[x]],2)] %>% mean
      OverlapSharingSums <- OverlapSharingSums + SubSums
    }
  }
  
  GridDF <- data.frame(
    Richness = OverlapSums,
    SharingSum = OverlapSharingSums,
    X = rep(1:ncol(RasterListb[[1]]), nrow(RasterListb[[1]])),
    Y = rep(nrow(RasterListb[[1]]):1, each = ncol(RasterListb[[1]]))
  ) %>%
    mutate(SharingMean = SharingSum/Richness) %>%
    mutate(SharingMean = ifelse(is.na(SharingMean), 0, SharingMean))
  
  GridDF$PredRep <- PredReps[x]
  
  EbolaGridList[[PredReps[x]]] <- GridDF
  
}

EbolaNEGridList <- list()

for(x in 1:length(NewEncountersList)){
  
  print(PredReps[x+1])
  
  Files <- list.files(paste0("Iceberg Input Files/Clipped/", PredReps[x+1]))
  
  NewEncounters <- EbolaEncounters[[x]]
  
  Files <- intersect(Files %>% str_replace(" ","_"), paste0(union(NewEncounters$Sp, NewEncounters$Sp2),".tif"))
  
  RasterLista <- lapply(Files, function(a){
    raster(paste(paste0("Iceberg Input Files/Clipped/",PredReps[x+1]), a, sep = '/'))
  })
  
  FutureRasters <- RasterLista
  
  names(FutureRasters) <- Files %>% 
    str_replace(" ", "_") %>% 
    str_replace("[.]", "_") %>% str_remove(".tif$") 
  
  NewIntersectsManual <- list()
  
  for(i in 1:nrow(NewEncounters)){
    
    if( i %% 1000 == 0) print(i)
    
    NewIntersectsManual[[paste(NewEncounters[i,c("Sp","Sp2")], collapse = ".")]] <- 
      raster::intersect(FutureRasters[[NewEncounters[i,"Sp"]]], 
                        FutureRasters[[NewEncounters[i,"Sp2"]]])
    
  }
  
  names(NewIntersectsManual) <- paste(NewEncounters[,c("Sp","Sp2")], sep = ".")
  
  OverlapSums <- rep(0, ncol(NewIntersectsManual[[1]])*nrow(NewIntersectsManual[[1]]))
  
  for(i in 1:length(NewIntersectsManual)){
    
    if(i %% 1000 == 0) print(i)
    SubSums <- raster::getValues(NewIntersectsManual[[i]])
    SubSums[is.na(SubSums)] <- 0
    if(sum(SubSums)==0) print(i, "Panic!!!")
    OverlapSums <- OverlapSums + SubSums
    
  }
  
  OverlapSharingSums <- rep(0, ncol(NewIntersectsManual[[1]])*nrow(NewIntersectsManual[[1]]))
  
  for(i in 1:length(NewIntersectsManual)){
    
    if( i %% 1000 == 0) print(i)
    SubSums <- raster::getValues(NewIntersectsManual[[i]])
    SubSums[is.na(SubSums)] <- 0
    SubSums[SubSums>0] <- NewEncounters[i,paste0(SharingVars[x+1],2)]
    OverlapSharingSums <- OverlapSharingSums + SubSums
    
  }
  
  GridDF <- data.frame(
    OverlapSum = OverlapSums,
    SharingSum = OverlapSharingSums,
    X = rep(1:ncol(NewIntersectsManual[[1]]), nrow(NewIntersectsManual[[1]])),
    Y = rep(nrow(NewIntersectsManual[[1]]):1, each = ncol(NewIntersectsManual[[1]]))
  ) %>%
    mutate(SharingMean = SharingSum/OverlapSum) %>%
    mutate(SharingMean = ifelse(is.na(SharingMean), 0, SharingMean))
  
  GridDF$PredRep <- PredReps[x+1]
  
  EbolaNEGridList[[PredReps[x+1]]] <- GridDF
  
  remove(NewIntersectsManual)
  
}

for(i in 1:length(EbolaNEGridList)) EbolaNEGridList[[i]] <- EbolaNEGridList[[i]][-Sea,]

save(EbolaNEGridList, file = "Iceberg Output Files/EbolaNEGridList.Rdata")




