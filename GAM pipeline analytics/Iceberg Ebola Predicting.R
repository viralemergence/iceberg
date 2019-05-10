
# Iceberg Ebola Hosts code ####

HP3EbolaHosts <- AssocsBase %>% filter(Virus == "Zaire_ebolavirus")
LauraEbolaHosts <- read.csv("Iceberg Input Files/ZEBOV hosts.csv")

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

EbolaNEGridList <- list()

for(x in 1:length(NewEncountersList)){
  
  print(PredReps[x+1])
  
  Files <- list.files(paste0("Iceberg Input Files/", PredReps[x+1]))
  
  NewEncounters <- EbolaEncounters[[x]]
  
  Files <- intersect(Files %>% str_replace(" ","_"), paste0(union(NewEncounters$Sp, NewEncounters$Sp2),".tif")) %>%
    str_replace("_"," ")
  
  Velox = T
  
  if(Velox){
    
    VeloxList <- lapply(Files, function(a){
      velox(paste(paste0("Iceberg Input Files/",PredReps[x+1]), a, sep = '/'))
    })
    
    RasterLista <- lapply(1:length(VeloxList), function(a){
      VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
    })
    
  }
  
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
    SubSums[SubSums>0] <- PredNetworkList[[x+1]][NewEncounters[i,"Sp"],NewEncounters[i,"Sp2"]]
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
