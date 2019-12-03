# Iceberg Bat-Primate ####

BPEncounters <- lapply(NewEncountersList, function(a){ a %>%
  mutate(BatPrimate = as.numeric(!hOrder.x==hOrder.y) + as.numeric(hOrder.x%in%c("Primates","Chiroptera")&
                                                                     hOrder.y%in%c("Primates","Chiroptera"))) %>% 
  filter(BatPrimate == 2)

  })
  
x = 1

BPNEGridList <- list()

for(x in 1:length(NewEncountersList)){
  
  print(PredReps[x+1])
  
  Files <- list.files(paste0("Iceberg Input Files/Clipped/", PredReps[x+1]))
  
  NewEncounters <- BPEncounters[[x]]
  
  Files <- intersect(Files %>% str_replace(" ","_"), paste0(union(NewEncounters$Sp, NewEncounters$Sp2),".tif")) 
  
  Velox = T
  
  if(Velox){
    
    VeloxList <- lapply(Files, function(a){
      velox(paste(paste0("Iceberg Input Files/Clipped/",PredReps[x+1]), a, sep = '/'))
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
  
  BPNEGridList[[PredReps[x+1]]] <- GridDF
  
  remove(NewIntersectsManual)
  
}

for(i in 1:length(BPNEGridList)) BPNEGridList[[i]] <- BPNEGridList[[i]][-Sea,]

save(BPNEGridList, file = "Iceberg Output Files/BPNEGridList.Rdata")
