
# Iceberg General Maps ####

GridList <- list()

for(x in 1:length(IcebergAdjList)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/", PredReps[x]))
  
  Velox = T
  
  if(Velox){
    
    VeloxList <- lapply(Files, function(a){
      velox(paste(paste0("Iceberg Input Files/",PredReps[x]), a, sep = '/'))
    })
    
    RasterLista <- lapply(1:length(VeloxList), function(a){
      VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
    })
    
  }
  
  
  if(x==1){
    
    RasterListb <- lapply(1:length(RasterLista), function(a){
      
      if(a %% 500==0) print(a)
      
      raster::resample(RasterLista[[a]], blank, method = 'ngb')
      
    })
    
  } else {
    
    RasterListb <- RasterLista
    
  }
  
  names(RasterListb) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
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
      SubSums[SubSums>0] <- AllMammaldf %>% filter(Sp==sp|Sp2==sp) %>% summarise(mean(paste0(SharingVars[PredReps[x]],2)))
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
  
  GridList[[PredReps[x]]] <- GridDF
  
  remove(RasterListb)
  
}

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

for(i in 1:length(GridList)) GridList[[i]] <- GridList[[i]][-Sea,]

save(GridList, file = "Iceberg Output Files/GridList.Rdata")
