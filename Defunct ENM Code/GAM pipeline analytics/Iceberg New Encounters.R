
# Iceberg New Overlaps ####

# Making the raster overlaps ####

NEGridList <- list()

for(x in 2:length(NewEncountersList)){
  
  print(PredReps[x+1])
  
  Files <- list.files(paste0("Iceberg Input Files/", PredReps[x+1]))
  
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
  
  NewEncounters = NewEncountersList[[x]]
  
  NewIntersectsManual <- list()
  
  for(i in 1:nrow(NewEncounters)){
    
    if( i %% 1000 == 0) print(i)
    
    NewIntersectsManual[[paste(NewEncounters[i,c("Sp","Sp2")], collapse = ".")]] <- 
      raster::intersect(FutureRasters[[NewEncounters[i,"Sp"]]], 
                        FutureRasters[[NewEncounters[i,"Sp2"]]])
    
  }
  
  names(NewIntersectsManual) <- paste(NewEncounters[1:length(NewIntersectsManual),c("Sp","Sp2")], sep = ".")
  
  OverlapSums <- rep(0, ncol(NewIntersectsManual[[1]])*nrow(NewIntersectsManual[[1]]))
  
  Zeroes <- list()
  
  for(i in 1:length(NewIntersectsManual)){
    
    if(i %% 1000 == 0) print(i)
    SubSums <- raster::getValues(NewIntersectsManual[[i]])
    SubSums[is.na(SubSums)] <- 0
    if(sum(SubSums)==0){
      print(paste(i, "Panic!!! Size = 0"))
      Zeroes[[length(Zeroes)+1]] <- i
    }
    OverlapSums <- OverlapSums + SubSums
    
  }
  
  OverlapSharingSums <- rep(0, ncol(NewIntersectsManual[[1]])*nrow(NewIntersectsManual[[1]]))
  
  for(i in 1:length(NewIntersectsManual)){
    
    if( i %% 1000 == 0) print(i)
    SubSums <- raster::getValues(NewIntersectsManual[[i]])
    SubSums[is.na(SubSums)] <- 0
    SubSums[SubSums>0] <- NewEncounters[i,SharingVars[x+1]]
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
  
  NEGridList[[PredReps[x+1]]] <- GridDF
  
  remove(NewIntersectsManual)
  
}

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

for(i in 1:length(NEGridList)) NEGridList[[i]] <- NEGridList[[i]][-Sea,]

save(NEGridList, file = "Iceberg Output Files/NEGridList.Rdata")

Compdf <- NEGridList %>% bind_rows() %>% 
  left_join(GridList %>% bind_rows, by = c("X", "Y", "PredRep")) %>%
  mutate(Overlapcorrect = OverlapSum/Richness^2)

NEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = OverlapSum)) + 
  coord_fixed() + 
  theme_void() +
  facet_wrap(~PredRep)

NEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = SharingSum)) + 
  coord_fixed() +
  theme_void() +
  facet_wrap(~PredRep)

NEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = SharingMean)) + 
  coord_fixed() +
  theme_void() +
  facet_wrap(~PredRep)

Compdf %>%  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = SharingMean.x)) + 
  coord_fixed() +
  theme_void() +
  facet_wrap(~PredRep)

Compdf %>% mutate(Overlapcorrect = ifelse(Overlapcorrect>100, 100, Overlapcorrect)) %>%  
  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = Overlapcorrect)) + 
  coord_fixed() +
  theme_void() +
  facet_wrap(~PredRep)
