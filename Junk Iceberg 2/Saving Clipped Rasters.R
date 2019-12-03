# Writing new clipped rasters ###


# Parsing Cory's maps ####

# Rscript "Iceberg Spatial.R"

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

DropSpecies <- c("Lobodon carcinophaga","Leptonychotes weddellii","Enhydra lutris","Histriophoca fasciata") %>% 
  str_replace(" ","_")

IcebergAdjList <- list()

PredReps <- c("Currents", paste0("Futures", 1:4))

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

CORES = 1

for(x in 1:length(PredReps)){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/Clipped/",PredReps[x])) %>% setdiff(DropSpecies)
  
  RasterListc <- lapply(Files, function(a){
    
    r1 = raster(paste(paste0("Iceberg Input Files/Clipped/",PredReps[x]), a, sep = '/'))
    
  })
  
  names(RasterListc) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
  
  RasterListc <- RasterListc[setdiff(names(RasterListc), DropSpecies)]
  
  ContinentClip = T
  
  if(ContinentClip){
    
    NAMES <- setdiff(names(RasterListc), KeepContinents)
    
    for(i in i:length(NAMES)){
      
      print(NAMES[i])
      
      Inhabited <- ContinentWhich[ContinentsInhabited[[NAMES[i]]]] %>% unlist
      
      if(!is.null(Inhabited)){
        
        values(RasterListc[[NAMES[i]]])[-unique(Inhabited)] <- NA
        
      }
      
    }
    
  }
  
  i = 1
  
  for(i in i:length(RasterListc)){
    print(i)
    writeRaster(RasterListc[[i]], file = paste0("Iceberg Input Files/Clipped/",PredReps[x],"/",names(RasterListc)[i],".tif"), overwrite = T)
  }
  
  IcebergAdj <- PairsWisely(RasterListc, Area = F)
  
  saveRDS(IcebergAdj, file = paste0(paste0("Iceberg Output Files/", PredReps[x]),"/IcebergRangeAdj.rds"))
  
  rownames(IcebergAdj) <- colnames(IcebergAdj) <- rownames(IcebergAdj) %>% str_replace('[.]',"_")
  
  IcebergAdjList[[PredReps[x]]] <- IcebergAdj
  
}

remove("RasterLista", "RasterListb", "RasterListc")

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

save(IcebergAdjList, file = "Iceberg Output Files/IcebergAdjList_Mercator_AreaCorrected.Rdata")

NewEncountersQuick <- lapply(IcebergAdjList[2:length(IcebergAdjList)], function(a){
  
  which(c(IcebergAdjList[[1]][lower.tri(IcebergAdjList[[1]])])==0&
          c(a[lower.tri(a)])>0)
  
})

lapply(NewEncountersQuick, length)

lapply(IcebergAdjList, range)

