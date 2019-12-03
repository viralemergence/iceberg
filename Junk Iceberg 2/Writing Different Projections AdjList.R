
# Rscript "Writing Different Projections AdjList.R"

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

ProjReps <- c("MollweideProj","MercatorProj","MercatorProjUnc")

IcebergAdjList <- IcebergAdjList2 <- list()

DropSpecies <- c("Lobodon carcinophaga","Leptonychotes weddellii","Enhydra lutris","Histriophoca fasciata") %>% 
  str_replace(" ","_")

PredReps <- c("Currents", paste0("Futures", 1:4))

for(y in 1:length(ProjReps)){
  
  for(x in 1:length(PredReps)){
    
    print(PredReps[x])
    
    Files <- list.files(paste0("Iceberg Output Files/",ProjReps[y],"/",PredReps[x])) %>% setdiff(DropSpecies)
    
    Velox = T
    
    if(Velox){
      
      VeloxList <- lapply(Files, function(a){
        if(which(Files==a) %% 500==0) print(a)
        r1 <- velox(paste(paste0("Iceberg Output Files/",ProjReps[y],"/",PredReps[x]), a, sep = '/'))
      })
      
      RasterLista <- lapply(1:length(VeloxList), function(a){
        
        if(a %% 500==0) print(Files[a])
        
        VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
        
      })
      
    } else{
      RasterLista <- lapply(Files, function(a){
        
        if(a %% 500==0) print(a)
        
        raster(paste(paste0("Iceberg Output Files/",ProjReps[y],"/",PredReps[x]), a, sep = '/'))
        
      })
      
    }
    
    names(RasterLista) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
    
    RasterLista <- RasterLista[setdiff(names(RasterLista), DropSpecies)]
    
    IcebergAdj <- PairsWisely(RasterLista, Area = T)
    
    rownames(IcebergAdj) <- colnames(IcebergAdj) <- rownames(IcebergAdj) %>% str_replace('[.]',"_")
    
    IcebergAdjList[[PredReps[x]]] <- IcebergAdj
    
    IcebergAdj <- PairsWisely(RasterLista, Area = F)
    
    rownames(IcebergAdj) <- colnames(IcebergAdj) <- rownames(IcebergAdj) %>% str_replace('[.]',"_")
    
    IcebergAdjList2[[PredReps[x]]] <- IcebergAdj
    
  }
  
  
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
  
  saveRDS(IcebergAdjList, file = paste0("Iceberg Output Files/IcebergAdjList_",ProjReps[y],".rds"))
  
  CurrentSpecies <- rownames(IcebergAdjList2[[1]])
  
  for(x in 2:length(IcebergAdjList2)){
    
    NewAdj <- IcebergAdjList2[[x]]
    InsertSpecies <- setdiff(CurrentSpecies, rownames(NewAdj))
    
    if(length(InsertSpecies)>0){
      
      NewAdj <- NewAdj %>% data.frame()
      NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
      NewAdj <- NewAdj %>% as.matrix
      
      IcebergAdjList2[[x]] <- NewAdj[CurrentSpecies, CurrentSpecies]
    }
  }
  
  saveRDS(IcebergAdjList2, file = paste0("Iceberg Output Files/IcebergAdjList_",ProjReps[y],"_NoArea",".rds"))
  
}