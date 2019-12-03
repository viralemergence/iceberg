# Iceberg Comparing IUCN and Range Bagging ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

# 1_Iceberg Spatial.R ####

# Doing the raster importing ####

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

CORES = 20

PredReps <- c("Rares", "Commons")

for(x in 1:2){
  
  print(PredReps[x])
  
  Files <- list.files(paste0("Iceberg Input Files/RangeBags/",PredReps[x]))
  
  FileImports <- paste0(paste0(paste0("Iceberg Input Files/RangeBags/",PredReps[x]),"/", Files))
  
  NoFiles <- list()
  
  RasterList <- lapply(FileImports, function(a){
    print(a)
    
    if(length(list.files(a))==0){
      NoFiles[[PredReps[x]]][[a]] <- a
    } else{
      r1 <- raster(paste0(a,"/",list.files(a)[1]))
    }
    
  })
  
  names(RasterList) <- Files
  NoFiles[[PredReps[x]]] <- RasterList[which(sapply(RasterList, class)=="character")]
  RasterList <- RasterList[-which(sapply(RasterList, class)=="character")]
  
  Condition = sapply(RasterList, function(a) any(!is.na(values(a))))
  
  RasterListb <- list()
  
  RasterListb <- mclapply(1:length(RasterList), function(a){
    #for(a in 1:length(RasterList)){
    
    print(a)
    
    #if(Condition[a]){
    
    testraster <- RasterList[[a]]
    testraster <- raster::resample(testraster, blank, method = 'ngb')*AreaRaster
    
    # return(testraster)
    RasterListb[[a]] <- testraster
    #}
  }  , mc.cores = CORES)
  
  names(RasterListb) <- names(RasterList)
  
  RasterListb %>% names %>% mclapply(function(a){
    #for(a in 1:length(RasterListb)){
    
    print(paste0(PredReps[x], a))
    writeRaster(RasterListb[[a]], file = paste0("Iceberg Output Files/ClippedBags/",PredReps[x],"/",a,".tif"))
  }  , mc.cores = CORES)
  
}

remove(RasterList)

CommonRasters <- "Iceberg Output Files/ClippedBags/Commons" %>%
  list.files %>%
  
  lapply(function(a){
    print(a)
    
    r1 <- raster(paste0("Iceberg Output Files/ClippedBags/Commons/",a))
    
  })

names(CommonRasters) <- 
  list.files("Iceberg Output Files/ClippedBags/Commons") %>% 
  str_remove(".tif")

RareRasters <- "Iceberg Output Files/ClippedBags/Rares" %>%
  list.files %>%
  
  lapply(function(a){
    print(a)
    
    r1 <- raster(paste0("Iceberg Output Files/ClippedBags/Rares/",a))
  
  })

names(RareRasters) <- 
  list.files("Iceberg Output Files/ClippedBags/Rares") %>% 
  str_remove(".tif")

RasterList <- combine(CommonRasters, RareRasters)

names(RasterList) <- c(names(CommonRasters),
                       names(RareRasters))

load("~/LargeFiles/MammalStackFullMercator.Rdata")

CommonSp <- intersect(names(RasterList),
                      names(MammalStackFull))

IUCN_RangeBag <- PairsWisening(RasterList, MammalStackFull)

IUCN_RangeBag <- IUCN_RangeBag %>% 
  mutate(Common = as.numeric(Sp%in%names(CommonRasters)))

ggplot(IUCN_RangeBag, aes(POverlap)) + geom_histogram() + 
  facet_wrap(~Common)

RangeBagOverlap <- PairsWisely(RasterList[CommonSp], Area = T)
save(RangeBagOverlap, file = "RangebagOverlap.Rdata")

IUCNOverlap <- PairsWisely(MammalStackFull[CommonSp], Area = T)
save(IUCNOverlap, file = "IUCNOverlap.Rdata")

NonZeroCommonSp <- reduce(list(rownames(RangeBagOverlap),
                               rownames(IUCNOverlap),
                               CommonSp),
                          intersect)

RangeBagOverlap <- RangeBagOverlap[NonZeroCommonSp,NonZeroCommonSp]
IUCNOverlap <- IUCNOverlap[NonZeroCommonSp,NonZeroCommonSp]

Overlapdf <- data.frame(
  
  RangeBag = c(RangeBagOverlap[lower.tri(RangeBagOverlap)]),
  IUCN = c(IUCNOverlap[lower.tri(IUCNOverlap)])
  
)

Overlapdf[,c("Sp", "Sp2")] <- 
  expand.grid(NonZeroCommonSp, NonZeroCommonSp)[which(lower.tri(IUCNOverlap)),]

ggplot(Overlapdf, aes(RangeBag, IUCN)) + geom_point() 

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

save(IcebergAdjList, file = "Iceberg Output Files/IcebergAdjList.Rdata")

lapply(IcebergAdjList, range)
lapply(IcebergAdjList, dim)
