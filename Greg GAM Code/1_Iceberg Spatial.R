# Rscript "Final Iceberg Code/1_Iceberg Spatial.R" ####

library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

CORES = 65

t1 <- Sys.time()

PredReps <- c("Currents", paste0("Futures", 1:4))

RCPs <- c("2.6","4.5","6.0","8.5")

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

Species <- 
  paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$")) %>% 
  sort

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Currents") %>% list.files(full.names = T)) ->
  CurrentFiles

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$")) ->
  names(CurrentFiles)

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Futures") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Futures") %>% list.files(full.names = T)) ->
  FutureFiles

paste0("Iceberg Input Files/","MaxEnt","/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/GretCDF/Futures") %>% list.files() %>% str_remove(".rds$")) ->
  names(FutureFiles)

CurrentCDFList <- FutureCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

PipelineReps <- LETTERS[1:4]

IcebergAdjList <- list()

if(file.exists(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "A",".rds"))&
   file.exists(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "B",".rds"))){
  
  print("Loading RangeAdj Files! Phew")
  
  IcebergAdjList$A <- IcebergAdjList$C <- IcebergAdjList$B <- IcebergAdjList$D <- list()
  
  IcebergAdjList$A$Currents <- IcebergAdjList$C$Currents <- readRDS(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "A",".rds"))
  IcebergAdjList$B$Currents <- IcebergAdjList$D$Currents <- readRDS(paste0("Iceberg Output Files/", "CurrentsRangeAdj", "B",".rds"))
  
} else {
  
  CurrentCDFList <- mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles[a])
    
    readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% as.data.frame() %>% dplyr::select(Climate, ClimateLandUse)
    
  }, mc.preschedule = F, mc.cores = CORES)
  
  object.size(CurrentCDFList)/(10^9)
  
  names(CurrentCDFList) <- Species
  
  # saveRDS(CurrentCDFList, file = paste0("Iceberg Output Files/","CurrentCDFList.rds"))
  
  for(Pipeline in LETTERS[1:4]){
    
    print(Pipeline)
    
    IcebergAdjList[[Pipeline]] <- list()
    
    if(Pipeline == "A"){
      
      CurrentVar <- "ClimateLandUse"
      FuturesVar <- "BufferClimateLandUse"
      
    }
    
    if(Pipeline == "B"){
      
      CurrentVar <- "Climate"
      FuturesVar <- "BufferClimate"
      
    }
    
    if(Pipeline == "C"){
      
      CurrentVar <- "ClimateLandUse"
      FuturesVar <- "ClimateLandUse"
      
    }
    
    if(Pipeline == "D"){
      
      CurrentVar <- "Climate"
      FuturesVar <- "Climate"
      
    }
    
    IcebergAdjList[[Pipeline]] <- list()
    
    if(Pipeline%in%LETTERS[c(1,2)]){
      
      print("Getting values!")
      
      CurrentCDFList %>% map(CurrentVar) %>% 
        map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
        ValueDF
      
      print("Calculating overlap!")
      
      RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
      
      saveRDS(RangeAdj, file = paste0("Iceberg Output Files/", "CurrentsRangeAdj", Pipeline,".rds"))
      
      IcebergAdjList[[Pipeline]]$Currents <- RangeAdj
      
    }else{
      
      IcebergAdjList[[Pipeline]]$Currents <- IcebergAdjList[[which(LETTERS == Pipeline)-2]]$Currents
      
    }
  }
}

# Futures ####
print("Doing the futures!")

FutureCDFList <- mclapply(1:length(Species), function(a){
  
  Sp = Species[a]
  
  print(FutureFiles[a])
  
  readRDS(FutureFiles[[Sp]]) %>% as.matrix %>% as.data.frame() %>% dplyr::select(contains("Futures"), -starts_with("LandUse"))
  
}, mc.preschedule = F, mc.cores = CORES)

names(FutureCDFList) <- Species

object.size(FutureCDFList)/(10^9)

for(Pipeline in LETTERS[1:4]){
  
  FuturesVar <- ifelse(Pipeline == "A", "BufferClimateLandUse", 
                       ifelse(Pipeline == "B", "BufferClimate", 
                              ifelse(Pipeline == "C","ClimateLandUse", "Climate")))
  
  PredReps[2:5] %>% lapply(function(a){
    
    FutureCDFList %>% map(paste0(FuturesVar,".", a)) %>% 
      map(function(b) b*AreaValues[-Sea]) %>% bind_cols() %>% as.data.frame() ->
      ValueDF
    
    RangeAdj <- PairsWisely(ValueDF, Area = T)
    
    saveRDS(RangeAdj, file = paste0("Iceberg Output Files/", a, "RangeAdj", Pipeline,".rds"))
    
    return(RangeAdj)
    
  }) -> IcebergAdjList[[Pipeline]][PredReps[2:5]]
}

IcebergAdjList <- IcebergAdjList[PipelineReps]

saveRDS(IcebergAdjList, file = paste0("Iceberg Output Files/","IcebergAdjList.rds"))

for(Pipeline in LETTERS[1:4]){
  
  print(Pipeline)
  
  CurrentSpecies <- rownames(IcebergAdjList[[Pipeline]][[1]])
  
  for(x in 2:length(IcebergAdjList[[Pipeline]])){
    
    NewAdj <- IcebergAdjList[[Pipeline]][[x]]
    InsertSpecies <- setdiff(CurrentSpecies, rownames(NewAdj))
    
    if(length(InsertSpecies)>0){
      
      NewAdj <- NewAdj %>% data.frame()
      NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
      NewAdj <- NewAdj %>% as.matrix
      
      IcebergAdjList[[Pipeline]][[x]] <- NewAdj[CurrentSpecies, CurrentSpecies]
    }
  }
  
}

saveRDS(IcebergAdjList, file = paste0("Iceberg Output Files/","IcebergAdjList.rds"))

lapply(IcebergAdjList, function(a) sapply(range(a)))
lapply(IcebergAdjList, function(a) sapply(dim(a)))

t2 <- Sys.time() 

t2 - t1

# Trying a comparison 

names(FutureCDFList[[1]])

FocalSp <- "Cervus_elaphus"

Pipeline = "D"
x = 5

NewEncounters <- NewEncountersList[[Pipeline]][[PredReps[x]]]

NewEncounters %>% filter(Sp == FocalSp|Sp2 == FocalSp) -> NewEncounters

NewEncounters %>% filter(Space.Futures1D == max(Space.Futures1D)) %>% pull(Sp2) -> FocalSp2

CurrentCDF <- readRDS(CurrentFiles[[FocalSp]]) %>% as.matrix %>% as.data.frame()

FutureCDF <- readRDS(FutureFiles[[FocalSp]]) %>% as.matrix %>% as.data.frame()

PracticeDF <- bind_cols(CurrentCDF, FutureCDF)

PracticeDF %>% dplyr::select(X,Y,contains(paste0("Climate.Futures"))) %>% 
  mutate(BufferClimate = PracticeDF$BufferClimate) %>%
  mutate(Climate = PracticeDF[,"Climate"]) %>%
  mutate_at(vars(contains(paste0("Climate.Futures"))), function(a){
    
    ifelse(a==1&PracticeDF$Climate==0, "Gain",
           ifelse(a==0&PracticeDF$Climate==1, "Loss",
                  ifelse(a==1&PracticeDF$Climate==1, "Kept", 0)))
    
  }) %>% 
  mutate_all(function(a) ifelse(a == 0, NA, a)) %>% gather("Key", "Value", 3:Climate) %>% na.omit ->
  LongDF

ggplot(LongDF, aes(X, Y, fill = Value)) + geom_tile() + facet_wrap(~Key) + coord_fixed()

FocalSp2 <- "Moschus_chrysogaster"
FocalSp2 <- "Capra_falconeri"

CurrentCDF <- readRDS(CurrentFiles[[FocalSp2]]) %>% as.matrix %>% as.data.frame()

FutureCDF <- readRDS(FutureFiles[[FocalSp2]]) %>% as.matrix %>% as.data.frame()

PracticeDF2 <- bind_cols(CurrentCDF, FutureCDF)

PracticeDF2 %>% dplyr::select(X,Y,contains(paste0("Climate.Futures"))) %>% 
  mutate(BufferClimate = PracticeDF2$BufferClimate) %>%
  mutate(Climate = PracticeDF2[,"Climate"]) %>%
  mutate_at(vars(contains(paste0("Climate.Futures"))), function(a){
    
    ifelse(a==1&PracticeDF2$Climate==0, "Gain",
           ifelse(a==0&PracticeDF2$Climate==1, "Loss",
                  ifelse(a==1&PracticeDF2$Climate==1, "Kept", 0)))
    
  }) %>% 
  mutate_all(function(a) ifelse(a == 0, NA, a)) %>% gather("Key", "Value", 3:Climate) %>% 
  na.omit ->
  LongDF2

ggplot(LongDF2, aes(X, Y, fill = Value)) + geom_tile() + facet_wrap(~Key) + coord_fixed()

LongDF %>% mutate(Sp = FocalSp) %>% bind_rows(LongDF2 %>% mutate(Sp = FocalSp2)) ->
  LongerDF

WideDF <- left_join(LongDF, LongDF2, by = c("X", "Y"), suffix = c(FocalSp, FocalSp2))

WideDF %>% filter(Climate.Futures1D.Cervus_elaphus%in%c("Gain", "Kept")&Climate.Futures1D.Capra_falconeri%in%c("Gain", "Kept"))

LongerDF %>% ggplot(aes(X, Y, fill = Value)) + geom_tile() + geom_tile() + facet_wrap(~Key) + coord_fixed()

PracticeDF %>% dplyr::select(contains(paste0("Climate.Futures"))) %>% 
  mutate(BufferClimate = PracticeDF$BufferClimate) %>%
  mutate(Climate = PracticeDF[,"Climate"]) %>%
  mutate_at(vars(contains(paste0("Climate.Futures"))), function(a){
    
    ifelse(a==1&PracticeDF$Climate==0, "Gain",
           ifelse(a==0&PracticeDF$Climate==1, "Loss",
                  ifelse(a==1&PracticeDF$Climate==1, "Kept", 0)))
    
  }) %>% 
  #mutate_all(function(a) ifelse(a == 0, NA, a)) %>%
  na.omit() %>% 
  apply(2, table) -> Totals

lapply(Totals, function(a) a/last(unlist(Totals))*100)

unlist(Totals)/last(unlist(Totals))
