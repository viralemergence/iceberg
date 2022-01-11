
# 6_Iceberg Submaps ####

library(tidyverse); library(raster); library(parallel); library(sf); 
library(Matrix); library(magrittr); library(SpRanger); library(cowplot);library(colorspace)

setwd(here::here())
setwd("Iceberg Files/CHELSA")

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Iceberg Ebola Hosts code ####

HP3EbolaHosts <- AssocsBase %>% filter(Virus == "Zaire_ebolavirus")
LauraEbolaHosts <- read.csv("Iceberg Input Files/ZEBOV hosts.csv")

EbolaHosts <- union(HP3EbolaHosts$Host, LauraEbolaHosts$bat_species) %>% 
  intersect(list.files("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% str_remove(".rds$"))

EbolaHosts %>% lapply(function(a){
  
  paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents/",a,".rds") %>% 
    readRDS %>% as.matrix %>% as.data.frame() ->
    Currents
  
  paste0("Iceberg Input Files/GretCDF/Futures/",a,".rds") %>% 
    readRDS %>% as.matrix %>% as.data.frame() ->
    Futures
  
  Currents %>% bind_cols(Futures[,setdiff(names(Futures), names(Currents))])
  
}) -> EbolaGridList

names(EbolaGridList) <- EbolaHosts

EbolaGridList %>% 
  bind_rows(.id = "Host") %>% 
  ggplot(aes(X,Y, fill = ClimateLandUse)) + geom_tile() + 
  scale_fill_continuous_sequential(palette = "Terrain") +
  facet_wrap(~Host)

saveRDS(EbolaGridList, file = "Iceberg Output Files/EbolaGridList.rds")

EbolaHosts %>% setdiff(c("Miniopterus_schreibersii","Pipistrellus_pipistrellus", "Hipposideros_pomona",
                         "Cynopterus_sphinx","Acerodon_jubatus", "Rousettus_leschenaultii")) ->
  
  EbolaHosts

# Continents ####
print("Continents!")

ContinentRaster <- raster("Iceberg Input Files/continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)[-Sea]==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# Doing the new encounters!

#AfricaEbolaEncounters <- readRDS("~/Albersnet/Iceberg Output Files/AfricaEbolaEncounters.rds")##

#1:length(AfricaEbolaEncounters) %>% lapply(function(a){

#  1:4 %>% lapply(function(b){

#    NewEncountersList[[a]][[b]] %>% 
#      filter((Sp%in%AfricaEbolaEncounters[[a]][[b]]|Sp2%in%AfricaEbolaEncounters[[a]][[b]]),
#             (Sp%in%EbolaHosts|Sp2%in%EbolaHosts),
#             !(Sp%in%EbolaHosts&Sp2%in%EbolaHosts),
#             !(Sp%in%setdiff(AfricaEbolaEncounters[[a]][[b]], EbolaHosts)&Sp2%in%setdiff(AfricaEbolaEncounters[[a]][[b]], EbolaHosts)))

#  })

#}) -> EbolaEncounters

EbolaEncounters <- lapply(NewEncountersList, function(a){
  
  a %>% lapply(function(b) b %>% 
                 filter(Sp%in%EbolaHosts|Sp2%in%EbolaHosts, 
                        !(Sp%in%EbolaHosts&Sp2%in%EbolaHosts)))
  
}) %>% unlist(recursive = F)

names(EbolaEncounters) <- paste0(rep(PredReps[2:length(PredReps)], 4), 
                                 rep(PipelineReps, each = length(PredReps[2:length(PredReps)])))

saveRDS(EbolaEncounters, file = "Iceberg Output Files/EbolaEncounters.rds")

names(EbolaEncounters) %>% lapply(function(a){
  
  list(
    Number = length(EbolaEncounters[[a]][,SharingVars[a]]),
    Mean.Sharing = mean(EbolaEncounters[[a]][,SharingVars[a]]),
    Sum.Sharing = sum(EbolaEncounters[[a]][,SharingVars[a]]),
    New.Sharing = sum(EbolaEncounters[[a]][,paste0("Delta",SharingVars[a])])
    
  ) %>% return
})

NEEbolaSpecies <- EbolaEncounters %>% lapply(function(a){
  
  a %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique
  
}) %>% reduce(union) %>% sort

EbolaFutureCDFList <- EbolaCurrentCDFList <- list()

for(i in 1:length(NEEbolaSpecies)){
  
  print(NEEbolaSpecies[i])
  
  EbolaFutureCDFList[[NEEbolaSpecies[i]]] <- 
    readRDS(paste0("Iceberg Input Files/GretCDF/Futures/", NEEbolaSpecies[i],".rds")) %>% 
    as.matrix %>% 
    as.data.frame()
  
  EbolaCurrentCDFList[[NEEbolaSpecies[i]]] <- 
    readRDS(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents/",
                   NEEbolaSpecies[i],".rds")) %>% 
    as.matrix %>% 
    as.data.frame()
  
}

FuturesBase <- c(A = "BufferClimateLandUse",
                 B = "BufferClimate",
                 C = "ClimateLandUse",
                 D = "Climate")

NewIntersectsManual <- EbolaFutureCDFList[[1]] %>% dplyr::select(X, Y)

EbolaNEGridList <- EbolaGridList <- AfricaEbolaHosts <- list()

for(Pipeline in PipelineReps){
  
  print(Pipeline)
  
  EbolaGridList[[Pipeline]] <- EbolaNEGridList[[Pipeline]] <- AfricaEbolaHosts[[Pipeline]] <- list()
  
  for(x in 2:length(PredReps)){
    
    print(PredReps[x])
    
    NewEncounters <- EbolaEncounters[[paste0(PredReps[x],Pipeline)]]
    
    EbolaFutureCDFList %>%
      map(paste0(FuturesBase[Pipeline], ".", PredReps[x])) %>% 
      bind_cols %>% as.data.frame() ->
      ValueDF
    
    ValueDF[-ContinentWhich$Africa,] <- 0
    
    which(colSums(ValueDF) == 0) %>% names %>% setdiff(NEEbolaSpecies, .) ->
      
      AfricaEbolaHosts[[Pipeline]][[PredReps[x]]]
    
    OverlapSums <- OverlapSharingSums <- DeltaOverlapSharing <- rep(0, nrow(ValueDF))
    
    for(y in 1:nrow(NewEncounters)){
      
      if(y %% 10000==0) print(y)
      
      Sp1 <- NewEncounters[y,"Sp"]
      Sp2 <- NewEncounters[y,"Sp2"]
      
      SubSums <- as.numeric(rowSums(ValueDF[,c(Sp1,Sp2)])>1)
      OverlapSums <- OverlapSums + SubSums
      # SubSumList[[paste0(Sp1,".",Sp2)]] <- SubSums
      
      OverlapSharingSums <- OverlapSharingSums +
        
        SubSums*NewEncounters[y, paste0("Sharing.",PredReps[x],Pipeline)]
      
      DeltaOverlapSharing <- DeltaOverlapSharing +
        
        SubSums*NewEncounters[y, paste0("DeltaSharing.",PredReps[x],Pipeline)]
      
    }
    
    NewIntersectsManual[,paste0("Overlap.",PredReps[x],Pipeline)] <-
      OverlapSums
    
    NewIntersectsManual[,paste0("OverlapSharing.",PredReps[x],Pipeline)] <-
      OverlapSharingSums
    
    NewIntersectsManual[,paste0("DeltaOverlapSharing.",PredReps[x],Pipeline)] <-
      DeltaOverlapSharing
    
    saveRDS(NewIntersectsManual, 
            file = "Iceberg Output Files/EbolaNewIntersects.rds")
  }
}

saveRDS(AfricaEbolaHosts, file = "Iceberg Output Files/AfricaEbolaEncounters.rds")

# Iceberg Bat-Primate code ####

BPHosts <- Panth1 %>% filter(hOrder %in% c("Primates", "Chiroptera")) %>% pull(Sp) %>%
  intersect(unlist(AllMammaldf[,c("Sp","Sp2")]))

# Doing the new encounters!

BPEncounters <- lapply(NewEncountersList, function(a){
  
  a %>% lapply(function(b){
    
    b %>% filter(hOrder.x %in%c("Primates", "Chiroptera"), hOrder.y %in%c("Primates", "Chiroptera"), 
                 !hOrder.x == hOrder.y)
    
  })
  
}) %>% unlist(recursive = F)

names(BPEncounters) <- paste0(rep(PredReps[2:length(PredReps)], 4), 
                              rep(PipelineReps, each = length(PredReps[2:length(PredReps)])))

saveRDS(BPEncounters, file = "Iceberg Output Files/BPEncounters.rds")

names(BPEncounters) %>% lapply(function(a){
  
  list(
    Number = length(BPEncounters[[a]][,SharingVars[a]]),
    Mean.Sharing = mean(BPEncounters[[a]][,SharingVars[a]]),
    Sum.Sharing = sum(BPEncounters[[a]][,SharingVars[a]]),
    New.Sharing = sum(BPEncounters[[a]][,paste0("Delta",SharingVars[a])])
    
  ) %>% return
})

NEBPSpecies <- BPEncounters %>% lapply(function(a){
  
  a %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique
  
}) %>% reduce(union) %>% sort

BPFutureCDFList <- BPCurrentCDFList <- list()

for(i in 1:length(NEBPSpecies)){
  
  print(NEBPSpecies[i])
  
  BPFutureCDFList[[NEBPSpecies[i]]] <- 
    readRDS(paste0("Iceberg Input Files/GretCDF/Futures/",NEBPSpecies[i],".rds")) %>% 
    as.matrix %>% 
    as.data.frame()
  
  BPCurrentCDFList[[NEBPSpecies[i]]] <- 
    readRDS(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents/",
                   NEBPSpecies[i],".rds")) %>% 
    as.matrix %>% 
    as.data.frame()
  
}

FuturesBase <- c(A = "BufferClimateLandUse",
                 B = "BufferClimate",
                 C = "ClimateLandUse",
                 D = "Climate")

NewIntersectsManual <- BPFutureCDFList[[1]] %>% dplyr::select(X, Y)

BPGridList <- BPNEGridList <- list()

for(Pipeline in PipelineReps){
  
  print(Pipeline)
  
  BPGridList[[Pipeline]] <- BPNEGridList[[Pipeline]] <- list()
  
  for(x in 2:length(PredReps)){
    
    NewEncounters <- BPEncounters[[paste0(PredReps[x],Pipeline)]]
    
    BPFutureCDFList %>%
      map(paste0(FuturesBase[Pipeline], ".", PredReps[x])) %>% 
      bind_cols %>% as.data.frame() ->
      ValueDF
    
    OverlapSums <- OverlapSharingSums <- DeltaOverlapSharing <- rep(0, nrow(ValueDF))
    
    for(y in 1:nrow(NewEncounters)){
      
      if(y %% 10000==0) print(y)
      
      Sp1 <- NewEncounters[y,"Sp"]
      Sp2 <- NewEncounters[y,"Sp2"]
      
      SubSums <- as.numeric(rowSums(ValueDF[,c(Sp1,Sp2)])>1)
      OverlapSums <- OverlapSums + SubSums
      
      OverlapSharingSums <- OverlapSharingSums +
        
        SubSums*NewEncounters[y, paste0("Sharing.",PredReps[x],Pipeline)]
      
      DeltaOverlapSharing <- DeltaOverlapSharing +
        
        SubSums*NewEncounters[y, paste0("DeltaSharing.",PredReps[x],Pipeline)]
      
    }
    
    NewIntersectsManual[,paste0("Overlap.",PredReps[x],Pipeline)] <-
      OverlapSums
    
    NewIntersectsManual[,paste0("OverlapSharing.",PredReps[x],Pipeline)] <-
      OverlapSharingSums
    
    NewIntersectsManual[,paste0("DeltaOverlapSharing.",PredReps[x],Pipeline)] <-
      DeltaOverlapSharing
    
    saveRDS(NewIntersectsManual, file = "Iceberg Output Files/BPNewIntersects.rds")
    
  }
}
