

# X_Checking CHELSA ####

library(tidyverse)

FocalSp <- "Catopuma_temminckii"
# FocalSp <- "Acinonyx_jubatus"

Currents <- readRDS(paste0("Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Currents/", FocalSp, ".rds"))

Futures <- readRDS(paste0("Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Futures/", FocalSp, ".rds"))

FullDF <- 
  Currents %>% as.matrix %>% as.data.frame %>% 
  bind_cols(
    
    Futures %>% as.matrix %>% as.data.frame %>% 
      dplyr::select(-c(X, Y, Continent), -matches("Buffer..Climate"))
    
  )

FullDF %>% names

FullDF %>% colSums

FullDF %>% dplyr::select(X, Y, matches("Buffer..Climate")) %>% 
  gather("Key", "Value", -c(X:Y)) %>% 
  ggplot(aes(X, Y, fill = Value)) + 
  geom_tile() + coord_sf() +
  facet_wrap(~Key)

FullDF %>% dplyr::select(X, Y, Climate, matches("BufferClimate.gf")) %>% 
  gather("Key", "Value", -c(X:Y)) %>% 
  ggplot(aes(X, Y, fill = Value)) + 
  geom_tile() + coord_sf() +
  facet_wrap(~Key)

# Importing adj list ####

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr)
library(fs); library(cowplot)

theme_set(theme_cowplot())

setwd(here::here())

IcebergAdjList <- "Iceberg Files/CHELSA/Iceberg Output Files" %>% dir_ls() %>% map(readRDS)

names(IcebergAdjList) <- "Iceberg Files/CHELSA/Iceberg Output Files" %>% list.files %>% 
  str_remove_all(".rds$|RangeAdj") %>% 
  str_remove("^[.]")

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

Witch <- IcebergAdjList[[x]] %>% lower.tri %>% which

PairsDF <- 
  # c("Currents", PredReps) %>% 
  names(IcebergAdjList) %>% 
  map(function(a){
    
    DF <- IcebergAdjList[[a]] %>% reshape2::melt() %>% slice(Witch)
    
    names(DF)[3] <- a
    
    DF %>% return
    
  }) %>% 
  reduce(~full_join(.x, .y)) %>% 
  rename(Sp1 = Var1, Sp2 = Var2)

NonEncounters <- PairsDF %>% filter(Currents == 0)

a <- 26

FocalSp2 <- 
  NonEncounters %>% filter(Sp1 == FocalSp | Sp2 == FocalSp, gf7011 > 0) %>% 
  slice(1) %>% 
  pull(Sp2)

# Next ####

Currents <- readRDS(paste0("Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Currents/", FocalSp2, ".rds"))

Futures <- readRDS(paste0("Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Futures/", FocalSp2, ".rds"))

FullDF2 <- 
  Currents %>% as.matrix %>% as.data.frame %>% 
  bind_cols(
    
    Futures %>% as.matrix %>% as.data.frame %>% 
      dplyr::select(-c(X, Y, Continent), -matches("Buffer..Climate"))
    
  )

CompDF <- 
  FullDF %>% dplyr::select(X, Y, 
                           CurrentSp1 = Climate,
                           FutureSp1 = BufferClimate.gf7011) %>% 
  bind_cols(FullDF2 %>% dplyr::select(FutureSp2 = BufferClimate.gf7011,
                                      CurrentSp2 = Climate))

CompDF %>% 
  gather("Key", "Value", -c(X:Y)) %>% 
  ggplot(aes(X, Y, fill = Value)) + 
  geom_tile() + 
  facet_wrap(~Key) +
  coord_sf()

CompDF %>% 
  mutate(Intersect = case_when(FutureSp1 == 1 & FutureSp2 == 0 ~ FocalSp,
                               FutureSp1 == 1 & FutureSp2 == 1 ~ "Both",
                               FutureSp1 == 0 & FutureSp2 == 1 ~ as.character(FocalSp2),
                               TRUE ~ "Neither")) %>% 
  ggplot(aes(X, Y, fill = Intersect)) + 
  geom_tile() + 
  coord_sf()

