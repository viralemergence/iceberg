
# 8_Display Output #####

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)

PredReps <- c("Currents", paste0("Futures", 1:4))

PipelineReps <- LETTERS[1:4]

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

# Overall trends ###

AllMammaldf %>% summarise_all(mean)
AllMammaldf %>% summarise_at(is.numeric, sum)

#Rasters: #####
#  Current richness
#Future richness A-D x 1-4
#New encounters A-D x 1-4
#Sum sharing A-D x 1-4

CurrentsGridDF <- readRDS("~/Albersnet/Iceberg Output Files/CurrentsGridDF.rds")
FuturesGridDF <- readRDS("~/Albersnet/Iceberg Output Files/FuturesGridDF.rds")

CurrentsGridDF %>% 
  bind_cols(FuturesGridDF[,setdiff(names(FuturesGridDF), names(CurrentsGridDF))]) %>%
  rename(CurrentsAC = Climate, CurrentsBD = ClimateLandUse,
         Sharing.CurrentsAC = Sharing.CurrentsA, Sharing.CurrentsBD = Sharing.CurrentsB) %>%
  rename(Futures1D = Climate.Futures1, Futures2D = Climate.Futures2, 
         Futures3D = Climate.Futures3, Futures4D = Climate.Futures4) %>%
  rename(Futures1C = ClimateLandUse.Futures1, Futures2C = ClimateLandUse.Futures2, 
         Futures3C = ClimateLandUse.Futures3, Futures4C = ClimateLandUse.Futures4) %>%
  rename(Futures1B = BufferClimate.Futures1, Futures2B = BufferClimate.Futures2, 
         Futures3B = BufferClimate.Futures3, Futures4B = BufferClimate.Futures4) %>%
  rename(Futures1A = BufferClimateLandUse.Futures1, Futures2A = BufferClimateLandUse.Futures2, 
         Futures3A = BufferClimateLandUse.Futures3, Futures4A = BufferClimateLandUse.Futures4) %>%
  dplyr::select(-starts_with("LandUse"),X,Y, starts_with("Currents"), 
                ends_with("A"), ends_with("B"), ends_with("C"), ends_with("D")) ->
  
  GretCDF

NewIntersects <- readRDS("~/Albersnet/Iceberg Output Files/NewIntersects.rds")

lapply(PipelineReps, function(a){
  
  lapply(PredReps[2:5], function(b){
    
    NewEncountersList[[a]][[b]]$SharingVar <- NewEncountersList[[a]][[b]][,SharingVars[paste0(b,a)]]
    NewEncountersList[[a]][[b]]$DeltaSharingVar <- NewEncountersList[[a]][[b]][,paste0("Delta",SharingVars[paste0(b,a)])]
    
    NewEncountersList[[a]][[b]] %>% dplyr::summarise(
      Number = n(),
      Mean.Sharing = mean(SharingVar),
      Sum.Sharing = sum(SharingVar),
      New.Sharing = sum(DeltaSharingVar)
    )
    
  }) %>% bind_rows(.id = "PredRep")
  
}) %>% bind_rows(.id = "Pipeline") %>% 
  mutate(Pipeline = PipelineReps[as.numeric(Pipeline)],
         PredRep = PredReps[as.numeric(PredRep)+1])

#Ebola: ####
#  Current richness

EbolaGridDF <- readRDS("~/Albersnet/Iceberg Output Files/EbolaGridList.rds")

#New encounters (Ebola-nonEbola)
#Total new encounters & estimated sharing events A-D x 1-4

EbolaEncounters <- readRDS("Iceberg Output Files/EbolaEncounters.rds")

names(EbolaEncounters) %>% lapply(function(a){
  
  list(
    Number = length(EbolaEncounters[[a]][,SharingVars[a]]),
    Mean.Sharing = mean(EbolaEncounters[[a]][,SharingVars[a]]),
    Sum.Sharing = sum(EbolaEncounters[[a]][,SharingVars[a]]),
    New.Sharing = sum(EbolaEncounters[[a]][,paste0("Delta",SharingVars[a])])
    
  ) %>% return
}) %>% bind_rows()

# Current probability of sharing matrix among _ONLY KNOWN and AFRICAN_ Ebola hosts

AllMammaldf %>% filter(Sp%in%EbolaHosts&Sp2%in%EbolaHosts) %>% summarise_if(is.numeric, mean)
AllMammaldf %>% filter(Sp%in%EbolaHosts&Sp2%in%EbolaHosts) %>% summarise_if(is.numeric, sum)

# New encounters and sharing probability matrix between Ebola-nonEbola, including metadata on host order, in 2.6 A

EbolaEncounters # Includes only one known African Ebola host and one other mammal, and possible data, for all 16 scenarios.

# For example, all hosts in new encounters
EbolaEncounters[[1]] %>% gather("Key", "Value", hOrder.x, hOrder.y) %>% 
  rename(hOrder = Value) %>%
  group_by(hOrder) %>% summarise(N = n()) %>% arrange(N) %>%
  as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
  ggplot(aes(hOrder, N, fill = hOrder)) + geom_col(colour = "black")

# Or if EXCLUDING KNOWN EBOLA HOSTS
EbolaEncounters[[1]] %>% gather("Key", "Value", Sp, Sp2) %>% filter(!Value%in%EbolaHosts) %>%
  mutate(hOrder = ifelse(Key == "Sp", 
                         as.character(hOrder.x), 
                         as.character(hOrder.y))) %>%
  group_by(hOrder) %>% summarise(N = n()) %>% arrange(N) %>%
  as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
  ggplot(aes(hOrder, N, fill = hOrder)) + geom_col(colour = "black")

# Then if looking at CHANGE IN SHARING of these hosts
EbolaEncounters[[1]] %>% gather("Key", "Value", Sp, Sp2) %>% filter(!Value%in%EbolaHosts) %>%
  mutate(hOrder = ifelse(Key == "Sp", 
                         as.character(hOrder.x), 
                         as.character(hOrder.y))) %>%
  group_by(hOrder) %>% summarise(SharingChange = sum(DeltaSharing.Futures1A)) %>% arrange(SharingChange) %>%
  as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
  ggplot(aes(hOrder, SharingChange, fill = hOrder)) + geom_col(colour = "black")

# Then if looking at CHANGE IN SHARING of these hosts
EbolaEncounters[[4]] %>% gather("Key", "Value", Sp, Sp2) %>% filter(!Value%in%EbolaHosts) %>%
  mutate(hOrder = ifelse(Key == "Sp", 
                         as.character(hOrder.x), 
                         as.character(hOrder.y))) %>%
  group_by(hOrder) %>% summarise(SharingChange = sum(DeltaSharing.Futures4A)) %>% arrange(SharingChange) %>%
  as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
  ggplot(aes(hOrder, SharingChange, fill = hOrder)) + geom_col(colour = "black")

#Tallies:
#  New encounters by scenario A-D x 1-4

EbolaEncounters %>% sapply(nrow) %>% matrix(ncol = 4) %>% as.data.frame() ->
  EbolaNETable

rownames(EbolaNETable) <- PredReps[2:5]
colnames(EbolaNETable) <- PipelineReps

# Estimated sharing events by scenario A-D x 1-4
# I think this is above

library(raster)

blank <- raster('~/Albersnet/Iceberg Input Files/UniversalBlank.tif')
rast <- function(x) {
  
  r <- blank 
  
  values(r)[-Sea] <- x 
  
  r
  
}

plot(rast(NewIntersects$Overlap.Futures1A))
plot(rast(NewIntersects$Overlap.Futures1D))

setwd('~/Albersnet/Iceberg Output Files/RasterDrafts')
writeRaster(rast(NewIntersects$Overlap.Futures1A), 'overlap1A.tif')
writeRaster(rast(NewIntersects$Overlap.Futures1D), 'overlap1D.tif')

mean1D <- NewIntersects$OverlapSharing.Futures1D / NewIntersects$Overlap.Futures1D
mean1D[NewIntersects$Overlap.Futures1D<500] <- NA
plot(rast(mean1D))

par(mfrow=c(2,1))
par(mar=c(0,0,0,4))
plot(rast(NewIntersects$Overlap.Futures1A),
     box=FALSE,
     axes=FALSE)
plot(rast(NewIntersects$Overlap.Futures1D),
     box=FALSE,
     axes=FALSE)
plot(rast(NewIntersects$Overlap.Futures4A),
     box=FALSE,
     axes=FALSE)
plot(rast(NewIntersects$Overlap.Futures4D),
     box=FALSE,
     axes=FALSE)



plot(rast(NewIntersects$OverlapSharing.Futures1A),
     box=FALSE,
     axes=FALSE)
plot(rast(NewIntersects$OverlapSharing.Futures1D),
     box=FALSE,
     axes=FALSE)
plot(rast(NewIntersects$OverlapSharing.Futures4A),
     box=FALSE,
     axes=FALSE)
plot(rast(NewIntersects$OverlapSharing.Futures4D),
     box=FALSE,
     axes=FALSE)



plot(rast(GretCDF$CurrentsAC),
     box=FALSE,
     axes=FALSE)
plot(rast(GretCDF$Futures1A)-rast(GretCDF$CurrentsAC),
     box=FALSE,
     axes=FALSE)
plot(rast(GretCDF$CurrentsBD),
     box=FALSE,
     axes=FALSE)
plot(rast(GretCDF$Futures1D)-rast(GretCDF$CurrentsAC),
     box=FALSE,
     axes=FALSE)


#### FIGURE 1 FIRST IMPRESSION
marker = c(color = colorRampPalette(brewer.pal(11,"Spectral"))(100))

par(mfrow=c(4,2))
par(mar=c(0,0,0,4.3))

plot(rast(GretCDF$CurrentsBD),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[10:100])
plot(rast(GretCDF$CurrentsAC),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[10:100])
plot(rast(GretCDF$Futures1D)-rast(GretCDF$CurrentsAC),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)), 
     zlim=c(-275, 275))
plot(rast(GretCDF$Futures1A)-rast(GretCDF$CurrentsAC),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)), 
     zlim=c(-160, 160))
plot(rast(NewIntersects$Overlap.Futures1D),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
plot(rast(NewIntersects$Overlap.Futures1A),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
plot(rast(NewIntersects$OverlapSharing.Futures1D),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
plot(rast(NewIntersects$OverlapSharing.Futures1A),
     box=FALSE,
     axes=FALSE,
     col=rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
