
# Flow chart processing figure ####

library(tidyverse); library(raster); library(parallel); library(sf); library(fasterize)

load("~/LargeFiles/MammalStackFullMercator.Rdata")
load("Iceberg Input Files/IUCNBuffers.Rdata")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

FocalSp <- "Cervus_elaphus"
FocalSp <- "Canis_lupus"
FocalSp <- "Vulpes_lagopus"
FocalSp <- "Eidolon_helvum"

FlowStages <- c("01_Raw/Currents", 
                "02_Resampled/Currents",
                "03_MaskedCurrents",
                "04_LandUseClipped/Currents", 
                "06_DispersalBuffersResampled",
                "Final/Currents")

FileLocations <- paste("Iceberg Input Files",
                       Method,
                       FlowStages, sep = "/")

FlowRasters <- list()
FlowRasters[1] <- paste0(FileLocations[1],"/",FocalSp) %>% list.files(full.names = T) %>% raster()
FileLocations[2:6] %>% lapply(function(a) paste0(a,"/",FocalSp,".tif") %>% raster()) ->
  FlowRasters[2:6]

FileLocations[c(2:4,6)] %>% lapply(function(a) paste0(a,"/",FocalSp,".tif") %>% raster()) ->
  FlowRasters[c(2:5)]

names(FlowRasters) <- FlowStages %>% str_split("/") %>% sapply(first)

# Addding IUCN ranges rasters ####

IUCN1 <- MammalStackFull[[FocalSp]]
raster::values(IUCN1)[!is.na(raster::values(IUCN1))] <- 1

sf1 <- sf::st_cast(IUCNBuffers[[FocalSp]], "MULTIPOLYGON")

IUCN2 <- fasterize(sf1, blank)
raster::values(IUCN2)[!is.na(raster::values(IUCN2))] <- 1
  
FinalFlow <- list(FlowRasters[[1]],FlowRasters[[2]],
                  IUCN1, IUCN2,
                  FlowRasters[[3]],FlowRasters[[4]],
                  FlowRasters[[5]],FlowRasters[[6]])

names(FinalFlow) <- c(names(FlowRasters)[1:2], 
                      c("IUCN", "IUCNBuffer"), 
                      names(FlowRasters)[3:6])


FlowDF <- data.frame(
  X = rep(1:ncol(FlowRasters[[2]]), nrow(FlowRasters[[2]])),
  Y = rep(nrow(FlowRasters[[2]]):1, each = ncol(FlowRasters[[2]])),
  Area = AreaValues
)

FlowDF[,names(FinalFlow)[2:length(FinalFlow)]] <- 
  FinalFlow[2:length(FinalFlow)] %>% lapply(function(a) raster::values(a))

# Comparing current and future ####

PredReps <- c("Currents", paste0("Futures", 1:4))

lapply(PredReps[2:5], function(a){
  
  "Iceberg Input Files/MaxEnt/Final/" %>% paste0(a,"/",FocalSp,".tif") %>% raster
  
}) -> FutureRasters

names(FutureRasters) <- PredReps[2:5]

FlowDF[,PredReps[2:5]] <- FutureRasters[PredReps[2:5]] %>% lapply(function(a) raster::values(a))

FlowDF %>% slice(-Sea) %>% 
  gather("Key", "Value", 
         3:ncol(FlowDF)) %>%
  mutate(Key = factor(Key, levels = unique(Key))) -> 
  LongFlowDF

LongFlowDF %>% na.omit %>% ggplot(aes(X, Y, fill = Value)) + 
  coord_fixed() +
  geom_tile() + 
  # lims(x = c(600,900), y = c(450,700)) +
  facet_wrap(~Key)

LongFlowDF %>% ggplot(aes(X, Y, fill = Value)) + 
  coord_fixed() +
  geom_tile()

LongFlowDF %>% na.omit %>% ggplot(aes(X, Y, fill = Key)) + 
  coord_fixed() +
  geom_tile()

