
# 

library(tidyverse); library(Matrix); library(raster)

PredReps <- c("Currents", paste0("Futures", 1:4))
PipelineReps <- LETTERS[1:4]

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# Blanks

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") -> Species ->
  names(FutureFiles)

ChangeList <- list()

list.files("Iceberg Input Files/GretCDF/Currents") %>% str_remove(".rds$") -> 
  
  Species

i <- 1

map(1:length(Species), function(i){
  
  print(Species[i])
  
  FocalSp <- Species[i]
  
  GretCDF <- 
    
    readRDS(paste0("Iceberg Input Files/GretCDF/Currents/", FocalSp, ".rds")) %>% as.matrix %>% as.data.frame() %>%
    full_join(readRDS(paste0("Iceberg Input Files/GretCDF/Futures/", FocalSp, ".rds")) %>% 
                as.matrix %>% as.data.frame() %>% dplyr::select(-c(Continent, BufferClimate)),
              by = c("X", "Y"), 
              suffix = c(".Currents", ".Futures"))
  
  GretCDF %>% dplyr::select(contains("LandUse"), -starts_with("LandUse"), -c(LandUse, ClimateLandUse)) %>%
    summarise_all(sum) %>% divide_by(sum(GretCDF$ClimateLandUse)) -> CLU
  
  GretCDF %>% dplyr::select(contains("Climate"), -contains("LandUse"), -c(X, Y, Climate, BufferClimate)) %>%
    summarise_all(sum) %>% divide_by(sum(GretCDF$Climate)) -> C
  
  bind_cols(CLU, C) %>% return
  
}) -> RangeChange

RangeChange %>% bind_rows %>% mutate(Species = Species) -> RangeChangeDF

saveRDS(RangeChangeDF, file = "Iceberg Output Files/RangeChangeDF.rds")

RangeChangeDF %>% 
  dplyr::select(ClimateLandUse.Futures1:BufferClimate.Futures4) %>%
  summarise_all(~mean(na.omit(.x[!.x == Inf]))) %>%
  multiply_by(100) ->
  PercentChanges

PercentChanges %>% reshape2::melt() %>% 
  mutate(Pipeline = rep(LETTERS[c(3, 1, 4, 2)], each = 4)) %>%
  mutate(RCP = rep(PredReps[2:5], 4)) %>% 
  rename(RangeChange = value) ->
  PercentChanges

Prev(RangeChangeDF$BufferClimateLandUse.Futures1>1)
Prev(RangeChangeDF$BufferClimateLandUse.Futures4>1)

PercentChanges[PercentChanges$variable == 
                 "BufferClimateLandUse.Futures1", 
               "RangeChange"] - 100

PercentChanges[PercentChanges$variable == 
                 "BufferClimateLandUse.Futures4", 
               "RangeChange"] - 100

RangeChangeDF %>% dplyr::select(-contains("LandUse"), -Species) %>% names ->
  Names
