
# Rscript "Final Iceberg Code/6c_Range Changes.R" ####

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

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") -> Species ->
  names(FutureFiles)

lapply(1:length(Species), function(i){
  
  Sp <- Species[i]
  
  print(Sp)
  
  CurrentCDF <- readRDS(CurrentFiles[Sp])
  FutureCDF <- readRDS(FutureFiles[Sp]) %>% as.matrix %>% as.data.frame() %>%
    mutate_at(vars(contains("Futures")), function(a) a*(AreaValues[-Sea]))
  
  Clim <- CurrentCDF[,"Climate"]*(AreaValues[-Sea])
  CLU <- CurrentCDF[,"ClimateLandUse"]*(AreaValues[-Sea])
  
  FutureCDF %>% dplyr::select(contains("Climate.Futures")) %>% lapply(function(a){
    
    data.frame(
      
      Loss = sum(Clim[(a==0&Clim>0)], na.rm = T),
      Gain = sum(a[(a>0&Clim==0)], na.rm = T),
      OverallChange = sum(a, na.rm = T)/sum(CLU, na.rm = T))
    
  }) %>% bind_rows(.id = "Rep")  %>% mutate(
    
    PercentLoss = Loss/sum(Clim, na.rm = T),
    PercentGain = Gain/sum(Clim, na.rm = T)
    
  )-> ClimDF
  
  FutureCDF %>% dplyr::select(contains("ClimateLandUse.Futures")) %>% lapply(function(a){
    
    data.frame(
      
      Loss = sum(CLU[(a==0&CLU>0)], na.rm = T),
      Gain = sum(a[(a>0&CLU==0)], na.rm = T),
      OverallChange = sum(a, na.rm = T)/sum(CLU, na.rm = T))
    
  }) %>% bind_rows(.id = "Rep") %>% mutate(
    
    PercentLoss = Loss/sum(CLU, na.rm = T),
    PercentGain = Gain/sum(CLU, na.rm = T)
    
  ) -> CLUDF
  
  ClimDF %>% bind_rows(CLUDF) %>% mutate(Sp = Sp)
  
}) -> ChangeList

saveRDS(ChangeList, file = "Iceberg Output Files/ChangeList.rds")