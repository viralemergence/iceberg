
# X_Iceberg Organisation ####

library(fs); library(tidyverse)

file_ls()

"Iceberg Revision/Currents/PPM/BinaryMaps" %>% #map(list.files)
  dir_ls(recurse = T, regexp = "[.]tif$") -> 
  
  ToCopyMaxEnt

"Iceberg Revision/Currents/RangeBag/BinaryMaps" %>% #map(list.files)
  dir_ls(recurse = T, regexp = "X0.165") -> 
  
  ToCopyRangeBag

file_copy(ToCopyMaxEnt, "Iceberg Files/Climate1/Iceberg Input Files/MaxEnt/01_Raw/Currents")
file_copy(ToCopyRangeBag, "Iceberg Files/Climate1/Iceberg Input Files/RangeBags/01_Raw/Currents")

ClimateReps <- glue::glue("Climate{1:9}")

ClimateReps[2:9] %>% map(~{
  
  print(.x)
  
  dir_copy("Iceberg Files/Climate1", paste0("Iceberg Files/", .x))
  
})

# Copying futures across #####

"Iceberg Revision" %>% list.files(full.names = T) %>% extract2(1) %>%
  
  paste0("/PPM") %>% list.files(full.names = T) %>% extract2(2) %>% 
  list.files %>%
  substr(1, 2) %>% 
  unique -> 
  
  CoryClimateReps

names(CoryClimateReps) <- ClimateReps

paste0("Iceberg Revision/VC_4_22_20_outputsbSmall/PPM/OtherEnvPred") %>% list.files(full.names = T) ->
  MaxEntFutures

ClimateReps[2:9] %>% 
  
  map(~{
    
    print(.x)
    
    SubReps <- MaxEntFutures[str_detect(MaxEntFutures, paste0("/", CoryClimateReps[.x]))] 
    
    1:length(SubReps) %>% 
      
      lapply(function(a){
        
        print(SubReps[a])
        
        SubReps[a] %>% dir_ls(recurse = T, regexp = "[.]tif$") -> 
          
          SubFiles
        
        SubFiles %>% lapply(function(b) file_copy(b, paste0("Iceberg Files/", .x, 
                                                            "/Iceberg Input Files/MaxEnt/01_Raw/Futures",
                                                            a), overwrite = T))
        
      })
    
  })

paste0("Iceberg Revision/VC_4_22_20_outputsbSmall/RangeBag/OtherEnvPred") %>% 
  list.files(full.names = T) ->
  RangeBagFutures

ClimateReps %>% 
  
  map(~{
    
    print(.x)
    
    SubReps <- RangeBagFutures[str_detect(RangeBagFutures, paste0("/", CoryClimateReps[.x]))] 
    
    1:length(SubReps) %>% 
      
      lapply(function(a){
        
        print(SubReps[a])
        
        SubReps[a] %>% dir_ls(recurse = T, regexp = "X0.165") -> 
          
          SubFiles
        
        SubFiles[str_detect(SubFiles, ".tif$")] -> SubFiles
        
        SubFiles %>% lapply(function(b) file_copy(b, paste0("Iceberg Files/", .x, 
                                                            "/Iceberg Input Files/RangeBags/01_Raw/Futures",
                                                            a), overwrite = T))
        
      })
  })

setwd("~/Albersnet/")

ClimateReps[2:9] %>%
  
  sapply(function(a){
    
    "Iceberg Files/Climate1/Iceberg Input Files/LandUseFutures3.gri" %>%
      file_copy(paste0("Iceberg Files/", a, 
                       "/Iceberg Input Files/LandUseFutures3.gri"), 
                overwrite = T)
    
  })

ClimateReps[2:9] %>%
  
  sapply(function(a){
    
    "Iceberg Files/Climate1/Iceberg Input Files/LandUseFutures3.grd" %>%
      file_copy(paste0("Iceberg Files/", a, 
                       "/Iceberg Input Files/LandUseFutures3.grd"), 
                overwrite = T)
    
  })

# Add something to copy across Futures2 to Futures3 for gf

ClimateReps[which(CoryClimateReps == "gf")]

# Adding something to remove non-.tif files in the input files?


# Removing faulty GretCDFs

ClimateReps[1:9] %>% map(~{
  
  print(.x)
  
  dir_ls(paste0("Iceberg Files/", .x, "/Iceberg Input Files/GretCDF/Futures")) %>%
    sapply(file_remove)
  
})


