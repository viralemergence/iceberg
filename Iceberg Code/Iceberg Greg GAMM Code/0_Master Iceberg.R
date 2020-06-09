
# Rscript "Iceberg Code/Iceberg Greg GAMM Code/0_Master Iceberg.R"
# Master Iceberg Code ####

rm(list = ls())

library(tidyverse); library(fs); library(magrittr)

setwd("~/Albersnet/")

CORES <- 60

# dir_copy("data/", "Iceberg Files/Climate1/data")

glue::glue("Climate{1:9}") ->
  
  ClimateReps

"Iceberg Revision" %>% list.files(full.names = T) %>% extract2(1) %>%
  
  paste0("/PPM") %>% list.files(full.names = T) %>% extract2(2) %>%
  list.files %>%
  substr(1, 2) %>%
  unique ->
  
  CoryClimateReps

names(CoryClimateReps) <- ClimateReps

CR <- 1

Spatial <- F

if(Spatial){
  
  for(CR in c(CR:length(ClimateReps))){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    print("Spatial!")
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/1_Iceberg Spatial.R")
    
    if(ClimateReps[CR] == "Climate10"){
      
      ClimateReps[2:9] %>%
        lapply(function(a){
          
          file_copy(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/",
                           "CurrentsRangeAdj", "A",".rds"),
                    
                    paste0("~/Albersnet/Iceberg Files/", a, "/Iceberg Output Files/",
                           "CurrentsRangeAdj", "A",".rds"), 
                    overwrite = T)
          
          file_copy(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/",
                           "CurrentsRangeAdj", "B",".rds"),
                    
                    paste0("~/Albersnet/Iceberg Files/", a, "/Iceberg Output Files/",
                           "CurrentsRangeAdj", "B",".rds"), 
                    overwrite = T)
          
        })
    }
  }
  
  stop()
  
}

CR <- 1

GAMRun <- F

if(GAMRun){
  
  for(CR in (CR:1)){# length(ClimateReps))){
    
    print("Data Import!")
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
    
    print("GAMs!")
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/3_Iceberg GAMs.R")
    
  }
}

print("Prediction!")

PredictRun <- F

if(PredictRun){
  
  CR <- 1
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
  
  CR <- 1
  
  for(CR in (CR:length(ClimateReps))){
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/4_Iceberg Prediction.R")
    
  }
}

print("Mapping!")

MappingRun <- F

if(MappingRun){
  
  CR <- 1
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
  
  #CR <- 4
  
  for(CR in (1)){
    
    # for(CR in (7:6)){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/5_Iceberg Mapping.R")
    
  }
}

SubMappingRun <- F

if(SubMappingRun){
  
  print("SubMapping!")
  
  CR <- 1
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
  
  CR <- 5
  
  for(CR in (5:9)){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/6_Iceberg Submaps.R")
    
  }
}

RasterGAMs <- F

if(RasterGAMs){
  
  print("RasterGAMs!")
  
  CR <- 1
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
  
  CR <- 1
  
  for(CR in (1:9)){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/7_New Encounter GAMs.R")
    
  }
}

Figures <- T

if(Figures){
  
  print("Figures!")
  
  CR <- 1
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
  
  CR <- 1
  
  for(CR in rev(1:9)){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/8_Iceberg Display Outputs.R")
    
  }
}

# Range Changes ####

RangeChanges <- T

if(RangeChanges){
  
  print("RangeChanges!")
  
  CR <- 1
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")
  
  CR <- 1
  
  for(CR in (1:9)){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/6c_Range Change.R")
    
  }
}


