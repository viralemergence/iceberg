
# Master Iceberg Code ####

rm(list = ls())

library(tidyverse); library(fs); library(magrittr)

setwd("~/Albersnet/")

CORES <- 30

# dir_copy("data/", "Iceberg Files/Climate1/data")

glue::glue("Climate{1:9}") ->
  
  ClimateReps

"~/Albersnet/Iceberg Code/Iceberg Greg ENM Code" %>% 
  
  list.files(full.names = T) ->
  
  ENMScripts

ENMScripts[!str_detect(ENMScripts, "00b")] -> ENMScripts

"Iceberg Revision" %>% list.files(full.names = T) %>% extract2(1) %>%
  
  paste0("/PPM") %>% list.files(full.names = T) %>% extract2(2) %>% 
  list.files %>%
  substr(1, 2) %>% 
  unique -> 
  
  CoryClimateReps

names(CoryClimateReps) <- ClimateReps

if(1){
  
  CR <- 1
  
  # for(CR in CR:length(ClimateReps)){
    
    print(ClimateReps[CR])
    
    setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
    
    for(Script in ENMScripts[1]){
      
      print(Script)
      
      source(Script)
      
    }
  }
# }

# GAMM Stuff ####

"~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code" %>% 
  list.files(full.names = T) ->
  
  GAMMScripts

GAMMScripts[!str_detect(GAMMScripts, "Master")] -> GAMMScripts

CR <- 1

for(CR in CR:length(ClimateReps)){
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  for(Script in GAMMScripts[11]){
    
    print(Script)
    
    source(Script)
    
  }
}
