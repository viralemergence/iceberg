
# Adding Deltas to the AllMammaldf's because I fucked it originally

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

PipelineReps <- LETTERS[1:4]

CR <- 4

for(CR in c(CR:length(ClimateReps))){
  
  print(ClimateReps[CR])
  
  AllMammaldf <- readRDS(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR], 
                                "/Iceberg Output Files/AllMammaldf.rds"))
  
  AllMammaldf %>% saveRDS(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR], 
                                 "/Iceberg Output Files/AllMammaldfSave.rds"))
  
  PredReps <- c("Currents", paste0("Futures", 1:4))
  
  if(CoryClimateReps[CR] == "gf"){
    
    PredReps <- c("Currents", paste0("Futures", 1:4))[c(1, 2, 4)]
    
  }
  
  SpaceVars <- paste0(paste("Space", PredReps, sep = "."),
                      rep(PipelineReps, each = length(PredReps)))
  
  SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), 
                        rep(PipelineReps, each = length(PredReps)))
  
  names(SpaceVars) <- names(SharingVars) <- paste0(PredReps, rep(PipelineReps, each = length(PredReps)))
  
  for(i in 1:length(PipelineReps)){
    
    AllMammaldf[, paste0("Delta", SpaceVars[2:length(PredReps) + (i-1)*length(PredReps)])] <-
      
      apply(AllMammaldf[, SpaceVars[2:length(PredReps) + (i-1)*length(PredReps)]], 2, function(a){
        a - AllMammaldf[, paste0("Space.Currents", PipelineReps[i])]
        
      })
  }
  
  for(i in 1:length(PipelineReps)){
    
    AllMammaldf[, paste0("Delta", SharingVars[2:length(PredReps) + (i-1)*length(PredReps)])] <-
      
      apply(AllMammaldf[, SharingVars[2:length(PredReps) + (i-1)*length(PredReps)]], 2, function(a){
        a - AllMammaldf[, paste0("Sharing.Currents", PipelineReps[i])]
        
      })
  }
  
  AllMammaldf %>% saveRDS(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR], 
                                 "/Iceberg Output Files/AllMammaldf.rds"))
  
  # Making new encounters ####
  
  NewEncountersList <- 
    
    lapply(PipelineReps, function(b){
      
      l1 <- lapply(PredReps[2:length(PredReps)], function(a){
        
        AllMammaldf[AllMammaldf[,paste0("Space.Currents",b)]==0&
                      AllMammaldf[,paste0("Space.", a, b)]>0,]
        
      })
      
      names(l1) <- PredReps[2:length(PredReps)]
      
      return(l1)
      
    })
  
  names(NewEncountersList) <- PipelineReps
  
  NewEncountersList %>% saveRDS(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR], 
                                       "/Iceberg Output Files/NewEncounters.rds"))
  
  # Making old encounters ####
  
  OldEncountersList <- 
    
    lapply(PipelineReps, function(b){
      
      l1 <- lapply(PredReps[2:length(PredReps)], function(a){
        
        AllMammaldf[AllMammaldf[,paste0("Space.Currents", b)]>0&
                      AllMammaldf[,paste0("Space.", a, b)]==0,]
        
      })
      
      names(l1) <- PredReps[2:length(PredReps)]
      
      return(l1)
      
    })
  
  names(OldEncountersList) <- PipelineReps
  
  OldEncountersList %>% saveRDS(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR], 
                                       "/Iceberg Output Files/OldEncounters.rds"))
  
}
