
# X_Chelsa Analysis ####

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

EncounterList <-
  c(26, 70, 85)[-3] %>% map(function(a){
    
    New11 <- NonEncounters[NonEncounters[,paste0("gf", a, 11)] > 0, c("Sp1", "Sp2")] %>% 
      dplyr::select(Sp1, Sp2)
    
    New41 <- NonEncounters[NonEncounters[,paste0("gf", a, 41)] > 0, c("Sp1", "Sp2")] %>% 
      dplyr::select(Sp1, Sp2)
    
    # New41 %<>% full_join(New11, by = c("Sp1", "Sp2"))
    
    New71 <- NonEncounters[NonEncounters[,paste0("gf", a, 71)] > 0, c("Sp1", "Sp2")] %>% 
      dplyr::select(Sp1, Sp2)
    
    # New71 %<>% full_join(New41, by = c("Sp1", "Sp2"))
    
    list(New11 = New11,
         New41 = New41,
         New71 = New71)
    
  })

EncounterList %>% saveRDS("Iceberg Output Files/EncounterList.rds")

NEncounterDF <- 
  EncounterList %>% 
  map(function(a){
    
    data.frame(
      
      New91 = 0,
      New11 = a$New11 %>% nrow,
      New41 = a$New11 %>% bind_rows(a$New41) %>% unique %>% 
        nrow,
      New71 = a$New11 %>% bind_rows(a$New41) %>% bind_rows(a$New71) %>% unique %>% 
        nrow
      
    )
  })

LongNEncounterDF <- 
  NEncounterDF %>% 
  bind_rows(.id = "RCP") %>% gather("Year", "Encounters", -RCP) %>% 
  mutate_at("Year", ~substr(.x, 4, 5) %>% as.numeric %>% add(2000)) %>% 
  mutate_at("Year", ~ifelse(.x > 2090, .x - 100, .x))

LongNEncounterDF %>% ggplot(aes(Year, Encounters)) + 
  geom_line(aes(colour = RCP)) + 
  geom_point() + 
  lims(y = c(0, NA), x = c(1990, NA))

